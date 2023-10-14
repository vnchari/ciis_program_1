//
// Created by Vivek Chari on 10/11/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_REGISTER_H
#define PROGRAMMING_ASSIGMENT_ONE_REGISTER_H

#include "Eigen/Dense"
#include "Eigen/SVD"

namespace Registration {

    uint64_t COHERENTPOINTDRIFT = 0;
    uint64_t PROCRUSTES = 1;

template<typename T>
class CoherentPointDrift {
    long D = 3; //dimension
    T w = 0.0;
    Eigen::Matrix<T, -1, 3, Eigen::RowMajor> M;
    Eigen::Matrix<T, -1, 3, Eigen::RowMajor> N;
    Eigen::Matrix<T, -1, -1, Eigen::ColMajor> prob;
    Eigen::JacobiSVD<Eigen::Matrix<T, -1, -1, Eigen::RowMajor>> EIGENSVD;
    T sigma_square = 0;
    T tol;
    long iter;
    double size_N, size_M; //number of points
public:
    Eigen::Matrix<T, 3, 3> B; //rotation
    Eigen::Vector<T, 3> t; //translation
    CoherentPointDrift(Eigen::Matrix<T, -1, 3, Eigen::RowMajor> M,
                       Eigen::Matrix<T, -1, 3, Eigen::RowMajor> N,
                       long niter = 250, T term_tol = 0.001);

    void register_points_cpd();

    void expectation();

    void maximization();
};

template<typename T>
void CoherentPointDrift<T>::register_points_cpd() {
  for (int i = 0; i < iter; i++) {
    auto v = sigma_square;
    expectation();
    maximization();
    if (abs(v - sigma_square) < tol) {
      break;
    }
  }
}

template<typename T>
void CoherentPointDrift<T>::maximization() {
  auto N_p = prob.sum();
  auto colwise_prob = (prob.colwise().sum()).transpose(); //M
  auto rowwise_prob = prob.rowwise().sum() * (1.0f / N_p); //N
  auto mu_m_hat = (M.transpose() * colwise_prob * (1.0f / N_p));
  auto mu_n_hat = (N.transpose() * rowwise_prob);
  auto Mhat = M.rowwise() - mu_m_hat.transpose();
  auto Nhat = N.rowwise() - mu_n_hat.transpose();

  auto A = Mhat.transpose() * prob.transpose() * Nhat;
  EIGENSVD.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  auto U = EIGENSVD.matrixU();
  auto V_T = EIGENSVD.matrixV().transpose();
  Eigen::Matrix<T, 3, 3, Eigen::RowMajor> C = Eigen::Matrix<T, 3, 3, Eigen::RowMajor>::Identity();
  C(2, 2) = (U * V_T).determinant();
  B = U * C * V_T;

  t = mu_m_hat - B * mu_n_hat;
  sigma_square = (1.0f / (N_p * D)) *
                 ((Mhat.transpose() * colwise_prob.asDiagonal() * Mhat).trace() - (A.transpose() * B).trace());
}


template<typename T>
void CoherentPointDrift<T>::expectation() {
  T inv_sigma_square = -0.5 / sigma_square;
  for (int m = 0; m < size_M; m++) {
    T denom = (size_M / size_N) * (w / (1 - w));//* std::pow(2.0 * PI * sigma_square, D/2.0);
    const auto x = M.row(m);
    for (int n = 0; n < size_N; n++) {
      const auto y = N.row(n);
      T probmass = exp(inv_sigma_square * (x.transpose() - (B * y.transpose() + t)).squaredNorm());
      prob(n, m) = probmass;
      denom += probmass;
    }
    prob.col(m) /= denom;
  }
}

template<typename T>
CoherentPointDrift<T>::CoherentPointDrift(Eigen::Matrix<T, -1, 3, Eigen::RowMajor> M_arg,
                                          Eigen::Matrix<T, -1, 3, Eigen::RowMajor> N_arg,
                                          long niter,
                                          T term_tol) {
  tol = term_tol;
  iter = niter;
  M = M_arg;
  N = N_arg;
  size_M = M.rows();
  size_N = N.rows();
  for (int n = 0; n < size_N; n++)
    for (int m = 0; m < size_M; m++)
      sigma_square += (N.row(n) - M.row(m)).squaredNorm();
  sigma_square /= (T) (D * size_M * size_N);

  prob = Eigen::Matrix<T, -1, -1, Eigen::ColMajor>::Zero(size_N, size_M);
  B = Eigen::Matrix<T, 3, 3>::Identity();
  t = Eigen::Vector<T, 3>::Zero();
}

template<typename T>
class Procrustes {
    Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> EIGENSVD;
public:
    Eigen::Matrix<T, 3, 3> B;
    Eigen::Vector<T, 3> t;

    Procrustes(Eigen::Matrix<T, -1, 3> M, Eigen::Matrix<T, -1, 3> N);
};

template<typename T>
Procrustes<T>::Procrustes(const Eigen::Matrix<T, -1, 3> M,
                          const Eigen::Matrix<T, -1, 3> N) {
  // N = FM = RM + t
  if (M.rows() != N.rows())
    throw std::invalid_argument("M and N must have the same number of rows.");
  auto num_points = M.rows();

  // compute the centroid of both point sets.
  auto center_M = M.colwise().sum() / num_points;
  auto center_N = N.colwise().sum() / num_points;

  // centering the point sets so that centroid lies at origin.
  auto centered_M = M.rowwise() - center_M;
  auto centered_N = N.rowwise() - center_N;

  //
  EIGENSVD.compute(centered_M.transpose() * centered_N, Eigen::ComputeFullU | Eigen::ComputeFullV);
  auto U_T = EIGENSVD.matrixU().transpose();
  auto V = EIGENSVD.matrixV();
  Eigen::Matrix<T, 3, 3> C = Eigen::Matrix<T, 3, 3>::Identity();
  C(2, 2) = (V * U_T).determinant();
  B = V * C * U_T;
  t = center_N - center_M * B.transpose();
}

}

#endif //PROGRAMMING_ASSIGMENT_ONE_REGISTER_H
