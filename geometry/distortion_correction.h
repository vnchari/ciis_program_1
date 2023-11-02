//
// Created by Vivek Chari on 11/1/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H
#define PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H

#include <iostream>


template<typename T>
struct DISTORTIONCOEFFICIENTS{
    Eigen::Matrix<T, 216, 3> coeffs;
    Eigen::Vector<T, 3> min_coeffs;
    Eigen::Vector<T, 3> max_coeffs;
};

template<typename T>
T fifth_degree_bernstein_polynomial(int k, T v) {
  T choose_vals[6] = {1, 5, 10, 10, 5, 1};
  return choose_vals[k] * pow(1 - v, 5 - k) * pow(v, k);
}

template<typename T>
Eigen::Matrix<T, -1, 3> scale_to_box(Eigen::Matrix<T, -1, 3> prescale, DISTORTIONCOEFFICIENTS<T> precomputed) {
  Eigen::Matrix<T, -1, 3> postscale(prescale.rows(), 3);
  postscale.col(0) = ((prescale.col(0).array() - precomputed.min_coeffs(0)) / (precomputed.max_coeffs(0) - precomputed.min_coeffs(0))).matrix();
  postscale.col(1) = ((prescale.col(1).array() - precomputed.min_coeffs(1)) / (precomputed.max_coeffs(1) - precomputed.min_coeffs(1))).matrix();
  postscale.col(2) = ((prescale.col(2).array() - precomputed.min_coeffs(2)) / (precomputed.max_coeffs(2) - precomputed.min_coeffs(2))).matrix();
  return postscale;
}

template<typename T>
struct DISTORTIONCOEFFICIENTS<T> compute_distortion_coeffs(Eigen::Matrix<T, Eigen::Dynamic, 3> measured,
        Eigen::Matrix<T, Eigen::Dynamic, 3> ground_truth) {
  DISTORTIONCOEFFICIENTS<T> c;
  c.min_coeffs = measured.colwise().minCoeff();
  c.max_coeffs = measured.colwise().maxCoeff();
  Eigen::Matrix<T, Eigen::Dynamic, 3> u = scale_to_box(measured, c);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(measured.rows(), 216);
  // compiler will flatten this
  for(int point = 0; point < measured.rows(); point++){
    for(int i = 0; i <= 5; i++) {
      for(int j = 0; j <= 5; j++) {
        for(int k = 0; k <= 5; k++) {
          A(point, 36 * i + 6 * j + k) = fifth_degree_bernstein_polynomial(i, u(point, 0))
                  * fifth_degree_bernstein_polynomial(j, u(point, 1))
                  * fifth_degree_bernstein_polynomial(k, u(point, 2));
        }
      }
    }
  }
  c.coeffs = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ground_truth);
  return c;
}

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 3> apply_distortion_correction(DISTORTIONCOEFFICIENTS<T> coeffs, Eigen::Matrix<T, Eigen::Dynamic, 3> points) {
  Eigen::Matrix<T, Eigen::Dynamic, 3> output = Eigen::Matrix<T, Eigen::Dynamic, 3>::Zero(points.rows(), 3);
  Eigen::Matrix<T, Eigen::Dynamic, 3> u = scale_to_box(points, coeffs);
  for(int point = 0; point < points.rows(); point++){
    for(int i = 0; i <= 5; i++) {
      for(int j = 0; j <= 5; j++) {
        for(int k = 0; k <= 5; k++) {
          output.row(point) += coeffs.coeffs.row(36 * i + 6 * j + k)
                  * fifth_degree_bernstein_polynomial(i, u(point, 0))
                  * fifth_degree_bernstein_polynomial(j, u(point, 1))
                  * fifth_degree_bernstein_polynomial(k, u(point, 2));
        }
      }
    }
  }
  return output;
}


template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 3> pack_points(std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> points, int num_frames, int num_markers) {
  Eigen::Matrix<T, -1, 3> packed = Eigen::Matrix<T, -1, 3>(num_frames * num_markers, 3);
  for(int i = 0; i < num_frames; i++){
    packed.template block(i * num_markers, 0, num_markers, 3) = points[i].transpose();
  }
  return packed;
}

template<typename T>
std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> unpack_points(Eigen::Matrix<T, Eigen::Dynamic, 3> packed, int num_frames, int num_markers) {
  std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> unpacked(num_frames);
  for(int i = 0; i < num_frames; i++){
    unpacked[i] = packed.template block(i * num_markers, 0, num_markers, 3).transpose();
  }
  return unpacked;
}


template<typename T>
std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> apply_distortion_correction(DISTORTIONCOEFFICIENTS<T> coeffs, std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> points) {
  int n_frames = points.size(); int num_markers = points[0].cols();
  return unpack_points<T>(apply_distortion_correction<T>(coeffs, pack_points<T>(points, n_frames, num_markers)), n_frames, num_markers);
}


template<typename T>
struct DISTORTIONCOEFFICIENTS<T> compute_distortion_coeffs(std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> measured,
                                                           std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>> ground_truth) {
  int n_frames = measured.size(); int num_markers = measured[0].cols();
  return compute_distortion_coeffs<T>(pack_points<T>(measured, n_frames, num_markers),
                                   pack_points<T>(ground_truth, n_frames, num_markers));
}

#endif //PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H
