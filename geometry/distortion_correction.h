//
// Created by Vivek Chari on 11/1/23.
//
// The following file contains methods related to bernstein polynomial distortion correction.
// It includes methods to both compute the coefficients, and apply a distortion correction.
// One can pass either a list of frames or a "packed" matrix to the aforementioned methods.
// Additionally, this file contains  unit tests at the component level for each of the
// methods and helper methods implement below.

#ifndef PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H
#define PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H

#include <map>
#include <iostream>
#include <Eigen/StdVector>
#include <Eigen/SVD>

std::map<int, std::vector<int>> binomial_coeffs;
template<typename T>
using VECTOROFXD3MATRICES = std::vector<Eigen::Matrix<T, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<T, -1, 3>>>;
template<typename T>
using VECTOROF3XDMATRICES = std::vector<Eigen::Matrix<T, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, -1>>>;

//helper data structure for storing values needed for distortion correction
template<int degree, typename T>
struct DISTORTIONCOEFFICIENTS{
    Eigen::Matrix<T, (int) pow(degree + 1, 3), 3> coeffs;
    Eigen::Vector<T, 3> min_coeffs;
    Eigen::Vector<T, 3> max_coeffs;
};

template<size_t degree, typename T>
struct DISTORTIONCOEFFICIENTS<degree, T> compute_distortion_coeffs(Eigen::Matrix<T, Eigen::Dynamic, 3> measured,
        Eigen::Matrix<T, Eigen::Dynamic, 3> ground_truth) {
  DISTORTIONCOEFFICIENTS<degree, T> output;
  //lambda function that compute the coefficient associated with the bernstein polynomial at a point in space.
  auto F = [](size_t index, Eigen::Vector3<T> vec) {
    return nth_degree_bernstein_polynomial<degree, T>((index / (int) pow(degree + 1, 2)) % (degree + 1), vec(0)) *
            nth_degree_bernstein_polynomial<degree, T>((index / (int) pow(degree + 1, 1)) % (degree + 1), vec(1)) *
                    nth_degree_bernstein_polynomial<degree, T>((index) % (degree + 1), vec(2));
  };
  //compute the min/max coefficient for each column. Save the coeffs for use when applying the correction.
  output.min_coeffs = measured.colwise().minCoeff(); output.max_coeffs = measured.colwise().maxCoeff();
  Eigen::Matrix<T, Eigen::Dynamic, 3> scaled_points = scale_to_unit_box<T>(measured, output.min_coeffs, output.max_coeffs);
  //create a (num_frames*num_markers) x ( (degree + 1)^3) matrix. The number of columns is the number of ways to
  // form a sequence of 3 numbers whose elements are of 0-degree inclusive
  size_t num_sequences = pow(degree + 1, 3);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(scaled_points.rows(), num_sequences);
  //fill the design matrix.
  for(int i = 0; i < measured.rows(); i++)
    for(size_t j = 0; j < num_sequences; j++)
      //accessing this way improves cache locality.
      A(i, j) = F(j, scaled_points.row(i));
  //compute M-P pseudoinverse to solve the lstsq problem
  output.coeffs = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ground_truth);
  return output;
}

template<size_t degree, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 3> apply_distortion_correction(DISTORTIONCOEFFICIENTS<degree, T> coeffs, Eigen::Matrix<T, Eigen::Dynamic, 3> points) {
  //lambda function that computes the coefficient associated with the bernstein polynomial at a point in space.
  auto F = [](size_t index, Eigen::Vector3<T> vec) {
      return nth_degree_bernstein_polynomial<degree, T>((index / (int) pow(degree + 1, 2)) % (degree + 1) , vec(0)) *
             nth_degree_bernstein_polynomial<degree, T>((index / (int) pow(degree + 1, 1)) % (degree + 1), vec(1)) *
             nth_degree_bernstein_polynomial<degree, T>((index) % (degree + 1), vec(2));
  };

  Eigen::Matrix<T, Eigen::Dynamic, 3> scaled_points = scale_to_unit_box<T>(points, coeffs.min_coeffs, coeffs.max_coeffs);
  //create a (num_frames*num_markers) x ( (degree + 1)^3) matrix. The number of columns is the number of ways to
  // form a sequence of 3 numbers whose elements are of 0-degree inclusive
  size_t num_sequences = pow(degree + 1, 3);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(scaled_points.rows(), num_sequences);
  //fill the design matrix.
  for(int i = 0; i < points.rows(); i++)
    for(size_t j = 0; j < num_sequences; j++)
      //accessing this way improves cache locality.
      A(i, j) = F(j, scaled_points.row(i));
  // right action of the coefficients on the design matrix computes the ground truth values.
  return A * coeffs.coeffs;
}


template<typename T>
void unit_test_compute_distortion_coeffs() {
  Eigen::Matrix<T, 4, 3> mat = Eigen::Matrix<T, 4, 3>::Random();
  Eigen::Matrix<T, 4, 3> distorted = mat; // Bernstein(degree=1, index=0)
  std::cout << compute_distortion_coeffs<1, T>(distorted, mat).coeffs  << std::endl << std::endl;

//  assert(is_zero);
//  std::cout << apply_distortion_correction<5, T>(compute_distortion_coeffs<5, T>(mat, mat), mat) - mat << std::endl << std::endl;
//  std::cout << apply_distortion_correction_2<T>(compute_distortion_coeffs_2<T>(mat, mat2), mat) - mat2 << std::endl << std::endl;
////  std::cout << apply_distortion_correction<T>(compute_distortion_coeffs_2<T>(mat, mat), mat) - mat << std::endl << std::endl;
//  Eigen::Matrix<T, 5, 3> mat2 = Eigen::Matrix<T, 5, 3>::Random();

}

//overloaded apply_distortion_correction to allow for passing a vector of frames, packing the frames, and then unpacking
//just for convenience.
template<int degree, typename T>
VECTOROF3XDMATRICES<T> apply_distortion_correction(DISTORTIONCOEFFICIENTS<degree, T> coeffs, VECTOROF3XDMATRICES<T> points) {
  int n_frames = points.size(); int num_markers = points[0].cols();
  return unpack_points<T>(apply_distortion_correction<degree, T>(coeffs, pack_points<T>(points, n_frames, num_markers)), n_frames, num_markers);
}

//overloaded compute_distortion_coeffs to allow for passing a vector of frames, packing the frames, and then unpacking.
//just for convenience.
template<int degree, typename T>
struct DISTORTIONCOEFFICIENTS<degree, T> compute_distortion_coeffs(VECTOROF3XDMATRICES<T> measured,
                                                           VECTOROF3XDMATRICES<T> ground_truth) {
  int n_frames = measured.size(); int num_markers = measured[0].cols();
  return compute_distortion_coeffs<degree, T>(pack_points<T>(measured, n_frames, num_markers),
                                   pack_points<T>(ground_truth, n_frames, num_markers));
}



//takes a single matrix of size (num_frames * num_markers) x 3 and unpacks the matrix into a
// vector of num_frames different 3 x num_markers matrices
template<typename T>
VECTOROF3XDMATRICES<T> unpack_points(Eigen::Matrix<T, Eigen::Dynamic, 3> packed, int num_frames, int num_markers) {
  VECTOROF3XDMATRICES<T> unpacked;
  for(int i = 0; i < num_frames; i++){
    unpacked.push_back(packed.template block(i * num_markers, 0, num_markers, 3).transpose());
  }
  return unpacked;
}

template<typename T>
void unit_test_unpack_points() {
  int num_frames = 4; int num_markers = 2;
  Eigen::Matrix<T, 8, 3> rand = Eigen::Matrix<T, 8, 3>::Random();
  VECTOROF3XDMATRICES<T> result = unpack_points<T>(rand, num_frames, num_markers);
  assert(result[0].transpose().row(0) == rand.row(0));
  assert(result[0].transpose().row(1) == rand.row(1));
  assert(result[1].transpose().row(0) == rand.row(2));
  assert(result[1].transpose().row(1) == rand.row(3));
  assert(result[2].transpose().row(0) == rand.row(4));
  assert(result[2].transpose().row(1) == rand.row(5));
  assert(result[3].transpose().row(0) == rand.row(6));
  assert(result[3].transpose().row(1) == rand.row(7));
}


//takes a vector of num_frames  3 x num_markers matrices and packs them
// into a single matrix of size (num_frames * num_markers) x 3
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 3> pack_points(VECTOROF3XDMATRICES<T> points, int num_frames, int num_markers) {
  Eigen::Matrix<T, -1, 3> packed = Eigen::Matrix<T, -1, 3>(num_frames * num_markers, 3);
  for(int i = 0; i < num_frames; i++){
    packed.template block(i * num_markers, 0, num_markers, 3) = points[i].transpose();
  }
  return packed;
}

template<typename T>
void unit_test_pack_points() {
  VECTOROF3XDMATRICES<T> a(3);
  a[0] = Eigen::Matrix<T, 3, 3>::Identity();
  a[1] = Eigen::Matrix<T, 3, 3>::Zero();
  a[2] = Eigen::Matrix<T, 3, 3>::Ones();
  auto test_mat = pack_points(a, a.size(), a[0].cols());
  for(int i = 0; i < 3; i ++)
    assert(test_mat.row(i, i) == 1 && test_mat.row(i).cwiseAbs().sum() == 1);
  for(int i = 3; i < 6; i ++)
    assert(test_mat.row(i).cwiseAbs().sum() == 0);
  for(int i = 6; i < 9; i ++)
    assert(test_mat(i, 0) == 1 && test_mat(i, 1) == 1 && test_mat(i, 0) == 1);
}


//helper function to recursively compute binomial coefficients
constexpr inline size_t binom(size_t n, size_t k) noexcept {
  if (k > n)
    return 0;
  if(k == 0 || k == n)
    return 1;
  if (k==1 || k==n-1)
    return n;
  if(k+k<n)
    return (binom(n - 1, k - 1) * n)/k;
  return (binom(n-1,k) * n)/(n-k);
}

//compute an n-th degree bernstein polynomial. If this is the first time calling this particular polynomial degree,
//compute the binomial coeffcicients and cache them in a global var for future use.
template<int degree, typename T>
T nth_degree_bernstein_polynomial(int k, T v) {
  if(!binomial_coeffs.contains(degree)){
    binomial_coeffs[degree] = std::vector<int>(degree + 1);
    for(int i = 0; i <= (degree + 1)/2; i++)
      binomial_coeffs[degree][i] = binom(degree, i);
    //we can just compute half of the array and relect it over because of how the binomial coefficients are computed (symmetric)
    for(int i = (degree + 1)/2; i < degree+1; i++)
      binomial_coeffs[degree][i] = binomial_coeffs[degree][degree - i];
  }
  //bernstein polynomial formula.
  return binomial_coeffs[degree][k] * pow(1 - v, degree - k) * pow(v, k);
}

template<typename T>
void unit_test_nth_degree_bernstein() {
  float x = 0.2;
  //check degree 5. coeffs from wolfram alpha.
  float arr[6] = {0.32768, 0.4096, 0.2048, 0.0512, 0.0064, 0.00032};
  for(int i = 0; i < 5 + 1; i++)
    assert(nth_degree_bernstein_polynomial<5>(i, x) == arr[i]);

  //check degree 4
  float arr2[6] = {0.4096, 0.4096, 0.1536, 0.0256, 0.0016};
  for(int i = 0; i < 4 + 1; i++)
    assert(nth_degree_bernstein_polynomial<4>(i, x) == arr2[i]);
}



//scale each column of a matrix between zero and one for use in a bernstein polynomial. Precompute the column-wise min and max values
// so that at "test time" we can use the same coefficents we did when generating the polynomial coefficients
template<typename T>
Eigen::Matrix<T, -1, 3> scale_to_unit_box(Eigen::Matrix<T, -1, 3> prescale, Eigen::Vector3<T> min_coeffs, Eigen::Vector3<T> max_coeffs) {
  Eigen::Matrix<T, -1, 3> postscale(prescale.rows(), 3);
  //make sure to check if min=max, so we dont divide by zero
  postscale.col(0) = ((prescale.col(0).array() - min_coeffs(0)) /
          ((abs(max_coeffs(0) - min_coeffs(0)) < 1e-5) ? 1 : (max_coeffs(0) - min_coeffs(0)))).matrix();
  postscale.col(1) = ((prescale.col(1).array() - min_coeffs(1)) /
          ((abs(max_coeffs(1) - min_coeffs(1)) < 1e-5) ? 1 : (max_coeffs(1) - min_coeffs(1)))).matrix();
  postscale.col(2) = ((prescale.col(2).array() - min_coeffs(2)) /
          ((abs(max_coeffs(2) - min_coeffs(2)) < 1e-5) ? 1 : (max_coeffs(2) - min_coeffs(2)))).matrix();
  return postscale;
}

template<typename T>
void unit_test_scale_to_unit_box() {
  DISTORTIONCOEFFICIENTS<5, T> a;
  Eigen::Matrix<T, 3, 3> prescale;
  prescale << -1, -2, -4,
               0, 0, 0,
               1, 2, 4;
  a.max_coeffs = prescale.colwise().maxCoeff();
  a.min_coeffs = prescale.colwise().minCoeff();
  Eigen::Matrix<T, 3, 3> postscale = scale_to_unit_box<T>(prescale, a.min_coeffs, a.max_coeffs);
  assert((postscale(0, 0) == 0) && (postscale(0, 1) == 0) && (postscale(0, 2) == 0));
  assert((postscale(1, 0) == 0.5) && (postscale(1, 1) == 0.5) && (postscale(1, 2) == 0.5));
  assert((postscale(2, 0) == 1) && (postscale(2, 1) == 1) && (postscale(2, 2) == 1));

  prescale = Eigen::Matrix<T, 3, 3>::Constant(5);
  a.max_coeffs = prescale.colwise().maxCoeff();
  a.min_coeffs = prescale.colwise().minCoeff();
  postscale = scale_to_unit_box<T>(prescale, a.min_coeffs, a.max_coeffs);
  for(int i = 0; i < 9; i++)
    assert(postscale(i) == 0);
}

#endif //PROGRAMMING_ASSIGMENT_ONE_DISTORTION_CORRECTION_H
