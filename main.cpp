#include <chrono>
#include <iostream>
#include "geometry/frame_lib.h"
#include "geometry/register.h"

int main() {
  auto clock = std::chrono::high_resolution_clock();
  constexpr size_t num_pts = 500;
  Eigen::Matrix<double, num_pts, 3> matA = Eigen::Matrix<double, num_pts, 3>::Random();
  Eigen::Matrix<double, num_pts, 3> matB = (matA * gen_random_orthogonal<double>().transpose()).rowwise() + Eigen::Vector<double, 3>::Random().transpose();
  Eigen::Matrix<double, num_pts, 3> matC = (matB * gen_random_orthogonal<double>().transpose()).rowwise() + Eigen::Vector<double, 3>::Random().transpose();
  Eigen::Matrix<double, num_pts, 3> matD = (matC * gen_random_orthogonal<double>().transpose()).rowwise() + Eigen::Vector<double, 3>::Random().transpose();

  auto begin = clock.now();
  FrameGraph<double> graph(Registration::EXTENDEDPROCRUSTES);
  graph.register_transform("A", matA, "B", matB);
  graph.register_transform("B", matB, "C", matC);
  graph.register_transform("C", matC, "D", matD);
  std::cout << (matD - graph.apply_transform("A", "D", matA)).cwiseAbs().mean() << std::endl;
//  std::cout << "Time elapsed (us): " << std::chrono::duration_cast<std::chrono::microseconds>(clock.now() - begin).count() << std::endl;
//  std::cout << (matD - graph.apply_transform("A", "D", matA)).cwiseAbs().mean();
  return 0;
}