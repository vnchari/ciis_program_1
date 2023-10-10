#include <chrono>
#include <iostream>
#include "geometry/frame_lib.h"
#include <chrono>

int main() {
  auto clock = std::chrono::high_resolution_clock();
  auto rand_orth = gen_random_orthogonal<double>();
  auto pos = Eigen::Vector3d{1, 1, 1};
  auto begin = clock.now();
  Frame<double> a = Frame<double> (rand_orth, pos);
  a.invert();
  a.forward(pos);
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(clock.now() - begin).count() << std::endl;
//  std::cout << a.get_frame() << std::endl;
//  std::cout << a.get_iframe() << std::endl;
  return 0;
}
