//
// Created by Vivek Chari on 10/12/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_CALIBRATION_H
#define PROGRAMMING_ASSIGMENT_ONE_CALIBRATION_H

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "frame_lib.h"

template<typename T>
Eigen::Vector<T, 6> solve_least_squares(const std::vector<Eigen::Matrix<T, 4, 4>,
        Eigen::aligned_allocator<Eigen::Matrix<T, 4, 4>>> &vec) {
  Eigen::Matrix<T, -1, 6> A = Eigen::Matrix<T, -1, 6>(vec.size() * 3,6);
  Eigen::Vector<T, -1> b = Eigen::Vector<T, -1>(vec.size() * 3);

  for(size_t i = 0; i < vec.size(); i++){
    A.template block<3, 3>(3 * i, 0) = vec[i].template topLeftCorner<3,3>();
    A.template block<3, 3>(3 * i, 3) = -1 * Eigen::Matrix<T, 3, 3>::Identity();
    b.template block<3, 1>(3 * i, 0) = -1 * vec[i].template topRightCorner<3,1>();
  }
  return A.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
}

template<typename T>
std::pair<Eigen::Vector3<T>, Eigen::Vector3<T>> pivot_calibration_routine(const std::vector<Eigen::Matrix<T, -1, 3>,
        Eigen::aligned_allocator<Eigen::Matrix<T, -1, 3>>> &obs){
  FrameGraph<T> graph(Registration::PROCRUSTES);
  std::vector<Eigen::Matrix<T, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<T, -1, 3>>>
          centered_obs = std::vector<Eigen::Matrix<T, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<T, -1, 3>>>(obs.size());
  Eigen::Matrix<T, -1, 3> g = obs[0].rowwise() - obs[0].colwise().mean(); //row vector

  std::vector<std::pair<std::string, std::string>> frame_pairs(obs.size());
  for(size_t i = 0; i < obs.size(); i++) {
    graph.register_transform("POINTER_LOCAL_FRAME", g,
                             "EM_TRACKER_" + std::to_string(i), obs[i]);
    frame_pairs[i] = std::pair<std::string, std::string> ("POINTER_LOCAL_FRAME", "EM_TRACKER_" + std::to_string(i));
  }
  auto sol = solve_least_squares( graph.retrieve_transform_matrices(frame_pairs));
  //first 3 coeffs of vector are the location of the tip of the pointer in local pointer frame.
  return std::pair<Eigen::Vector3<T>, Eigen::Vector3<T>>(sol.template head<3>(), sol.template tail<3>());
}

#endif //PROGRAMMING_ASSIGMENT_ONE_CALIBRATION_H
