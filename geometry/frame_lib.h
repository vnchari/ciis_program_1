//
// Created by Vivek Chari on 10/10/23.
//
// This file contains methods and classes that contruct "frame graphs"
// that describe the relationship between different frames. It contains
// a class FrameTransformation that can be used to store frame transformations
// in a packed 4x4 matrix, along with methods to apply the transform to vectors
// and matrices, as well as methods to compute (and apply) the inverse transform.
// The FrameGraph class provides an interface to easily compute and apply
// transforms between frames, and has an internal BFS that computes the shortest path
// between any frames in the graph.

#ifndef PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
#define PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H

#include <map>
#include <set>
#include <list>
#include <iostream>
#include "register.h"
#include "Eigen/Dense"
#include "Eigen/StdVector"

// class to house frame, inverse frame, and frame application methods.
template<typename T>
class FrameTransformation {
    Eigen::Matrix<T, 4, 4> frame;
    Eigen::Matrix<T, 4, 4> iframe;
public:
    FrameTransformation(Eigen::Matrix<T, 3, 3> rot, Eigen::Vector<T, 3> pos);
    void compute_inverse();
    Eigen::Vector<T,3> apply_forward_transform(Eigen::Vector<T,3> vec);
    Eigen::Vector<T,3> apply_inverse_transform(Eigen::Vector<T,3> vec);
    Eigen::Matrix<T,-1, 3> apply_forward_transform(Eigen::Matrix<T,-1, 3> mat);
    Eigen::Matrix<T,-1, 3> apply_inverse_transform(Eigen::Matrix<T,-1, 3> packed_pts);
    Eigen::Matrix<T, 4, 4> get_frame() {return frame;};
    Eigen::Matrix<T, 4, 4> get_iframe() {return iframe;};
};

// this is not very efficient. it will be refactored into a single packed_pts multiply later on
template<typename T>
Eigen::Matrix<T,-1, 3> FrameTransformation<T>::apply_inverse_transform(const Eigen::Matrix<T, -1, 3> packed_pts) {
  Eigen::Matrix<T, -1, 3> tmp = Eigen::Matrix<T, -1, 3>::Ones(packed_pts.rows(), 3);
  //for every point in packed_pts, compute the inverse transform and store it
  for(int i = 0; i < packed_pts.rows(); i++) {
    const auto vec = packed_pts.row(i);
    Eigen::Vector<T, 4> tmpvec =  {vec[0], vec[1], vec[2], 1};
    tmp.row(i) = (iframe * tmpvec).template head<3>();
  }
  return tmp;
}

// this is not very efficient. it will be refactored into a single mat multiply later on
template<typename T>
Eigen::Matrix<T,-1, 3> FrameTransformation<T>::apply_forward_transform(const Eigen::Matrix<T, -1, 3> mat) {
  Eigen::Matrix<T, -1, 3> tmp = Eigen::Matrix<T, -1, 3>::Ones(mat.rows(), 3);
  //for every point in packed_pts, compute the transform and store it
  for(int i = 0; i < mat.rows(); i++) {
    const auto vec = mat.row(i);
    Eigen::Vector<T, 4> tmpvec =  {vec[0], vec[1], vec[2], 1};
    tmp.row(i) = (frame * tmpvec).template head<3>();
  }
  return tmp;
}

template<typename T>
Eigen::Vector<T, 3> FrameTransformation<T>::apply_inverse_transform(const Eigen::Vector<T, 3> vec) {
  Eigen::Vector<T, 4> tmp {vec[0], vec[1], vec[2], 1};
  return (iframe * tmp).template head<3>();
}

template<typename T>
Eigen::Vector<T, 3> FrameTransformation<T>::apply_forward_transform(const Eigen::Vector<T, 3> vec) {
  Eigen::Vector<T, 4> tmp {vec[0], vec[1], vec[2], 1};
  return (frame * tmp).template head<3>();
}

template<typename T>
void FrameTransformation<T>::compute_inverse() {
  //the 4x4 matrix "frame" stores the transform
  //transpose the top left 3x3 block (recall that R^{-1} = R^T)
  iframe = frame.transpose();
  //compute t^{-1} = -R^{-1}t
  iframe.template block<3,1>(0, 3) =
          -1 * iframe.template block<3,3>(0, 0) * iframe.template block<1, 3>(3, 0).transpose();
  //make sure to zero the bottom row of the packed transform. it is now [0,0,0,1]
  iframe.template block<1, 3>(3, 0).setZero();
}


template<typename T>
FrameTransformation<T>::FrameTransformation(Eigen::Matrix<T, 3, 3> rot, Eigen::Vector<T, 3> pos) {
  frame.setIdentity();
  //store the frame transform packed in a 4x4 matrix. this just sets the appropriate blocks in the matrix.
  frame.template block<3, 3>(0, 0)= rot;
  frame.template block<3, 1>(0, 3) = pos;
  compute_inverse();
}

//helper method to test implementations. this generates a random 3x3 orthogonal matrix.
template<typename T>
Eigen::Matrix<T,3, 3> gen_random_orthogonal(){
  return Eigen::Quaternion<T>::UnitRandom().toRotationMatrix();
}

//class to hold a graph of frame transforms, and call appropriate ones when needed.
template<typename T>
class FrameGraph {
    //edge class to store information on a transform.
    struct edge {
        std::string dest;
        FrameTransformation<T> transform;
        //we need to know whether this edge is an inverse, because the same transform object representes both
        // A->B and B->A
        bool is_inverse = false;
    };
    std::map<std::string, std::vector<edge>> graph;
    uint64_t registration_method = -1;
    //helper vars for BFS.
    std::map<std::string, std::pair<std::string, uint64_t>> prev_and_dist;
    std::set<std::string> set;
  public:
      explicit FrameGraph(uint64_t reg_method) { registration_method = reg_method; }
      void _register_transform_procrustes(const std::string& frame_A_name, Eigen::Matrix<T, -1, 3> points_in_A, const std::string& frame_B_name, Eigen::Matrix<T, -1, 3> points_in_B);
      void _register_transform_cpd(const std::string& frame_A_name, Eigen::Matrix<T, -1, 3> points_in_A, const std::string& frame_B_name, Eigen::Matrix<T, -1, 3> points_in_B);
      void register_transform(std::string frame_A_name, Eigen::Matrix<T, -1, 3> points_in_A, std::string frame_B_name, Eigen::Matrix<T, -1, 3> points_in_B);
      //returns points in B
      Eigen::Matrix<T, -1, 3> apply_direct_transform(const std::string& frame_A_name, const std::string& frame_B_name, Eigen::Matrix<T, -1, 3> points_in_A);
      Eigen::Matrix<T, -1, 3> apply_transform(const std::string& frame_A_name, const std::string& frame_B_name, Eigen::Matrix<T, -1, 3> points_in_A);
      //probably RVO-d out.
      std::vector<Eigen::Matrix<T, 4, 4>, Eigen::aligned_allocator<Eigen::Matrix<T, 4, 4>>> retrieve_transform_matrices(std::vector<std::pair<std::string, std::string>> &transforms);
      void clear() {graph.clear(); set.clear(); prev_and_dist.clear();}
};


//given a vector of edges (A,B),(C, D),.... retrieve the packed transform matrices for each of these edges, respecting
// whether the edge represents a forward or inverse transform. Mainly for use in pivot calibration.
template<typename T>
std::vector<Eigen::Matrix<T, 4, 4>, Eigen::aligned_allocator<Eigen::Matrix<T, 4, 4>>>
FrameGraph<T>::retrieve_transform_matrices(std::vector<std::pair<std::string, std::string>> &transforms) {
  std::vector<Eigen::Matrix<T, 4, 4>, Eigen::aligned_allocator<Eigen::Matrix<T, 4, 4>>> res(transforms.size());
  for(size_t i = 0; i < transforms.size(); i++){
    for(edge e : graph[transforms[i].first]) {
      if (e.dest == transforms[i].second) {
        res[i] = (e.is_inverse) ? e.transform.get_iframe() : e.transform.get_frame();
        break;
      }
    }
  }
  return res;
}

//given points in frame "A" and a frame "B" compute the transform "B from A". The points returned will now be in the
// coordinate frame "B". This methods traverses the internal graph using Djikstra's method to find the shortest
//path connecting the two frames, and then computes the appropriate transform by chaining transforms along
// traversed edges.
template<typename T>
Eigen::Matrix<T, -1, 3> FrameGraph<T>::apply_transform(const std::string& frame_A_name, const std::string& frame_B_name,
                                                       Eigen::Matrix<T, -1, 3> points_in_A) {
  //helper vars for Djikstra's
  prev_and_dist.clear();
  set.clear();
  for(auto a : graph) {
    prev_and_dist[a.first] = std::pair<std::string, uint64_t>("", INFINITY);
    set.insert(a.first);
  }
  prev_and_dist[frame_A_name] = std::pair<std::string, uint64_t>("", 0);
  //do Djikstra BFS. Not the most efficient implementation, but it works.
  //TODO: Refactor with focus on speed
  while (!set.empty()) {
      std::string best; T dist = INFINITY;
      for(const auto& node : set) {
        auto vals = prev_and_dist[node];
        if (vals.second < dist) {
          best = node;
          dist = vals.second;
        }
      }
      set.erase(best);
      for(edge neighbor : graph[best]){
        if(set.contains(neighbor.dest)){
          auto alt = prev_and_dist[best].second + 1;
          if (alt < prev_and_dist[neighbor.dest].second){
            prev_and_dist[neighbor.dest].second = alt;
            prev_and_dist[neighbor.dest].first = best;
          }
        }
        if(neighbor.dest == frame_B_name) {
          set.clear();
          break;
        }
      }
  }
  //reconstruct the shortest path taken
  std::list<std::string> transform_chain;
  auto target = frame_B_name;
  while(target != frame_A_name) {
    transform_chain.push_front(target);
    target = prev_and_dist[target].first;
  }
  transform_chain.push_front(frame_A_name);
  Eigen::Matrix<T, -1, 3> res = points_in_A;
  size_t sz = transform_chain.size();
  //traverse the shortest path, applying transforms as we go along the path.
  for(size_t i = 0; i < sz; i++) {
    auto cur = transform_chain.front();
    transform_chain.pop_front();
    for (edge e : graph[cur])
      if(e.dest == transform_chain.front()) {
        res = (e.is_inverse) ? e.transform.apply_inverse_transform(res) :
              e.transform.apply_forward_transform(res);
        break;
      }
  }
  //points are now in frame B.
  return res;
}

//given two frames, and corresponding point sets, compute a transform and register it in the graph.
//call into different registration subroutines depending on user specification.
template<typename T>
void FrameGraph<T>::register_transform(std::string frame_A_name, Eigen::Matrix<T, -1, 3> points_in_A,
                                       std::string frame_B_name, Eigen::Matrix<T, -1, 3> points_in_B) {
  if(registration_method == Registration::PROCRUSTES)
    return _register_transform_procrustes(frame_A_name, points_in_A, frame_B_name, points_in_B);
  if (registration_method == Registration::COHERENTPOINTDRIFT)
    return _register_transform_cpd(frame_A_name, points_in_A, frame_B_name, points_in_B);
  throw std::invalid_argument("Registration Method invalid.");
}


//given points in frame "A" and a frame "B" compute the tranform "B from A". The points returned will now be in the
// coordinate frame "B".  This method assumes a direct connection between A and B
template<typename T>
Eigen::Matrix<T, -1, 3> FrameGraph<T>::apply_direct_transform(const std::string& frame_A_name, const std::string& frame_B_name,
                                                             Eigen::Matrix<T, -1, 3> points_in_A) {
  if(!graph.contains(frame_A_name) || !graph.contains(frame_B_name))
    throw std::invalid_argument("either frame A or frame B doesn't exist");
  const auto A_edges = graph[frame_A_name];
  for(edge edge: A_edges)
    if(edge.dest == frame_B_name)
      return (edge.is_inverse) ? edge.transform.apply_inverse_transform(points_in_A) :
      edge.transform.apply_forward_transform(points_in_A);
  throw std::invalid_argument("Frame A and B don't have a direct edge");
}

template<typename T>
void FrameGraph<T>::_register_transform_procrustes(const std::string& frame_A_name, Eigen::Matrix<T, -1, 3> points_in_A,
                                                   const std::string& frame_B_name,
                                                   Eigen::Matrix<T, -1, 3> points_in_B) {
  Registration::Procrustes<T> registration(points_in_A, points_in_B);
  if(!graph.contains(frame_A_name)) {
    graph[frame_A_name] = std::vector<edge>();
  } else {
    for(edge edge : graph[frame_A_name])
      if(edge.dest == frame_B_name)
        throw std::invalid_argument("Transform A<->B already exists.");
  }
  if(!graph.contains(frame_B_name))
    graph[frame_B_name] =  std::vector<edge>();

  FrameTransformation<T> t(registration.B, registration.t);
  graph.at(frame_A_name).emplace_back(frame_B_name, t, false);
  graph.at(frame_B_name).emplace_back(frame_A_name, t, true);
}

template<typename T>
void FrameGraph<T>::_register_transform_cpd(const std::string& frame_A_name,
                                            Eigen::Matrix<T, -1, 3> points_in_A,
                                            const std::string& frame_B_name,
                                            Eigen::Matrix<T, -1, 3> points_in_B) {
  Registration::CoherentPointDrift<T> registration(points_in_A, points_in_B);
  if(!graph.contains(frame_A_name))
    graph[frame_A_name] =  std::vector<edge>();
  if(!graph.contains(frame_B_name))
    graph[frame_B_name] =  std::vector<edge>();

  FrameTransformation<T> t(registration.B, registration.t);
  graph.at(frame_A_name).emplace_back(frame_B_name, t, false);
  graph.at(frame_B_name).emplace_back(frame_A_name, t, true);
}




#endif //PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
