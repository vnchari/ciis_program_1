//
// Created by Vivek Chari on 10/10/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
#define PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H

#include <iostream>
#include <list>
#include "Eigen/Dense"

template<typename T>
class Frame {
    Eigen::Matrix<T, 4, 4> frame;
    Eigen::Matrix<T, 4, 4> iframe;
public:
    Frame(Eigen::Matrix<T, 3, 3> rot, Eigen::Vector<T, 3> pos);
    void invert();
    Eigen::Vector<T,3> apply_forward(Eigen::Vector<T,3>);
    Eigen::Vector<T,3> apply_inverse(Eigen::Vector<T,3>);
    Eigen::Matrix<T, 4, 4> get_frame() {return frame;};
    Eigen::Matrix<T, 4, 4> get_iframe() {return iframe;};
};

template<typename T>
Eigen::Vector<T, 3> Frame<T>::apply_inverse(Eigen::Vector<T, 3> vec) {
  Eigen::Vector<T, 4> tmp {vec[0], vec[1], vec[2], 1};
  return (iframe * tmp).template head<3>();
}

template<typename T>
Eigen::Vector<T, 3> Frame<T>::apply_forward(Eigen::Vector<T, 3> vec) {
  Eigen::Vector<T, 4> tmp {vec[0], vec[1], vec[2], 1};
  return (frame * tmp).template head<3>();
}

template<typename T>
void Frame<T>::invert() {
  iframe = frame;
  iframe.transposeInPlace();
  iframe.template block<3,1>(0, 3) =
          -1 * iframe.template block<3,3>(0, 0) * iframe.template block<1, 3>(3, 0).transpose();
  iframe.template block<1, 3>(3, 0).setZero();
}


template<typename T>
Frame<T>::Frame(Eigen::Matrix<T, 3, 3> rot, Eigen::Vector<T, 3> pos) {
  frame.setIdentity();
  frame.template block<3,3>(0, 0) = rot;
  frame.template block<3,1>(0, 3) = pos;
}


template<typename T>
Eigen::Matrix<T,3, 3> gen_random_orthogonal(){
  return Eigen::Quaternion<T>::UnitRandom().toRotationMatrix();
}

#endif //PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
