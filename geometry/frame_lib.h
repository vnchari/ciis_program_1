//
// Created by Vivek Chari on 10/10/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
#define PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H

#include <Eigen/Dense>

template<typename T>
class Frame {
    Eigen::Matrix<T, 4, 4> frame;
};


template<typename T>
Eigen::Matrix<T,3, 3> gen_random_orthogonal(){
  return Eigen::Quaternion<T>::UnitRandom().toRotationMatrix();
}




#endif //PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
