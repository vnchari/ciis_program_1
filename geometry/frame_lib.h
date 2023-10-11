//
// Created by Vivek Chari on 10/10/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
#define PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H

#include <iostream>
#include "Eigen/Dense"

template<typename T>
class Frame {
    /* Frames will be defined as 4 x 4 matrices consisting of the rotational component composing the top left 3x3 elements
        and the translational component composing the elements (1,4), (2,4), and (3,4). the values of the last row
        are not relevant, but element (4,4) must be 1 to properly apply frame transformations to positions.
    */
    Eigen::Matrix<T, 4, 4> frame; // 4x4 frame matrix
    Eigen::Matrix<T, 4, 4> inverse_frame; // 4x4 inverse frame matrix

public:
    Frame(Eigen::Matrix<T, 3, 3> rotation, Eigen::Matrix<T, 3, 1> translation); //public frame constructor
    void invert(); // Inverts the frame
    Eigen::Matrix<T, 4, 4> get_frame() {return frame;}; // Public frame getter
    Eigen::Matrix<T, 4, 4> get_inverse_frame() {return inverse_frame;}; // Public inverse_Frame getter
    Eigen::Vector<T, 3> forward_transform(Eigen::Vector<T, 3> startingPos); // Method to apply a forward transformation
    Eigen::Vector<T, 3> inverse_transform(Eigen::Vector<T, 3> startingPos); // Method to apply an invese transformation
};


/*
 * To apply the frame transformation using the 4 x 4 Frame matrix, we simply convert the startingPos 3x1 vector into a
 * 4x1 vector, multiply the frame transformation by this vector, and remove the last element.
 */
template<typename T>
Eigen::Vector<T, 3> Frame<T>::inverse_transform(Eigen::Vector<T, 3> startingPos) {
    Eigen::Vector<T, 4> temp = (startingPos[0], startingPos[1], startingPos[2], 1);
    return (this->inverse_frame * temp).template head<3>();
}


/*
 * To apply the frame transformation using the 4 x 4 Frame matrix, we simply convert the startingPos 3x1 vector into a
 * 4x1 vector, multiply the frame transformation by this vector, and remove the last element.
 */
template<typename T>
Eigen::Vector<T, 3> Frame<T>::forward_transform(Eigen::Vector<T, 3> startingPos) {
    Eigen::Vector<T, 4> temp = (startingPos[0], startingPos[1], startingPos[2], 1);
    return (this->frame * temp).template head<3>();
}

template<typename T>
void Frame<T>::invert() {
    this->inverse_frame = this->frame; // Start by setting the inverse frame to frame. This allows efficient operations

    /* Because the rotation matrix is an SO3 matrix, its inverse it its transpose
       we can transpose the inverse matrix and then recompute the translational component. */
    this->inverse_frame.transposeInPlace(); // This sets the rotation component to its inverse already

    // Compute the inverse translational component as -1*R^-1*p.
    this->inverse_frame.template block<3,1>(0,3) =
            -1 * this->inverse_frame.template block<3,3>(0,0)*this->inverse_frame.template block<1,3>(3,0).transpose();
    this->inverse_frame.template block<1,3>(3,0).setZero(); //sets first 3 elements in last row to 0.
}


template<typename T>
Frame<T>::Frame(Eigen::Matrix<T, 3, 3> rotation, Eigen::Matrix<T, 3, 1> translation) {
    this->frame.setIdentity();
    this->frame.template block<3,3>(0,0) = rotation;
    this->frame.template block<3,1>(0,3) = translation;

    std::cout << this->frame << std::endl;
}


template<typename T>
Eigen::Matrix<T,3, 3> gen_random_orthogonal(){
  return Eigen::Quaternion<T>::UnitRandom().toRotationMatrix();
}




#endif //PROGRAMMING_ASSIGMENT_ONE_FRAME_LIB_H
