//
// Created by Vivek Chari on 10/12/23.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_READ_INPUT_FILES_H
#define PROGRAMMING_ASSIGMENT_ONE_READ_INPUT_FILES_H

#include <list>
#include <iterator>
#include "Eigen/StdVector"
#include "Eigen/Dense"
#include <sstream>
#include <fstream>
#include <iostream>

template<typename T>
struct OPTPIVOTDATA {
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_D_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_H_vals;
};

template<typename T>
struct CALREADINGSDATA {
    int n_frames, n_C;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_D_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_A_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_C_vals;
};

template<typename T>
struct EMPIVOTDATA {
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> n_G_vals;
};

template<typename T>
struct CALBODYDATA {
    Eigen::Matrix<T, 3, Eigen::Dynamic> n_d_vals;
    Eigen::Matrix<T, 3, Eigen::Dynamic> n_a_vals;
    Eigen::Matrix<T, 3, Eigen::Dynamic> n_c_vals;
};

template<typename T>
struct CALBODYDATA<T> read_calbody_file(std::string_view file_name) {
  std::ifstream in(file_name.data());
  if (!in.is_open()) {
    throw std::runtime_error{"File does not exist"};
  }

  std::string line;
  char tmp;
  int n_D, n_A, n_C;
  std::string file_name_for_output;

  std::getline(in, line);
  std::istringstream iss_line(line);
  iss_line >> n_D >> tmp;
  iss_line >> n_A >> tmp;
  iss_line >> n_C >> tmp;
  std::getline(iss_line, file_name_for_output);

  Eigen::Matrix<T, 3, Eigen::Dynamic> n_d_vals(3, n_D);
  n_d_vals.setZero();
  Eigen::Matrix<T, 3, Eigen::Dynamic> n_a_vals(3, n_A);
  n_a_vals.setZero();
  Eigen::Matrix<T, 3, Eigen::Dynamic> n_c_vals(3, n_C);
  n_c_vals.setZero();

  for (int i = 0; i < n_D; i++) {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> n_d_vals(0, i) >> tmp;
    iss_line >> n_d_vals(1, i) >> tmp;
    iss_line >> n_d_vals(2, i);
  }
  for (int i = 0; i < n_A; i++) {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> n_a_vals(0, i) >> tmp;
    iss_line >> n_a_vals(1, i) >> tmp;
    iss_line >> n_a_vals(2, i);
  }
  for (int i = 0; i < n_C; i++) {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> n_c_vals(0, i) >> tmp;
    iss_line >> n_c_vals(1, i) >> tmp;
    iss_line >> n_c_vals(2, i);
  }

  return CALBODYDATA<T> {n_d_vals, n_a_vals, n_c_vals};
}

template<typename T>
OPTPIVOTDATA<T> read_optpivot_data(std::string_view file_name) {
  OPTPIVOTDATA<T> data;
  std::ifstream in(file_name.data());
  if (!in.is_open()) {
    throw std::runtime_error{"File does not exist"};
  }

  std::string line;
  char tmp;
  int n_D, n_H, n_Frames;
  std::string file_name_for_output;

  std::getline(in, line);
  std::istringstream iss_line(line);
  iss_line >> n_D >> tmp;
  iss_line >> n_H >> tmp;
  iss_line >> n_Frames >> tmp;
  std::getline(iss_line, file_name_for_output);

  for (int i = 0; i < n_Frames; i++) {
    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpD;
    tmpD.resize(3, n_D);
    tmpD.setZero();
    for (int i = 0; i < n_D; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpD(0, i) >> tmp;
      iss_line >> tmpD(1, i) >> tmp;
      iss_line >> tmpD(2, i);
    }
    data.n_D_vals.push_back(tmpD);

    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpH;
    tmpH.resize(3, n_H);
    tmpH.setZero();
    for (int i = 0; i < n_H; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpH(0, i) >> tmp;
      iss_line >> tmpH(1, i) >> tmp;
      iss_line >> tmpH(2, i);
    }
    data.n_H_vals.push_back(tmpH);
  }

  return data;
}

template<typename T>
CALREADINGSDATA<T> read_calreadings_data(std::string_view file_name) {
  CALREADINGSDATA<T> data;

  std::ifstream in(file_name.data());
  if (!in.is_open()) {
    throw std::runtime_error{"File does not exist"};
  }

  std::string line;
  char tmp;
  int n_D, n_A, n_C, n_Frames;
  std::string file_name_for_output;

  std::getline(in, line);
  std::istringstream iss_line(line);
  iss_line >> n_D >> tmp;
  iss_line >> n_A >> tmp;
  iss_line >> n_C >> tmp;
  iss_line >> n_Frames >> tmp;
  std::getline(iss_line, file_name_for_output);

  for (int i = 0; i < n_Frames; i++) {
    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpD;
    tmpD.resize(3, n_D);
    tmpD.setZero();
    for (int i = 0; i < n_D; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpD(0, i) >> tmp;
      iss_line >> tmpD(1, i) >> tmp;
      iss_line >> tmpD(2, i);
    }
    data.n_D_vals.push_back(tmpD);

    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpA;
    tmpA.resize(3, n_A);
    tmpA.setZero();
    for (int i = 0; i < n_A; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpA(0, i) >> tmp;
      iss_line >> tmpA(1, i) >> tmp;
      iss_line >> tmpA(2, i);
    }
    data.n_A_vals.push_back(tmpA);

    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpC;
    tmpC.resize(3, n_C);
    tmpC.setZero();
    for (int i = 0; i < n_C; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpC(0, i) >> tmp;
      iss_line >> tmpC(1, i) >> tmp;
      iss_line >> tmpC(2, i);
    }
    data.n_C_vals.push_back(tmpC);
  }

  data.n_frames = n_Frames;
  data.n_C = n_C;
  return data;
}

template<typename T>
EMPIVOTDATA<T> read_empivot_data(std::string_view file_name) {
  EMPIVOTDATA<T> data;

  std::ifstream in(file_name.data());
  if (!in.is_open()) {
    throw std::runtime_error{"File does not exist"};
  }

  std::string line;
  char tmp;
  int n_G, n_Frames;
  std::string file_name_for_output;

  std::getline(in, line);
  std::istringstream iss_line(line);
  iss_line >> n_G >> tmp;
  iss_line >> n_Frames >> tmp;
  std::getline(iss_line, file_name_for_output);

  for (int i = 0; i < n_Frames; i++) {
    Eigen::Matrix<T, 3, Eigen::Dynamic> tmpG;
    tmpG.resize(3, n_G);
    tmpG.setZero();
    for (int i = 0; i < n_G; i++) {
      std::getline(in, line);
      std::istringstream iss_line(line);
      iss_line >> tmpG(0, i) >> tmp;
      iss_line >> tmpG(1, i) >> tmp;
      iss_line >> tmpG(2, i);
    }
    data.n_G_vals.push_back(tmpG);
  }

  return data;
}

#endif //PROGRAMMING_ASSIGMENT_ONE_READ_INPUT_FILES_H
