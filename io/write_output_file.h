//
// Created by Aabha on 10/11/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_WRITE_OUTPUT_FILE_H
#define PROGRAMMING_ASSIGMENT_ONE_WRITE_OUTPUT_FILE_H


#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <iterator>
#include <list>
#include <iomanip>

template<typename T>
void write_output_file(const std::string& filename,
                       size_t n_C, size_t n_frames,
                       Eigen::Vector<T, 3> em_prob_pos, Eigen::Vector<T, 3> opt_prob_pos,
                       std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<T, 3, Eigen::Dynamic>>> c_expected_vals) {
  std::ofstream output_file;
  output_file.open(filename);
  output_file << n_C << ", " << n_frames << ", NAME-OUTPUT1.TXT" << std::endl;
  output_file << std::fixed;
  output_file << std::setprecision(2);
  output_file << "  " << em_prob_pos(0) << ",   " << em_prob_pos(1) << ",   " << em_prob_pos(2) << "\n";
  output_file << "  " << opt_prob_pos(0) << ",   " << opt_prob_pos(1) << ",   " << opt_prob_pos(2) << "\n";
  Eigen::Matrix<T, 3, Eigen::Dynamic> tmp;
  tmp.resize(3, n_C);
  for (auto elem : c_expected_vals)
    for (size_t i = 0; i < n_C; i++)
      output_file << "  " << elem(0, i) << ",   " << elem(1, i) << ",   " << elem(2, i) << "\n";
  output_file.close();
}


#endif //PROGRAMMING_ASSIGMENT_ONE_WRITE_OUTPUT_FILE_H
