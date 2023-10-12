//
// Created by Aabha on 10/10/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_CAL_BODY_PARSER_H
#define PROGRAMMING_ASSIGMENT_ONE_CAL_BODY_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Dense"

//takes in a file name, spits out a list of matrices
template<typename T>
class Cal_Body_Parser {
    std::ifstream in;
    std::string line;
    char tmp;
    int n_D, n_A, n_C;
    std::string file_name_for_output;

    Eigen::Matrix<T, 3, Eigen::Dynamic> n_d_vals;
    Eigen::Matrix<T, 3, Eigen::Dynamic> n_a_vals;
    Eigen::Matrix<T, 3, Eigen::Dynamic> n_c_vals;

    void initialize_matrices();

public:
    Cal_Body_Parser(std::string_view file_name);
    ~Cal_Body_Parser();
    Eigen::Matrix<T, 3, Eigen::Dynamic> get_n_d_vals() {return n_d_vals;};
    Eigen::Matrix<T, 3, Eigen::Dynamic> get_n_a_vals() {return n_a_vals;};
    Eigen::Matrix<T, 3, Eigen::Dynamic> get_n_c_vals() {return n_c_vals;};
};

template<typename T>
void Cal_Body_Parser<T>::initialize_matrices() {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> this->n_D >> tmp;
    iss_line >> this->n_A >> tmp;
    iss_line >> this->n_C >> tmp;
    std::getline(iss_line, file_name_for_output);
    n_d_vals.resize(3, n_D);
    n_d_vals.setZero();
    n_a_vals.resize(3, n_A);
    n_a_vals.setZero();
    n_c_vals.resize(3, n_C);
    n_c_vals.setZero();
}

template<typename T>
Cal_Body_Parser<T>::~Cal_Body_Parser() {
    in.close();
}

template<typename T>
Cal_Body_Parser<T>::Cal_Body_Parser(std::string_view file_name) {
    in = std::ifstream(file_name.data());
    if (!in.is_open()) {
        throw std::runtime_error{"File does not exist"};
    }

    initialize_matrices();

    for (int i = 0; i < n_D; i++) {
        std::getline(in, line);
        std::istringstream iss_line(line);
        iss_line >> n_d_vals(0,i) >> tmp;
        iss_line >> n_d_vals(1, i) >> tmp;
        iss_line >> n_d_vals(2,i);
    }
    for (int i = 0; i < n_A; i++) {
        std::getline(in, line);
        std::istringstream iss_line(line);
        iss_line >> n_a_vals(0,i) >> tmp;
        iss_line >> n_a_vals(1, i) >> tmp;
        iss_line >> n_a_vals(2,i);
    }
    for (int i = 0; i < n_C; i++) {
        std::getline(in, line);
        std::istringstream iss_line(line);
        iss_line >> n_c_vals(0,i) >> tmp;
        iss_line >> n_c_vals(1, i) >> tmp;
        iss_line >> n_c_vals(2,i);
    }
}


#endif //PROGRAMMING_ASSIGMENT_ONE_CAL_BODY_PARSER_H
