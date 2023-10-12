//
// Created by Aabha on 10/11/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_CAL_READINGS_PARSER_H
#define PROGRAMMING_ASSIGMENT_ONE_CAL_READINGS_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Dense"
#include <iterator>
#include <list>
#include <vector>

template<typename T>
class Cal_Readings_Parser {
    std::ifstream in;
    std::string line;
    char tmp;
    int n_D, n_A, n_C, n_Frames;
    std::string file_name_for_output;

    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_D_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_A_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_C_vals;
    void initialize_matrices();

public:
    Cal_Readings_Parser(std::string_view file_name);
    ~Cal_Readings_Parser();
    void show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName);
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_D_vals() {return n_D_vals;};
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_A_vals() {return n_A_vals;};
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_C_vals() {return n_C_vals;};

};

template<typename T>
void Cal_Readings_Parser<T>::show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName) {
    typename std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>>::iterator it;
    for (it = listName.begin(); it != listName.end(); ++it) {
        std::cout << '\n' << *it;
    }
}

template<typename T>
void Cal_Readings_Parser<T>::initialize_matrices() {

    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> this->n_D >> tmp;
    iss_line >> this->n_A >> tmp;
    iss_line >> this->n_C >> tmp;
    iss_line >> this->n_Frames >> tmp;
    std::getline(iss_line, file_name_for_output);
    //std::cout << n_D << " " << n_A << " " << n_C << " " << n_Frames << " " << file_name_for_output << std::endl;
}

template<typename T>
Cal_Readings_Parser<T>::~Cal_Readings_Parser() {
    in.close();
}

template<typename T>
Cal_Readings_Parser<T>::Cal_Readings_Parser(std::string_view file_name) {
    in = std::ifstream(file_name.data());
    if (!in.is_open()) {
        throw std::runtime_error{"File does not exist"};
    }

    initialize_matrices();

    for (int i = 0; i < n_Frames; i++) {
        Eigen::Matrix<T, 3, Eigen::Dynamic> tmpD;
        tmpD.resize(3, this->n_D);
        tmpD.setZero();
        for (int i = 0; i < n_D; i++) {
            std::getline(in, line);
            std::istringstream iss_line(line);
            iss_line >> tmpD(0, i) >> tmp;
            iss_line >> tmpD(1, i) >> tmp;
            iss_line >> tmpD(2, i);
        }
        n_D_vals.push_back(tmpD);

        Eigen::Matrix<T, 3, Eigen::Dynamic> tmpA;
        tmpA.resize(3, this->n_A);
        tmpA.setZero();
        for (int i = 0; i < n_A; i++) {
            std::getline(in, line);
            std::istringstream iss_line(line);
            iss_line >> tmpA(0, i) >> tmp;
            iss_line >> tmpA(1, i) >> tmp;
            iss_line >> tmpA(2, i);
        }
        n_A_vals.push_back(tmpA);

        Eigen::Matrix<T, 3, Eigen::Dynamic> tmpC;
        tmpC.resize(3, this->n_C);
        tmpC.setZero();
        for (int i = 0; i < n_C; i++) {
            std::getline(in, line);
            std::istringstream iss_line(line);
            iss_line >> tmpC(0, i) >> tmp;
            iss_line >> tmpC(1, i) >> tmp;
            iss_line >> tmpC(2, i);
        }
        n_C_vals.push_back(tmpC);
    }
}


#endif //PROGRAMMING_ASSIGMENT_ONE_CAL_READINGS_PARSER_H
