//
// Created by Aabha on 10/11/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_OPT_PIVOT_PARSER_H
#define PROGRAMMING_ASSIGMENT_ONE_OPT_PIVOT_PARSER_H


#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/dense"
#include <iterator>
#include <list>
#include <vector>

template<typename T>
class Opt_Pivot_Parser {
    std::ifstream in;
    std::string line;
    char tmp;
    int n_D, n_H, n_Frames;
    std::string file_name_for_output;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_D_vals;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_H_vals;
    void initialize_matrices();

public:
    Opt_Pivot_Parser(std::string_view file_name);
    ~Opt_Pivot_Parser();
    void show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName);
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_D_vals() {return n_D_vals;};
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_H_vals() {return n_H_vals;};

};

template<typename T>
void Opt_Pivot_Parser<T>::show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName) {
    typename std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>>::iterator it;
    for (it = listName.begin(); it != listName.end(); ++it) {
        std::cout << '\n' << *it << '\n';
    }
}

template<typename T>
void Opt_Pivot_Parser<T>::initialize_matrices() {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> this->n_D >> tmp;
    iss_line >> this->n_H >> tmp;
    iss_line >> this->n_Frames >> tmp;
    std::getline(iss_line, file_name_for_output);
    //std::cout << n_D << " "  << n_H << " " << n_Frames << " " << file_name_for_output << std::endl;

}

template<typename T>
Opt_Pivot_Parser<T>::~Opt_Pivot_Parser() {
    in.close();
}

template<typename T>
Opt_Pivot_Parser<T>::Opt_Pivot_Parser(std::string_view file_name) {
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

        Eigen::Matrix<T, 3, Eigen::Dynamic> tmpH;
        tmpH.resize(3, this->n_H);
        tmpH.setZero();
        for (int i = 0; i < n_H; i++) {
            std::getline(in, line);
            std::istringstream iss_line(line);
            iss_line >> tmpH(0, i) >> tmp;
            iss_line >> tmpH(1, i) >> tmp;
            iss_line >> tmpH(2, i);
        }
        n_H_vals.push_back(tmpH);
    }
    //show_values(n_D_vals);
    //show_values(n_H_vals);
}

#endif //PROGRAMMING_ASSIGMENT_ONE_OPT_PIVOT_PARSER_H
