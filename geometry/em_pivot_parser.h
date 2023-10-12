//
// Created by Aabha on 10/11/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_EM_PIVOT_PARSER_H
#define PROGRAMMING_ASSIGMENT_ONE_EM_PIVOT_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/dense"
#include <iterator>
#include <list>
#include <vector>

template<typename T>
class Em_Pivot_Parser {
    std::ifstream in;
    std::string line;
    char tmp;
    int n_G, n_Frames;
    std::string file_name_for_output;
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> n_G_vals;
    void initialize_matrices();

public:
    Em_Pivot_Parser(std::string_view file_name);
    ~Em_Pivot_Parser();
    void show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName);
    std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> get_n_G_vals() {return n_G_vals;};

};

template<typename T>
void Em_Pivot_Parser<T>::show_values(std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> listName) {
    typename std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>>::iterator it;
    for (it = listName.begin(); it != listName.end(); ++it) {
        std::cout << '\n' << *it << '\n';
    }
}

template<typename T>
void Em_Pivot_Parser<T>::initialize_matrices() {
    std::getline(in, line);
    std::istringstream iss_line(line);
    iss_line >> this->n_G >> tmp;
    iss_line >> this->n_Frames >> tmp;
    std::getline(iss_line, file_name_for_output);
    //std::cout << n_G << " "  << n_Frames << " " << file_name_for_output << std::endl;

}

template<typename T>
Em_Pivot_Parser<T>::~Em_Pivot_Parser() {
    in.close();
}

template<typename T>
Em_Pivot_Parser<T>::Em_Pivot_Parser(std::string_view file_name) {
    in = std::ifstream(file_name.data());
    if (!in.is_open()) {
        throw std::runtime_error{"File does not exist"};
    }

    initialize_matrices();

    for (int i = 0; i < n_Frames; i++) {
        Eigen::Matrix<T, 3, Eigen::Dynamic> tmpG;
        tmpG.resize(3, this->n_G);
        tmpG.setZero();
        for (int i = 0; i < n_G; i++) {
            std::getline(in, line);
            std::istringstream iss_line(line);
            iss_line >> tmpG(0, i) >> tmp;
            iss_line >> tmpG(1, i) >> tmp;
            iss_line >> tmpG(2, i);
        }
        n_G_vals.push_back(tmpG);
    }
    //show_values(n_G_vals);
}

#endif //PROGRAMMING_ASSIGMENT_ONE_EM_PIVOT_PARSER_H
