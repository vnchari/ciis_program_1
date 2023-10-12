//
// Created by Aabha on 10/11/2023.
//

#ifndef PROGRAMMING_ASSIGMENT_ONE_OUTPUT_FILE_CREATOR_H
#define PROGRAMMING_ASSIGMENT_ONE_OUTPUT_FILE_CREATOR_H


#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/dense"
#include <iterator>
#include <list>
#include <iomanip>

template<typename T>
class Output_File_Creator {

public:
    Output_File_Creator(int n_C, int n_frames,
                        Eigen::Vector<T, 3> em_prob_pos, Eigen::Vector<T, 3> opt_prob_pos,
                        std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> c_expected_vals);
};

template<typename T>
Output_File_Creator<T>::Output_File_Creator(int n_C, int n_frames,
                                            Eigen::Vector<T, 3> em_prob_pos, Eigen::Vector<T, 3> opt_prob_pos,
                                            std::vector<Eigen::Matrix<T, 3, Eigen::Dynamic>> c_expected_vals) {
    std::ofstream output_file;
    output_file.open("C:\\Users\\Aabha\\Desktop\\CIS\\PA1\\ciis_program_1\\PA1-OUTPUT-1.txt");
    output_file << n_C << ", " << n_frames << ", NAME-OUTPUT1.TXT \n";
    output_file << "  " << em_prob_pos(0) << ",   " << em_prob_pos(1) << ",   " << em_prob_pos(2) << "\n";
    output_file << "  " << opt_prob_pos(0) << ",   " << opt_prob_pos(1) << ",   " << opt_prob_pos(2) << "\n";

    Eigen::Matrix<T, 3, Eigen::Dynamic> tmp;
    tmp.resize(3, n_C);
        output_file << std::fixed;
        output_file << std::setprecision(2);

        for (auto elem : c_expected_vals) {
            tmp = c_expected_vals[elem];
            for (int i = 0; i < n_C; i++) {
                output_file << "  " << tmp(0, i) << ",   " << tmp(1, i) << ",   " << tmp(2, i) << "\n";
            }
        }
    output_file.close();

}


#endif //PROGRAMMING_ASSIGMENT_ONE_OUTPUT_FILE_CREATOR_H
