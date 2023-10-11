#include <chrono>
#include <iostream>
#include "geometry/frame_lib.h"
#include "geometry/cal_body_parser.h"
#include "geometry/cal_readings_parser.h"
#include "geometry/em_pivot_parser.h"
#include "geometry/opt_pivot_parser.h"
#include "geometry/output_file_creator.h"

int main() {
    std::string_view cal_body =
            R"(C:\Users\Aabha\Desktop\CIS\PA1\ciis_program_1\PA1-Student-Data\pa1-debug-b-calbody.txt)";
    auto cal_body_parser = Cal_Body_Parser<double> (cal_body);
    std::string_view cal_readings =
            R"(C:\Users\Aabha\Desktop\CIS\PA1\ciis_program_1\PA1-Student-Data\pa1-debug-b-calreadings.txt)";
    auto cal_readings_parser = Cal_Readings_Parser<double> (cal_readings);
    std::string_view em_pivot =
            R"(C:\Users\Aabha\Desktop\CIS\PA1\ciis_program_1\PA1-Student-Data\pa1-debug-b-empivot.txt)";
    auto em_pivot_parser = Em_Pivot_Parser<double> (em_pivot);
    std::string_view opt_pivot =
            R"(C:\Users\Aabha\Desktop\CIS\PA1\ciis_program_1\PA1-Student-Data\pa1-debug-b-optpivot.txt)";
    auto opt_pivot_parser = Opt_Pivot_Parser<double> (opt_pivot);
    //std::cout << opt_pivot_parser.get_n_D_vals().front() << '\n';
    //std::cout << opt_pivot_parser.get_n_H_vals().front();

    Output_File_Creator<double> (cal_body_parser.get_n_c_vals().cols(),cal_readings_parser.get_n_C_vals().size(),
                                 Eigen::Vector3d(0,0,1),
                                 Eigen::Vector3d(0,0,1), cal_readings_parser.get_n_C_vals());
    return 0;
}
