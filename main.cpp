#include <chrono>
#include <iostream>
#include "Eigen/StdVector"
#include "geometry/frame_lib.h"
#include "geometry/register.h"
#include "geometry/calibration.h"
#include "io/read_input_files.h"
#include "io/write_output_file.h"

#define PRECISION double

int main(int argc, char ** argv) {
  if(argc != 3){
    std::cerr << "Malformed arguments. arg 1: dir and name prefix of input files (e.g. 'input_files/pa1-debug-a'). "
                 "arg 2: dir of output_file";
    exit(1);
  }
  auto testname = std::string(argv[1]);
  auto outname = std::string(argv[2]);
  FrameGraph<PRECISION> graph(Registration::PROCRUSTES);
  const auto calbodydata = read_calbody_file<PRECISION>( testname + "-calbody.txt");
  const auto optpivotdata = read_optpivot_data<PRECISION>( testname + "-optpivot.txt");
  const auto empivotdata = read_empivot_data<PRECISION>(testname + "-empivot.txt");
  const auto calreadingsdata = read_calreadings_data<PRECISION>(testname + "-calreadings.txt");

  // Q4
  for(int i = 0; i < calreadingsdata.n_frames; i++){
    graph.register_transform("EMTRACKER", calbodydata.n_d_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.n_D_vals[i].transpose());
  }
  for(int i = 0; i < calreadingsdata.n_frames; i++){
    graph.register_transform("CALBODY_LOCAL", calbodydata.n_a_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.n_A_vals[i].transpose());
  }
  std::vector<Eigen::Matrix<PRECISION, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, 3, -1>>> c_expected_vals(calreadingsdata.n_frames);
  for(int i = 0; i < calreadingsdata.n_frames; i++) {
    auto tmp = graph.apply_direct_transform(
            "CALBODY_LOCAL", "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
            calbodydata.n_c_vals.transpose());
    tmp = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                                 "EMTRACKER", tmp);
    c_expected_vals[i] = tmp.transpose();
  }
  // Q5
  std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>> G_vals_transposed =
          std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>>(empivotdata.n_G_vals.size());
  for(size_t i = 0; i < empivotdata.n_G_vals.size(); i++)
    G_vals_transposed[i] = empivotdata.n_G_vals[i].transpose();
  auto empivot_results = pivot_calibration_routine(G_vals_transposed);
  // Q6
  graph.clear();
  std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>> optpivot_H_vals_transposed =
          std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>>(optpivotdata.n_H_vals.size());
  for(size_t i = 0; i < optpivotdata.n_H_vals.size(); i++) {
    graph.register_transform("EMTRACKER", calbodydata.n_d_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             optpivotdata.n_D_vals[i].transpose());
    optpivot_H_vals_transposed[i] = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                                                                 "EMTRACKER", optpivotdata.n_H_vals[i].transpose());
  }
  auto optpivot_results = pivot_calibration_routine(optpivot_H_vals_transposed);
  //write files
  outname = (outname.back() == '/') ? outname : outname + "/";
  write_output_file<PRECISION>(outname + "-output_us.txt", calreadingsdata.n_C, calreadingsdata.n_frames,
                            empivot_results.second, optpivot_results.second, c_expected_vals);
  return 0;
}