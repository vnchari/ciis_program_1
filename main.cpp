#include <chrono>
#include <iostream>
#include <filesystem>
#include "geometry/frame_lib.h"
#include "geometry/register.h"
#include "geometry/calibration.h"
#include "io/read_input_files.h"
#include "io/write_output_file.h"

#define PRECISION double //all methods/classses/procedures work with both float and double precision.

int main(int argc, char ** argv) {
  if(argc != 3){
    std::cerr << "Malformed arguments. arg 1: dir of output_file (e.g OUTPUTS/ ). "
                 "arg 2: dir and name prefix of input files (e.g. input_files/pa1-debug-a ).";
    exit(1);
  }
  auto testprefix = std::string(argv[2]);
  auto testname = testprefix.substr(testprefix.find_last_of('/') + 1,std::string::npos);
  auto outdir = std::string(argv[1]);
  outdir = (outdir.back() == '/') ? outdir : outdir + "/";
  auto outpath = std::filesystem::path(outdir);
  if (!exists(outpath)) {
    std::cerr << "Output Dir Not Found. Creating " + outdir << std::endl;
    std::filesystem::create_directory(outpath);
  }

  bool debug = testprefix.find("debug") != std::string::npos;
  //load data from files.
  FrameGraph<PRECISION> graph(Registration::PROCRUSTES);
  CALBODYDATA<PRECISION> calbodydata; OPTPIVOTDATA<PRECISION> optpivotdata;
  EMPIVOTDATA<PRECISION> empivotdata; CALREADINGSDATA<PRECISION> calreadingsdata;
  try {
    calbodydata = read_calbody_file<PRECISION>(testprefix + "-calbody.txt");
    optpivotdata = read_optpivot_data<PRECISION>(testprefix + "-optpivot.txt");
    empivotdata = read_empivot_data<PRECISION>(testprefix + "-empivot.txt");
    calreadingsdata = read_calreadings_data<PRECISION>(testprefix + "-calreadings.txt");
  } catch (FileNotFound& e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  //setup some containers for use later.
  std::vector<Eigen::Matrix<PRECISION, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, 3, -1>>> c_expected_vals(calreadingsdata.n_frames);
  std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>> G_vals_transposed =
          std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>>(empivotdata.n_G_vals.size());
  std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>> optpivot_H_vals_transposed =
          std::vector<Eigen::Matrix<PRECISION, -1, 3>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, -1, 3>>>(optpivotdata.n_H_vals.size());
  for(size_t i = 0; i < empivotdata.n_G_vals.size(); i++)
    G_vals_transposed[i] = empivotdata.n_G_vals[i].transpose();
  auto begin = std::chrono::high_resolution_clock::now();
  // Q4
  for(int i = 0; i < calreadingsdata.n_frames; i++)
    graph.register_transform("EMTRACKER", calbodydata.n_d_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.n_D_vals[i].transpose());

  for(int i = 0; i < calreadingsdata.n_frames; i++)
    graph.register_transform("CALBODY_LOCAL", calbodydata.n_a_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.n_A_vals[i].transpose());

  for(int i = 0; i < calreadingsdata.n_frames; i++) {
    auto tmp = graph.apply_direct_transform(
            "CALBODY_LOCAL", "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
            calbodydata.n_c_vals.transpose());
    tmp = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                                 "EMTRACKER", tmp);
    c_expected_vals[i] = tmp.transpose();
  }
  // Q5
  auto empivot_results = pivot_calibration_routine(G_vals_transposed);
  // Q6
  graph.clear();
  for(size_t i = 0; i < optpivotdata.n_H_vals.size(); i++) {
    graph.register_transform("EMTRACKER", calbodydata.n_d_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             optpivotdata.n_D_vals[i].transpose());
    optpivot_H_vals_transposed[i] = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                                                                 "EMTRACKER", optpivotdata.n_H_vals[i].transpose());
  }
  auto optpivot_results = pivot_calibration_routine(optpivot_H_vals_transposed);
  std::cout << "Executed all computations in: " <<
      std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - begin).count() << " us" << std::endl;
  //write files
  if(debug) {
    try {
      const auto debug_data = read_debug_file<PRECISION>(testprefix + "-output1.txt");
      std::ofstream debug_error_file;
      debug_error_file.open(outdir + testname + "-error-analysis.txt");
      debug_error_file << "Mean Error for EM Pivot Post Pos. : " << (empivot_results.second - debug_data.em_pivot_post_pos).cwiseAbs().mean() << " mm" << std::endl;
      debug_error_file << "Mean Error for Optical Pivot Post Pos. : " << (optpivot_results.second - debug_data.opt_pivot_post_pos).cwiseAbs().mean() << " mm" << std::endl;
      PRECISION error = 0.0;
      for(size_t i = 0; i <  debug_data.n_frames; i++)
        error += (debug_data.C_vals[i] - c_expected_vals[i]).cwiseAbs().mean();
      debug_error_file << "Mean Error in C: " << error / (PRECISION) debug_data.n_frames << " mm" << std::endl;
      debug_error_file.close();
      std::cout << "Successfully wrote " << outdir << testname << "-error-analysis.txt" << std::endl;
    } catch (FileNotFound &e) {
      std::cout << e.what() << " Continuing without error check." << std::endl;
    }
  }
  write_output_file<PRECISION>(outdir + testname + "-output-us.txt", calreadingsdata.n_C, calreadingsdata.n_frames,
                               empivot_results.second, optpivot_results.second, c_expected_vals);
  std::cout << "Successfully wrote " << outdir << testname << "-output-us.txt" << std::endl;
  return 0;
}