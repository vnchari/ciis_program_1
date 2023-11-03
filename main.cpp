#include <chrono>
#include <iostream>
#include <filesystem>
//#include "geometry/frame_lib.h"
//#include "geometry/register.h"
//#include "geometry/calibration.h"
#include "io/read_input_files.h"
//#include "io/write_output_file.h"
#include "geometry/distortion_correction.h"

#define PRECISION double //all methods/classses/procedures work with both float and double precision.

#define BERNSTEIN_DEGREE 5 //all methods/classses/procedures work with both float and double precision.


int main(int argc, char ** argv) {
  std::cout << "BERNSTEIN DEGREE: " << BERNSTEIN_DEGREE << std::endl;
//  if(argc != 3){
//    std::cerr << "Malformed arguments. arg 1: dir of output_file (e.g OUTPUTS/ ). "
//                 "arg 2: dir and name prefix of input files (e.g. input_files/pa1-debug-a ).";
//    exit(1);
//  }
  std::string testprefix = "input_data_pa2/pa2-debug-c"; //std::string(argv[2]);
  auto testname = testprefix.substr(testprefix.find_last_of('/') + 1,std::string::npos);
//  auto outdir = std::string(argv[1]);
//  outdir = (outdir.back() == '/') ? outdir : outdir + "/";
//  auto outpath = std::filesystem::path(outdir);
//  if (!exists(outpath)) {
//    std::cerr << "Output Dir Not Found. Creating " + outdir << std::endl;
//    std::filesystem::create_directory(outpath);
//  }
  const auto debug_data = read_debug_file_pa2<PRECISION>(testprefix + "-output2.txt");
  std::cout << debug_data.ct_nav_vals << std::endl;
//  bool debug = testprefix.find("debug") != std::string::npos;
//  //load data from files.
//  FrameGraph<PRECISION> graph(Registration::PROCRUSTES);
//  CALBODYDATA<PRECISION> calbodydata; EMDATA<PRECISION> empivotdata;
//  CALREADINGSDATA<PRECISION> calreadingsdata;   EMDATA<PRECISION> emnavdata;
//  CTFIDUCIALSDATA<PRECISION> ctfiducialsdata; EMDATA<PRECISION> emfiducialsdata;
//  try {
//    calbodydata = read_calbody_file<PRECISION>(testprefix + "-calbody.txt");
//    empivotdata = read_em_data<PRECISION>(testprefix + "-empivot.txt");
//    calreadingsdata = read_calreadings_data<PRECISION>(testprefix + "-calreadings.txt");
//    ctfiducialsdata = read_ctfiducials_file<PRECISION>(testprefix + "-ct-fiducials.txt");
//    emfiducialsdata = read_em_data<PRECISION>(testprefix + "-em-fiducialss.txt");
//    emnavdata = read_em_data<PRECISION>(testprefix + "-EM-nav.txt");
//  } catch (FileNotFound& e) {
//    std::cerr << e.what() << std::endl;
//    exit(1);
//  }
//  //setup some containers for use later.
//  std::vector<Eigen::Matrix<PRECISION, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, 3, -1>>> c_expected_vals(calreadingsdata.n_frames);
//  VECTOROFXD3MATRICES<PRECISION> distortion_corrected_G_vals_transposed = VECTOROFXD3MATRICES<PRECISION>(empivotdata.n_G_vals.size());
//
//  auto begin = std::chrono::high_resolution_clock::now();
//
//  for(int i = 0; i < calreadingsdata.n_frames; i++)
//    graph.register_transform("EMTRACKER", calbodydata.n_d_vals.transpose(),
//                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                             calreadingsdata.n_D_vals[i].transpose());
//
//  for(int i = 0; i < calreadingsdata.n_frames; i++)
//    graph.register_transform("CALBODY_LOCAL", calbodydata.n_a_vals.transpose(),
//                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                             calreadingsdata.n_A_vals[i].transpose());
//
//  for(int i = 0; i < calreadingsdata.n_frames; i++) {
//    auto tmp = graph.apply_direct_transform(
//            "CALBODY_LOCAL", "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//            calbodydata.n_c_vals.transpose());
//    tmp = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                                       "EMTRACKER", tmp);
//    c_expected_vals[i] = tmp.transpose();
//  }
//
//  DISTORTIONCOEFFICIENTS<BERNSTEIN_DEGREE, PRECISION> coeffs = compute_distortion_coeffs<BERNSTEIN_DEGREE, PRECISION>(calreadingsdata.n_C_vals,
//                                                                       c_expected_vals);
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs, empivotdata.n_G_vals);
//
//  for(size_t i = 0; i < empivotdata.n_G_vals.size(); i++)
//    distortion_corrected_G_vals_transposed[i] = distortion_corrected_G_vals[i].transpose();
//
//  auto empivot_results = pivot_calibration_routine(distortion_corrected_G_vals_transposed);
//
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_fiducial_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,emfiducialsdata.n_G_vals);
//
//  Eigen::Matrix<PRECISION, -1, 3> points_in_em(emfiducialsdata.n_G_vals.size(), 3);
//  Eigen::Matrix<PRECISION, -1, 3> frame_zero_local_coords = distortion_corrected_G_vals_transposed[0].rowwise() - distortion_corrected_G_vals_transposed[0].colwise().mean();
//  for (size_t i = 0; i < emfiducialsdata.n_G_vals.size(); i++){
//    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
//                             "EM_PROBE_" + std::to_string(i), distortion_corrected_fiducial_G_vals[i].transpose());
//    points_in_em.template block<1, 3>(i, 0) = graph.apply_direct_transform("FRAME_0_PIVOTCAL", "EM_PROBE_" + std::to_string(i), empivot_results.first.transpose());
//  }
//
//  graph.register_transform("EM_TRACKER", points_in_em, "CT_SCAN", ctfiducialsdata.n_b_vals.transpose());
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_em_nav_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,emnavdata.n_G_vals);
//  Eigen::Matrix<PRECISION,-1, 3> em_nav_G_vals(emnavdata.n_G_vals.size(), 3);
//  for (size_t i = 0; i < emnavdata.n_G_vals.size(); i++){
//    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
//                             "CURRENT_EM_PROBE_" + std::to_string(i), distortion_corrected_em_nav_G_vals[i].transpose());
//    em_nav_G_vals.row(i) = graph.apply_direct_transform("EM_TRACKER", "CT_SCAN",
//          graph.apply_direct_transform("FRAME_0_PIVOTCAL", "CURRENT_EM_PROBE_" + std::to_string(i), empivot_results.first.transpose()));
//  }
//
//  std::cout << "Executed all computations in: " <<
//            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - begin).count() << " us" << std::endl;
////  //write files
////  if(debug) {
////    try {
////      const auto debug_data = read_debug_file<PRECISION>(testprefix + "-output1.txt");
////      std::ofstream debug_error_file;
////      debug_error_file.open(outdir + testname + "-error-analysis.txt");
////      debug_error_file << "Mean Error for EM Pivot Post Pos. : " << (empivot_results.second - debug_data.em_pivot_post_pos).cwiseAbs().mean() << " mm" << std::endl;
////      debug_error_file << "Mean Error for Optical Pivot Post Pos. : " << (optpivot_results.second - debug_data.opt_pivot_post_pos).cwiseAbs().mean() << " mm" << std::endl;
////      PRECISION error = 0.0;
////      for(size_t i = 0; i <  debug_data.n_frames; i++)
////        error += (debug_data.C_vals[i] - c_expected_vals[i]).cwiseAbs().mean();
////      debug_error_file << "Mean Error in C: " << error / (PRECISION) debug_data.n_frames << " mm" << std::endl;
////      debug_error_file.close();
////      std::cout << "Successfully wrote " << outdir << testname << "-error-analysis.txt" << std::endl;
////    } catch (FileNotFound &e) {
////      std::cout << e.what() << " Continuing without error check." << std::endl;
////    }
////  }
////  write_output_file<PRECISION>(outdir + testname + "-output-us.txt", calreadingsdata.n_C, calreadingsdata.n_frames,
////                               empivot_results.second, optpivot_results.second, c_expected_vals);
////  std::cout << "Successfully wrote " << outdir << testname << "-output-us.txt" << std::endl;
//  return 0;
}