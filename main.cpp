#include <chrono>
#include <iostream>
#include <filesystem>
#include "geometry/frame_lib.h"
#include "geometry/register.h"
#include "geometry/calibration.h"
#include "io/read_input_files.h"
#include "io/write_output_file.h"
#include "geometry/distortion_correction.h"

#define PRECISION double //all methods/classses/procedures work with both float and double precision.

#define BERNSTEIN_DEGREE 5 //degree of bernstein polynomials used in distortion correction

#define APPLY_DISTORTION_CORRECTION false

int main(int argc, char ** argv) {
  // ARGUMENT HANDLING

  //argument check
  if(argc != 3){
    std::cerr << "Malformed arguments. arg 1: dir of output_file (e.g OUTPUTS/ ). "
                 "arg 2: dir and name prefix of input files (e.g. input_files/pa1-debug-a ).";
    exit(1);
  }
  //file prefix to add e.g. "calbody.txt" to
  std::string testprefix = std::string(argv[2]);
  auto testname = testprefix.substr(testprefix.find_last_of('/') + 1,std::string::npos);
  auto outdir = std::string(argv[1]);
  outdir = (outdir.back() == '/') ? outdir : outdir + "/";
  auto outpath = std::filesystem::path(outdir);
  if (!exists(outpath)) {
    std::cerr << "Output Dir Not Found. Creating " + outdir << std::endl;
    std::filesystem::create_directory(outpath);
  }
  //are we in debug mode? if so we will output error analysis for PA2
  bool debug = testprefix.find("debug") != std::string::npos;

  // READING DATA
  //load data from files.
  FrameGraph<PRECISION> graph(Registration::PROCRUSTES);
  CALBODYDATA<PRECISION> calbodydata; EMDATA<PRECISION> empivotdata;
  CALREADINGSDATA<PRECISION> calreadingsdata;   EMDATA<PRECISION> emnavdata;
  CTFIDUCIALSDATA<PRECISION> ctfiducialsdata; EMDATA<PRECISION> emfiducialsdata;
  try {
    calbodydata = read_calbody_file<PRECISION>(testprefix + "-calbody.txt");
    empivotdata = read_em_data<PRECISION>(testprefix + "-empivot.txt");
    calreadingsdata = read_calreadings_data<PRECISION>(testprefix + "-calreadings.txt");
    ctfiducialsdata = read_ctfiducials_file<PRECISION>(testprefix + "-ct-fiducials.txt");
    emfiducialsdata = read_em_data<PRECISION>(testprefix + "-em-fiducialss.txt");
    emnavdata = read_em_data<PRECISION>(testprefix + "-EM-nav.txt");
  } catch (FileNotFound& e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  //setup some containers for use later.
  std::vector<Eigen::Matrix<PRECISION, 3, -1>, Eigen::aligned_allocator<Eigen::Matrix<PRECISION, 3, -1>>> c_expected_vals(calreadingsdata.n_frames);
  VECTOROFXD3MATRICES<PRECISION> distortion_corrected_G_vals_transposed = VECTOROFXD3MATRICES<PRECISION>(empivotdata.G_vals.size());
  if(APPLY_DISTORTION_CORRECTION)
    std::cout << "APPLYING DISTORTION CORRECTION WITH DEGREE " << BERNSTEIN_DEGREE << " POLYNOMIALS." << std::endl;


  auto begin = std::chrono::high_resolution_clock::now();
  // FRAME COMPUTATIONS, CALIBRATION, ETC.
  for(int i = 0; i < calreadingsdata.n_frames; i++)
    //for each frame of data, compute the (bidirectional) frame transform from the EM TRACKER to the OPTICAL TRACKER
    graph.register_transform("EMTRACKER", calbodydata.d_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.D_vals[i].transpose());

  for(int i = 0; i < calreadingsdata.n_frames; i++)
    //for each frame of data, compute the (bidirectional) frame transform from the CALIBRATION BODY to the OPTICAL TRACKER
    graph.register_transform("CALBODY_LOCAL", calbodydata.a_vals.transpose(),
                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                             calreadingsdata.A_vals[i].transpose());

  for(int i = 0; i < calreadingsdata.n_frames; i++) {
    //for each frame of data, transform the points on the CALIBRATION BODY to the EM TRACKER coordinate frame
    auto tmp = graph.apply_direct_transform(
            "CALBODY_LOCAL", "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
            calbodydata.c_vals.transpose());
    tmp = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
                                       "EMTRACKER", tmp);
    c_expected_vals[i] = tmp.transpose();
  }

  DISTORTIONCOEFFICIENTS<BERNSTEIN_DEGREE, PRECISION> coeffs;
  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_G_vals;
  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_fiducial_G_vals;
  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_em_nav_G_vals;

  if(APPLY_DISTORTION_CORRECTION) {
    //compute distortion from c_expected and apply distortion correction to pivot calibration data
    coeffs = compute_distortion_coeffs<BERNSTEIN_DEGREE, PRECISION>(calreadingsdata.C_vals, c_expected_vals);
    distortion_corrected_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs, empivotdata.G_vals);
    //transpose G vals
    for (size_t i = 0; i < empivotdata.G_vals.size(); i++)
      distortion_corrected_G_vals_transposed[i] = distortion_corrected_G_vals[i].transpose();
  } else {
    for (size_t i = 0; i < empivotdata.G_vals.size(); i++)
      distortion_corrected_G_vals_transposed[i] = empivotdata.G_vals[i].transpose();
  }

  //perform pivot calibration
  auto empivot_results = pivot_calibration_routine(distortion_corrected_G_vals_transposed);
  if(APPLY_DISTORTION_CORRECTION){
    //apply distortion correction to em fiducials data.
    distortion_corrected_fiducial_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,emfiducialsdata.G_vals);
  } else {
    distortion_corrected_fiducial_G_vals = emfiducialsdata.G_vals;
  }


  Eigen::Matrix<PRECISION, -1, 3> points_in_em(emfiducialsdata.G_vals.size(), 3);
  // the pointer tip location is with respect to frame_zero_local_coords so we need to align the fiducial EM tracker vals to frame_zero_local_coords
  //before we can compute the location of the tip.
  Eigen::Matrix<PRECISION, -1, 3> frame_zero_local_coords = distortion_corrected_G_vals_transposed[0].rowwise() - distortion_corrected_G_vals_transposed[0].colwise().mean();

  for (size_t i = 0; i < emfiducialsdata.G_vals.size(); i++){
    //compute the transform from the coordinates in which the tip location is known, to the current EM probe location (at a fiducial marker)
    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
                             "EM_PROBE_" + std::to_string(i), distortion_corrected_fiducial_G_vals[i].transpose());
    //then use the transform computed above to find the location of the tip for the current EM probe location (at a fiducial marker)
    points_in_em.template block<1, 3>(i, 0) = graph.apply_direct_transform("FRAME_0_PIVOTCAL", "EM_PROBE_" + std::to_string(i), empivot_results.first.transpose());
  }
  // apply procrustes to determine the transform from EM to CT_SCAN
  graph.register_transform("EM_TRACKER", points_in_em, "CT_SCAN", ctfiducialsdata.b_vals.transpose());
  // apply distortion correction to all values read from the EM, including the EM Nav vals.
  if (APPLY_DISTORTION_CORRECTION) {
    distortion_corrected_em_nav_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,
                                                                                                  emnavdata.G_vals);
  }  else {
    distortion_corrected_em_nav_G_vals = emnavdata.G_vals;
  }  //this will store the given EM Probe locations in CT coords
  Eigen::Matrix<PRECISION,-1, 3> ct_nav_G_vals(emnavdata.G_vals.size(), 3);
  for (size_t i = 0; i < emnavdata.G_vals.size(); i++){
    //compute the transform from the coordinates in which the tip location is known (taken from pivot calibration),
    // to the current EM probe location (at a nav loc)
    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
                             "CURRENT_EM_PROBE_" + std::to_string(i), distortion_corrected_em_nav_G_vals[i].transpose());
    //transform the location of the tip in EM coords to the ct coord system
    ct_nav_G_vals.row(i) = graph.apply_direct_transform("EM_TRACKER", "CT_SCAN",
                      //compute the location of the tip in EM coords.
          graph.apply_direct_transform("FRAME_0_PIVOTCAL", "CURRENT_EM_PROBE_" + std::to_string(i), empivot_results.first.transpose()));
  }

  std::cout << "Executed all computations in: " <<
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - begin).count() << " us" << std::endl;

  //DEBUG FILE ERROR ANALYSIS
  //write files
  if(debug) {
    try {
      const auto debug_data = read_debug_file_pa2<PRECISION>(testprefix + "-output2.txt");
      std::ofstream debug_error_file;
      debug_error_file.open(outdir + testname + "-error-analysis.txt");
      debug_error_file << "Mean Error for Nav Vals. : " << (ct_nav_G_vals - debug_data.ct_nav_vals).cwiseAbs().mean() << " mm" << std::endl;
      debug_error_file.close();
      std::cout << "Successfully wrote " << outdir << testname << "-error-analysis.txt" << std::endl;
    } catch (FileNotFound &e) {
      std::cout << e.what() << " Continuing without error check." << std::endl;
    }
  }

  write_output_file_pa2<PRECISION>(outdir + testname + "-output2-us.txt", emnavdata.G_vals.size(), ct_nav_G_vals);
  std::cout << "Successfully wrote " << outdir << testname << "-output2-us.txt" << std::endl;
  return 0;
}

//int main() {
//  //file prefix to add e.g. "calbody.txt" to
//  std::string testprefix = "input_data_pa2/pa2-debug-c"; //std::string(argv[2]);
//  auto testname = testprefix.substr(testprefix.find_last_of('/') + 1,std::string::npos);
//
//  // READING DATA
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
//  VECTOROFXD3MATRICES<PRECISION> distortion_corrected_G_vals_transposed = VECTOROFXD3MATRICES<PRECISION>(empivotdata.G_vals.size());
//
//  std::cout << "APPLYING DISTORTION CORRECTION WITH DEGREE " << BERNSTEIN_DEGREE << " POLYNOMIALS." << std::endl;
//
//
//  auto begin = std::chrono::high_resolution_clock::now();
//  // FRAME COMPUTATIONS, CALIBRATION, ETC.
//  for(int i = 0; i < calreadingsdata.n_frames; i++)
//    //for each frame of data, compute the (bidirectional) frame transform from the EM TRACKER to the OPTICAL TRACKER
//    graph.register_transform("EMTRACKER", calbodydata.d_vals.transpose(),
//                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                             calreadingsdata.D_vals[i].transpose());
//
//  for(int i = 0; i < calreadingsdata.n_frames; i++)
//    //for each frame of data, compute the (bidirectional) frame transform from the CALIBRATION BODY to the OPTICAL TRACKER
//    graph.register_transform("CALBODY_LOCAL", calbodydata.a_vals.transpose(),
//                             "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                             calreadingsdata.A_vals[i].transpose());
//
//  for(int i = 0; i < calreadingsdata.n_frames; i++) {
//    //for each frame of data, transform the points on the CALIBRATION BODY to the EM TRACKER coordinate frame
//    auto tmp = graph.apply_direct_transform(
//            "CALBODY_LOCAL", "OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//            calbodydata.c_vals.transpose());
//    tmp = graph.apply_direct_transform("OPTICAL_TRACKER_FRAME_" + std::to_string(i),
//                                       "EMTRACKER", tmp);
//    c_expected_vals[i] = tmp.transpose();
//  }
//
//  //compute distortion from c_expected and apply distortion correction to pivot calibration data
//  DISTORTIONCOEFFICIENTS<BERNSTEIN_DEGREE, PRECISION> coeffs = compute_distortion_coeffs<BERNSTEIN_DEGREE, PRECISION>(calreadingsdata.C_vals, c_expected_vals);
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs, empivotdata.G_vals);
//
//  //transpose G vals
//  for(size_t i = 0; i < empivotdata.G_vals.size(); i++)
//    distortion_corrected_G_vals_transposed[i] = distortion_corrected_G_vals[i].transpose();
//
//  //perform pivot calibration
//  auto empivot_results = pivot_calibration_routine(distortion_corrected_G_vals_transposed);
//
//  //apply distortion correction to em fiducials data.
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_fiducial_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,emfiducialsdata.G_vals);
//
//  Eigen::Matrix<PRECISION, -1, 3> points_in_em(emfiducialsdata.G_vals.size(), 3);
//  // the pointer tip location is with respect to frame_zero_local_coords so we need to align the fiducial EM tracker vals to frame_zero_local_coords
//  //before we can compute the location of the tip.
//  Eigen::Matrix<PRECISION, -1, 3> frame_zero_local_coords = distortion_corrected_G_vals_transposed[0].rowwise() - distortion_corrected_G_vals_transposed[0].colwise().mean();
//
//  for (size_t i = 0; i < emfiducialsdata.G_vals.size(); i++){
//    //compute the transform from the coordinates in which the tip location is known, to the current EM probe location (at a fiducial marker)
//    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
//                             "EM_PROBE_" + std::to_string(i), distortion_corrected_fiducial_G_vals[i].transpose());
//    //then use the transform computed above to find the location of the tip for the current EM probe location (at a fiducial marker)
//    points_in_em.template block<1, 3>(i, 0) = graph.apply_direct_transform("FRAME_0_PIVOTCAL", "EM_PROBE_" + std::to_string(i), empivot_results.first.transpose());
//  }
//  // apply procrustes to determine the transform from EM to CT_SCAN
//  graph.register_transform("EM_TRACKER", points_in_em, "CT_SCAN", ctfiducialsdata.b_vals.transpose());
//  // apply distortion correction to all values read from the EM, including the EM Nav vals.
//  VECTOROF3XDMATRICES<PRECISION> distortion_corrected_em_nav_G_vals = apply_distortion_correction<BERNSTEIN_DEGREE, PRECISION>(coeffs,emnavdata.G_vals);
//  //this will store the given EM Probe locations in CT coords
//  Eigen::Matrix<PRECISION,-1, 3> ct_nav_G_vals(emnavdata.G_vals.size(), 3);
//  for (size_t i = 0; i < emnavdata.G_vals.size(); i++){
//    //compute the transform from the coordinates in which the tip location is known (taken from pivot calibration),
//    // to the current EM probe location (at a nav loc)
//    graph.register_transform("FRAME_0_PIVOTCAL", frame_zero_local_coords,
//                             "CURRENT_EM_PROBE_" + std::to_string(i), distortion_corrected_em_nav_G_vals[i].transpose());
//    //transform the location of the tip in EM coords to the ct coord system
//    ct_nav_G_vals.row(i) = graph.apply_direct_transform("EM_TRACKER", "CT_SCAN",
//            //compute the location of the tip in EM coords.
//                                                        graph.apply_direct_transform("FRAME_0_PIVOTCAL", "CURRENT_EM_PROBE_" + std::to_string(i), empivot_results.first.transpose()));
//  }
//
//  std::cout << "Executed all computations in: " <<
//            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - begin).count() << " us" << std::endl;
//
//  //DEBUG FILE ERROR ANALYSIS
//  //write files
//  const auto debug_data = read_debug_file_pa2<PRECISION>(testprefix + "-output2.txt");
//  std::cout << "Mean Error for Nav Vals. : " << (ct_nav_G_vals - debug_data.ct_nav_vals).cwiseAbs().mean() << " mm" << std::endl;
//  return 0;
//}