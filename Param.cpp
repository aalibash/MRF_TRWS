/*
// Author: Abhilash Srikantha (MPI Tuebingen, UniBonn)
// Email: abhilash.srikantha@tue.mpg.de
*/

#include "Param.h"

#include <libconfig.h++>
#include <iostream>

using namespace libconfig;
using namespace std;

// parses parameters for detection
bool StructParam::loadConfigVideoTrain(const char* filename) {

    Config configFile;

    //try to read configuration file
    try {
        configFile.readFile(filename);
    }
    catch(const FileIOException &fioex) {
        cerr << "Could not read config file " << filename << endl;
        return false;
    }
    catch(ParseException &pex) {
        cerr << "Parse error at " << filename << ":" << pex.getLine()
        << " - " << pex.getError() << endl;
        return false;
    }

    // look up parameter settings / values
    try{

        // root path of the data
        data_path = (const char*)configFile.lookup("data_path");
        // location of file that contains list of actors and corresponding sequences to consider
        action_list_file = (const char*)configFile.lookup("action_list_file");
        // paht for storing potentials
        pot_file_path = (const char*)configFile.lookup("pot_file_path");
        // size of rgb image
        rgb_image_size.w = configFile.lookup("rgb_image_size")[0];
        rgb_image_size.h = configFile.lookup("rgb_image_size")[1];
        // size of depth image
        dep_image_size.w = configFile.lookup("dep_image_size")[0];
        dep_image_size.h = configFile.lookup("dep_image_size")[1];
        // identifier for optical flow
        opt_flow_idn = (const char*)configFile.lookup("opt_flow_idn");
        // identifier for superpixel images
        superpixel_idn = (const char*)configFile.lookup("superpixel_idn");
        // identifier for rgb image
        rgb_image_idn = (const char*)configFile.lookup("rgb_image_idn");
        // identifier for depth image
        dep_image_idn = (const char*)configFile.lookup("dep_image_idn");
        // path for pose-w.r.t-rgb-image file
        joints_rgb_idn = (const char*)configFile.lookup("joints_rgb_idn");
        joints_dep_idn = (const char*)configFile.lookup("joints_dep_idn");

#if ETHZ_ACTION_DATASET
        // RGB camera calibration file
        rgb_calib_path = (const char*)configFile.lookup("rgb_calib_path");
#endif

        // size priors
        size_prior_mu.w  = 40;
        size_prior_mu.h  = 20;
        size_prior_var.w = 50*50;
        size_prior_var.h = 50*50;

        // number of randomly selected frames
        num_rand_frames = configFile.lookup("num_rand_frames");
        // number of sample points per frame
        num_rand_points = configFile.lookup("num_rand_points");
        // distance around hand to sample
        perc_dist_hand_head = configFile.lookup("perc_dist_hand_head");
        // number of frames for tube quality evaluation
        num_frames_qual_eval = configFile.lookup("num_frames_qual_eval");
        // number of images for SURF BoW dictionary
        num_frames_surf_bow = configFile.lookup("num_frames_surf_bow");
        // surf bow file path
        surf_bow_file_path = (const char*)configFile.lookup("surf_bow_file_path");

        // path to save the tubes
        tube_file_path = (const char*)configFile.lookup("tube_file_path");

#if CAD_120_DATASET
        gt_path = (const char*)configFile.lookup("gt_path");
        bmf_file_path = (const char*)configFile.lookup("bmf_file_path");
        jnt_path = (const char*)configFile.lookup("joint_file_path");
#endif

#if MPII_COOKING_DATASET
        gt_path = (const char*)configFile.lookup("gt_path");
        tube_len = configFile.lookup("tube_length");
#endif

//        cout<<"debug : "<<endl;
//        cout<<"data_path : "<<data_path<<endl;

    }
    catch(const SettingNotFoundException &nfex) {
        cerr << "Not found in configuration file!" << endl;
        return false;
    }

    return true;

}
