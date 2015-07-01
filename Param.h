/*
// Author: Abhilash Srikantha (MPI Tuebingen, UniBonn)
// Email: abhilash.srikantha@tue.mpg.de
*/

#ifndef ParamH
#define ParamH

#include <string>
#include <vector>
#include <fstream>

#define ETHZ_ACTION_DATASET         0
#define CAD_120_DATASET             1
#define MPII_COOKING_DATASET        0

struct img_size {
    int h;
    int w;
};

// program parameters
struct StructParam {

    std::string data_path;          // root path of the data
    std::string tube_file_path;     // path for saving tubes
    std::string pot_file_path;      // path for storing potentials
    std::string surf_bow_file_path; // path for saving surf bow file
    std::string action_list_file;   // location of file that contains list of actors and corresponding sequences to consider

    img_size rgb_image_size;        // size of rgb image
    img_size dep_image_size;        // size of depth image
    img_size size_prior_mu;         // size prior mean
    img_size size_prior_var;        // size prior variance

    std::string opt_flow_idn;       // identifier for optical flow
    std::string superpixel_idn;     // identifier for superpixel images
    std::string rgb_image_idn;      // identifier for rgb image
    std::string dep_image_idn;      // identifier for depth image
    std::string joints_rgb_idn;     // path for pose-w.r.t-rgb_un-image file
    std::string joints_dep_idn;     // path for pose-w.r.t-dep-image file
    std::string rgb_calib_path;     // path for rgb calibration file

    bool loadConfigVideoTrain(const char* filename);
    int  numberOfActions();

    int num_rand_frames;            // number of randomly selected frames
    int num_rand_points;            // number of sample points per frame
    int perc_dist_hand_head;        // distance around hand to sample
    int num_frames_qual_eval;       // number of frames for tube quality evaluation
    int num_frames_surf_bow;        // number of images for SURF BoW dictionary

#if CAD_120_DATASET || MPII_COOKING_DATASET
    std::string gt_path;            // path to ground truth
    std::string bmf_file_path;
    std::string jnt_path;
#endif

#if MPII_COOKING_DATASET
    int tube_len;
#endif
};

template<typename T1, typename T2>
std::ostream& operator << (std::ostream& out, const std::pair<T1,T2>& inpair){
    out << "[" << inpair.first << "," << inpair.second <<"]" << "\t";
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last)
            out << ", ";
    }
    out << "]";
    return out;
}

#endif
