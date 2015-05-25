#ifndef __CUSTOMFUNCTIONS_H__
#define __CUSTOMFUNCTIONS_H__

#include <opencv2/core/core.hpp>
#include "MRFEnergy.h"
#include "Param.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#define WRITE_FLAG                  1
#define REPLACE_LAST_TUBE_AS_NULL   1
#define AUXILLARY_UNARY_VAL_MAX     1
#define AUXILLARY_BINARY_VAL        1
#define NUM_ITER                    10

using namespace std;
using namespace cv;

void loadPotentials(char* configFile, const int paramIter, vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary,
                    vector<TypeGeneral::REAL*>& skin_unary, vector<TypeGeneral::REAL*>& saliency_unary, vector<TypeGeneral::REAL*>& size_unary,
                    vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary, vector<vector<TypeGeneral::REAL*> >& phog_binary,
                    int& nodeNum, int& numLabels, string& writePath, vector<vector<vector<Rect_<int> > > >& gtTubes, vector<vector<vector<Rect_<int> > > >& detTubes);

void destroyPotentials(vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary, vector<TypeGeneral::REAL*>& skin_unary,
                       vector<TypeGeneral::REAL*>& saliency_unary,vector<TypeGeneral::REAL*>& size_unary,
                       vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary,
                       vector<vector<TypeGeneral::REAL*> >& phog_binary, const int& nodeNum);


/* video annotations related */
struct VideoAnno{
    int    start_frame;
    int    end_frame;

    string object_name;
    string seq_name;

    string bmf_file_path;
    string rgb_file_path;
    string dep_file_path;
    string opt_file_path;
    string spx_file_path;
    string jnt_rgb_un_file_path;
    string jnt_dep_file_path;
    string pos_3d_file_path;
    const  StructParam* param;

#if CAD_120_DATASET
    string subseq_name;
    int actor_idx;
#endif

};

struct ObjectAnno{
    vector<VideoAnno> vVideoAnno;
    int loadVideoAnno(StructParam* par);
    int size();
};

int loadTubes(const StructParam* param, const VideoAnno& anno, vector<vector<Rect> >& vvRect);
int loadGtTubes(const StructParam* param, const VideoAnno& anno, vector<vector<Rect> >& vvRect);
double calculate_solution_gt_overlap(const vector<vector<vector<Rect_<int> > > >& gtTubes, const vector<vector<vector<Rect_<int> > > >& detTubes, const vector<int>& solution_id);
double calculate_precision_recall(const vector<vector<vector<Rect_<int> > > >& refTubes, const vector<vector<vector<Rect_<int> > > >& tstTubes, const vector<int>& indices,
                                     double& precision, double& recall);

#endif
