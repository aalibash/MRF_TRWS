#ifndef __CUSTOMFUNCTIONS_H__
#define __CUSTOMFUNCTIONS_H__

#include <opencv2/core/core.hpp>"
#include "MRFEnergy.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#define WRITE_FLAG                1
#define REPLACE_LAST_TUBE_AS_NULL 1
#define MININUM_UNARY_VAL         0.5
#define MININUM_BINARY_VAL        1

using namespace std;
using namespace cv;

void loadPotentials(char* configFile, vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary,
                    vector<TypeGeneral::REAL*>& skin_unary, vector<TypeGeneral::REAL*>& saliency_unary, vector<TypeGeneral::REAL*>& size_unary,
                    vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary, vector<vector<TypeGeneral::REAL*> >& phog_binary,
                    int& nodeNum, int& numLabels, string& writePath);

void destroyPotentials(vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary, vector<TypeGeneral::REAL*>& skin_unary,
                       vector<TypeGeneral::REAL*>& saliency_unary,vector<TypeGeneral::REAL*>& size_unary,
                       vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary,
                       vector<vector<TypeGeneral::REAL*> >& phog_binary, const int& nodeNum);

#endif
