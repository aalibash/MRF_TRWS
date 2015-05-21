#include "customFunctions.h"

void loadPotentials(char* configFile, vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary,
                    vector<TypeGeneral::REAL*>& skin_unary, vector<TypeGeneral::REAL*>& saliency_unary, vector<TypeGeneral::REAL*>& size_unary,
                    vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary, vector<vector<TypeGeneral::REAL*> >& phog_binary,
                    int& nodeNum, int& numLabels, string& writePath){

    string handpath_unary_path, skin_unary_path, saliency_unary_path, size_unary_path, loc_binary_path, phog_binary_path, filtered_tube_ind_path;
    string gt_tubes_path, det_tubes_path;
    vector<vector<Rect_<int> > > gt_tubes, det_tubes;
    TypeGeneral::REAL temp;
    int row, col;
    float maxWeight=1000;

    // read the potential file names
    ifstream iStream(configFile);
    iStream >> gt_tubes_path;
    iStream >> det_tubes_path;
    iStream >> handpath_unary_path;
    iStream >> skin_unary_path;
    iStream >> saliency_unary_path;
    iStream >> size_unary_path;
    iStream >> filtered_tube_ind_path;
    iStream >> loc_binary_path;
    iStream >> phog_binary_path;
    iStream >> writePath;
    iStream.close();

    // read the tubes that must be selected
    vector<vector<bool> > isChosen;
    iStream.open(filtered_tube_ind_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    isChosen.resize(nodeNum);
    int tempIdx;
    for(unsigned int rowidx=0; rowidx<isChosen.size(); ++rowidx){
        isChosen[rowidx].resize(numLabels,0);
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> tempIdx;
            isChosen[rowidx][tempIdx-1] = 1;
        }
    }
    iStream.close();

    // read handpath_unary potentials
    iStream.open(handpath_unary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    handpath_unary.resize(nodeNum);
    unary.resize(nodeNum);
    for(unsigned int rowidx=0; rowidx<handpath_unary.size(); ++rowidx){
        handpath_unary[rowidx] = new TypeGeneral::REAL [numLabels];
        unary[rowidx] = new TypeGeneral::REAL [numLabels];
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> handpath_unary[rowidx][idx];
            if(handpath_unary[rowidx][idx]>1000){
                handpath_unary[rowidx][idx]=1000;
            }
            handpath_unary[rowidx][idx]/=1000;
//            if(!isChosen[rowidx][idx]) handpath_unary[rowidx][idx]=maxWeight;
        }
    }
    iStream.close();

    // read skin_unary potentials
    iStream.open(skin_unary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    skin_unary.resize(nodeNum);
    for(unsigned int rowidx=0; rowidx<skin_unary.size(); ++rowidx){
        skin_unary[rowidx] = new TypeGeneral::REAL [numLabels];
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> skin_unary[rowidx][idx];
//            if(!isChosen[rowidx][idx]) skin_unary[rowidx][idx]=maxWeight;
        }
    }
    iStream.close();

    // read size_unary potentials
    iStream.open(size_unary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    size_unary.resize(nodeNum);
    for(unsigned int rowidx=0; rowidx<size_unary.size(); ++rowidx){
        size_unary[rowidx] = new TypeGeneral::REAL [numLabels];
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> size_unary[rowidx][idx];
            size_unary[rowidx][idx] = 1-size_unary[rowidx][idx];
//            if(!isChosen[rowidx][idx]) size_unary[rowidx][idx]=maxWeight;
        }
    }
    iStream.close();

    // read saliency_unary potentials
    iStream.open(saliency_unary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    saliency_unary.resize(nodeNum);
    for(unsigned int rowidx=0; rowidx<saliency_unary.size(); ++rowidx){
        saliency_unary[rowidx] = new TypeGeneral::REAL [numLabels];
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> saliency_unary[rowidx][idx];
            saliency_unary[rowidx][idx] = 1-saliency_unary[rowidx][idx];
//            if(!isChosen[rowidx][idx]) saliency_unary[rowidx][idx]=maxWeight;
        }
    }
    iStream.close();

    // read loc_binary potentials
    iStream.open(loc_binary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    loc_binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            iStream >> row;
            iStream >> col;
            loc_binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
            for(int idx=0; idx<numLabels*numLabels; ++idx){
                iStream >> loc_binary[row][col][idx];
            }
        }
    }
    iStream.close();

    // read phog_binary potentials
    iStream.open(phog_binary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    phog_binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
    binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            iStream >> row;
            iStream >> col;
            phog_binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
            binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
            for(int idx=0; idx<numLabels*numLabels; ++idx){
                iStream >> phog_binary[row][col][idx];
            }
        }
    }
    iStream.close();

#if REPLACE_LAST_TUBE_AS_NULL
    for(unsigned int rowidx=0; rowidx<saliency_unary.size(); ++rowidx){
        handpath_unary[rowidx][numLabels-1] = MININUM_UNARY_VAL;
        skin_unary[rowidx][numLabels-1] = MININUM_UNARY_VAL;
        size_unary[rowidx][numLabels-1] = MININUM_UNARY_VAL;
        saliency_unary[rowidx][numLabels-1] = MININUM_UNARY_VAL;
    }
    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            for(int idx=(numLabels-1)*numLabels; idx<numLabels*numLabels; ++idx){
                phog_binary[rowidx][colidx][idx]=MININUM_BINARY_VAL;
                loc_binary[rowidx][colidx][idx]=MININUM_BINARY_VAL;
            }
            for(int idx=numLabels-1; idx<numLabels*numLabels; idx+=numLabels){
                phog_binary[rowidx][colidx][idx]=MININUM_BINARY_VAL;
                loc_binary[rowidx][colidx][idx]=MININUM_BINARY_VAL;
            }
        }
    }
#endif

#if WRITE_FLAG
    // clear the write file file
    ofstream oStream(writePath.c_str());
    oStream.close();
#endif

}

void destroyPotentials(vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary, vector<TypeGeneral::REAL*>& skin_unary,
                       vector<TypeGeneral::REAL*>& saliency_unary,vector<TypeGeneral::REAL*>& size_unary,
                       vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary,
                       vector<vector<TypeGeneral::REAL*> >& phog_binary, const int& nodeNum){

    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        delete [] unary[rowidx];
        delete [] handpath_unary[rowidx];
        delete [] skin_unary[rowidx];
        delete [] saliency_unary[rowidx];
        delete [] size_unary[rowidx];
    }

    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            delete [] binary[rowidx][colidx];
            delete [] loc_binary[rowidx][colidx];
            delete [] phog_binary[rowidx][colidx];
        }
    }
}
