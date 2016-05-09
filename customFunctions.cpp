#include "customFunctions.h"
#include <libconfig.h++>

using namespace libconfig;

void loadPotentials(char* configFile, const int paramIter, vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& handpath_unary,
                    vector<TypeGeneral::REAL*>& skin_unary, vector<TypeGeneral::REAL*>& saliency_unary, vector<TypeGeneral::REAL*>& size_unary, vector<TypeGeneral::REAL*>& objness_unary,
                    vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary, vector<vector<TypeGeneral::REAL*> >& phog_binary,
                    int& nodeNum, int& numLabels, string& writePath, vector<vector<vector<Rect_<int> > > >& gtTubes, vector<vector<vector<Rect_<int> > > >& detTubes,
                    vector<vector<pair<int,int> > >& tubeBounds){

    string handpath_unary_path, skin_unary_path, saliency_unary_path, size_unary_path, loc_binary_path, phog_binary_path, filtered_tube_ind_path, objness_unary_path;
    string configFilePath;
    char buffer[1024];
    TypeGeneral::REAL temp;
    int row, col;
    float maxWeight=1000;

    // read the potential file names
    ifstream iStream(configFile);
    iStream >> configFilePath;
    iStream >> handpath_unary_path;
    iStream >> skin_unary_path;
    iStream >> saliency_unary_path;
    iStream >> size_unary_path;
    iStream >> filtered_tube_ind_path;
    iStream >> loc_binary_path;
    iStream >> phog_binary_path;
    iStream >> writePath;

#if INCLUDE_OBJNESS_POTENTIAL
    iStream >> objness_unary_path;
#endif
    iStream.close();

    // update save path
    stringstream ss;
    ss << paramIter;
    sprintf(buffer,(writePath).c_str(),(ss.str()).c_str());
    writePath = string(buffer);

    // load the ground truth and detected tubes
    StructParam param;
    param.loadConfigVideoTrain(configFilePath.c_str());

    ObjectAnno oAnno;
    oAnno.loadVideoAnno(&param);

    gtTubes.resize(oAnno.size());
    detTubes.resize(oAnno.size());
    tubeBounds.resize(oAnno.size());
    for(unsigned int idx=0; idx<oAnno.size(); ++idx){
        loadGtTubes(&param, oAnno.vVideoAnno[idx], gtTubes[idx]);
        loadTubes(&param, oAnno.vVideoAnno[idx], detTubes[idx]);
#if MPII_COOKING_DATASET
        loadTubeBounds(&param,oAnno.vVideoAnno[idx],tubeBounds[idx]);
#endif
    }

#if ETHZ_ACTION_DATASET || CAD_120_DATASET
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
#endif

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

    // read objness_unary potentials
    iStream.open(objness_unary_path.c_str());
    iStream >> nodeNum;
    iStream >> numLabels;
    objness_unary.resize(nodeNum);
    for(unsigned int rowidx=0; rowidx<objness_unary.size(); ++rowidx){
        objness_unary[rowidx] = new TypeGeneral::REAL [numLabels];
        for(int idx=0; idx<numLabels; ++idx){
            iStream >> objness_unary[rowidx][idx];
            if(objness_unary[rowidx][idx]>10){
                objness_unary[rowidx][idx]=10;
            }
            objness_unary[rowidx][idx]/=10;
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
        handpath_unary[rowidx][numLabels-1] = (double)paramIter*AUXILLARY_UNARY_VAL_MAX/NUM_ITER;
        skin_unary[rowidx][numLabels-1] = (double)paramIter*AUXILLARY_UNARY_VAL_MAX/NUM_ITER;
        size_unary[rowidx][numLabels-1] = (double)paramIter*AUXILLARY_UNARY_VAL_MAX/NUM_ITER;
        saliency_unary[rowidx][numLabels-1] = (double)paramIter*AUXILLARY_UNARY_VAL_MAX/NUM_ITER;
    }
    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            for(int idx=(numLabels-1)*numLabels; idx<numLabels*numLabels; ++idx){
                phog_binary[rowidx][colidx][idx]=AUXILLARY_BINARY_VAL;
                loc_binary[rowidx][colidx][idx]=AUXILLARY_BINARY_VAL;
            }
            for(int idx=numLabels-1; idx<numLabels*numLabels; idx+=numLabels){
                phog_binary[rowidx][colidx][idx]=AUXILLARY_BINARY_VAL;
                loc_binary[rowidx][colidx][idx]=AUXILLARY_BINARY_VAL;
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
                       vector<TypeGeneral::REAL*>& saliency_unary,vector<TypeGeneral::REAL*>& size_unary,vector<TypeGeneral::REAL*>& objness_unary,
                       vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& loc_binary,
                       vector<vector<TypeGeneral::REAL*> >& phog_binary, const int& nodeNum){

    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        delete [] unary[rowidx];
        delete [] handpath_unary[rowidx];
        delete [] skin_unary[rowidx];
        delete [] saliency_unary[rowidx];
        delete [] size_unary[rowidx];
#if INCLUDE_OBJNESS_POTENTIAL
        delete [] objness_unary[rowidx];
#endif
    }

    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
            delete [] binary[rowidx][colidx];
            delete [] loc_binary[rowidx][colidx];
            delete [] phog_binary[rowidx][colidx];
        }
    }
}

double calculate_solution_gt_overlap(const vector<vector<vector<Rect_<int> > > >& refTubes, const vector<vector<vector<Rect_<int> > > >& tstTubes,
                                     const vector<vector<pair<int,int> > >& tubeBounds, const vector<int>& indices){
    double goodness=0, tempGoodness;
    int num_videos=0;
    int startIdx, endIdx;
    int   tstIdx;
    float areaI, areaU;
    Rect  refRect, tstRect;

    for(unsigned int shotIdx=0; shotIdx<refTubes.size(); ++shotIdx){
        tstIdx = indices[shotIdx];
        if(indices[shotIdx]>=0){
            num_videos +=1;
            tempGoodness = 0;

            startIdx=0;
            endIdx=tstTubes[shotIdx][tstIdx].size();
#if MPII_COOKING_DATASET
            startIdx = tubeBounds[shotIdx][tstIdx].first;
            endIdx = tubeBounds[shotIdx][tstIdx].second;
#endif

            for(unsigned int idx=startIdx; idx<endIdx; ++idx){

                refRect = refTubes[shotIdx][0][idx];
                tstRect = tstTubes[shotIdx][tstIdx][idx];

                Rect inter = refRect & tstRect;
                areaI = inter.width*inter.height;
                areaU = refRect.width*refRect.height + tstRect.width*tstRect.height + FLT_MIN;

//                cout << areaI << " " << areaU << "\t" << flush;
                tempGoodness += (double)areaI/areaU;
            }
            tempGoodness /= tstTubes[shotIdx][tstIdx].size();
            goodness += tempGoodness;
        }
    }

    if(num_videos>0) goodness/=num_videos;
    return goodness;
}

double calculate_precision_recall(const vector<vector<vector<Rect_<int> > > >& refTubes, const vector<vector<vector<Rect_<int> > > >& tstTubes, const vector<int>& indices,
                                     double& precision, double& recall){

    int tp=0, fp=0, fn=0;
    Rect refBb, tstBb, inter;
    double iouTh = 0.5, iou;

    precision = 0;
    recall = 0;

    for(unsigned int vidx=0; vidx<refTubes.size(); ++vidx){

        // if no tube is selected
        if(indices[vidx]<0){
            fn += tstTubes[vidx][0].size();
        }
        // if a tube is selected
        else{
            for(unsigned int fidx=0; fidx<tstTubes[vidx][indices[vidx]].size(); ++fidx){
                refBb = refTubes[vidx][0][fidx];
                tstBb = tstTubes[vidx][indices[vidx]][fidx];
                inter = refBb & tstBb;
                iou = (double)(inter.width*inter.height)/(refBb.width*refBb.height+tstBb.width*tstBb.height-inter.width*inter.height);
                if(iou>=iouTh){
                    tp += 1;
                }
                else{
                    fp += 1;
                    fn += 1;
                }
            }
        }
        // end of processing this video
    }

    if(tp + fp > 0)precision = tp/((double)tp + fp);
    if(tp + fn > 0) recall = tp/((double)tp + fn);

}

/** Videoshots related **/
// read the anno file
int ObjectAnno::loadVideoAnno(StructParam* par){

    Config configFile;
    try {
        configFile.readFile(par->action_list_file.c_str());
    }
    catch(const FileIOException &fioex) {
        cerr << "Could not read action list file " << par->action_list_file << endl;
        return -1;
    }
    catch(const ParseException &pex) {
        cerr << "Parse error: " << par->action_list_file << endl;
        return -1;
    }

#if ETHZ_ACTION_DATASET
    // set paths for each action sequence
    try{

        const Setting &action_files = configFile.getRoot()["action_files"];
        vVideoAnno.resize(action_files.getLength());
        for(unsigned int i=0;i<vVideoAnno.size();++i) {

            string object_name = (const char*)action_files[i]["object_name"];
            string seq_name    = (const char*)action_files[i]["seq_name"];

            vVideoAnno[i].object_name          = object_name;
            vVideoAnno[i].seq_name             = seq_name;
            vVideoAnno[i].param                = par;
            vVideoAnno[i].start_frame          = action_files[i]["start_frame"];
            vVideoAnno[i].end_frame            = action_files[i]["end_frame"];
            vVideoAnno[i].bmf_file_path        = par->data_path+"table-object"+object_name+"/data/s0"+seq_name+".bmf";
            vVideoAnno[i].rgb_file_path        = par->data_path+"processed_data/table-object"+object_name+"/for_opticalflow";
            vVideoAnno[i].dep_file_path        = par->data_path+"processed_data/table-object"+object_name+"/for_superpixel";
            vVideoAnno[i].opt_file_path        = par->data_path+"processed_data/table-object"+object_name+"/opticalflow";
            vVideoAnno[i].spx_file_path        = par->data_path+"processed_data/table-object"+object_name+"/superpixel";
            vVideoAnno[i].jnt_rgb_un_file_path = par->data_path+"processed_data/table-object"+object_name+"/joints_rgb/"+par->joints_rgb_idn+"_0"+seq_name+".txt";
            vVideoAnno[i].jnt_dep_file_path    = par->data_path+"processed_data/table-object"+object_name+"/joints_rgb/"+par->joints_dep_idn+"_0"+seq_name+".txt";
            vVideoAnno[i].pos_3d_file_path     = par->data_path+"/table-object"+object_name;

//            cout<<vVideoAnno[i].bmf_file_path<<endl;
//            cout<<vVideoAnno[i].spx_file_path<<endl;
//            cout<<vVideoAnno[i].jnt_rgb_un_file_path<<endl;
        }
    }
#elif CAD_120_DATASET
    // set paths for each action sequence
    try{

        const Setting &action_files = configFile.getRoot()["action_files"];
        vVideoAnno.resize(action_files.getLength());
        for(unsigned int i=0;i<vVideoAnno.size();++i) {

            string object_name = (const char*)action_files[i]["object_name"];
            string seq_name    = (const char*)action_files[i]["seq_name"];
            string subseq_name = (const char*)action_files[i]["subseq_name"];

            vVideoAnno[i].object_name          = object_name;
            vVideoAnno[i].seq_name             = seq_name;
            vVideoAnno[i].subseq_name          = subseq_name;
            vVideoAnno[i].param                = par;
            vVideoAnno[i].start_frame          = action_files[i]["start_frame"];
            vVideoAnno[i].start_frame          -= 1;
            vVideoAnno[i].end_frame            = action_files[i]["end_frame"];
            vVideoAnno[i].end_frame            -= 2;

            vVideoAnno[i].bmf_file_path        = par->bmf_file_path+"/"+subseq_name+".txt";
            vVideoAnno[i].rgb_file_path        = par->data_path+object_name+"/"+seq_name+"/"+subseq_name;
            vVideoAnno[i].dep_file_path        = par->data_path+object_name+"/"+seq_name+"/"+subseq_name;
            vVideoAnno[i].opt_file_path        = par->data_path+"opticalflow/"+object_name+"/"+seq_name+"/"+subseq_name;;
            vVideoAnno[i].spx_file_path        = par->data_path+"superpixels/"+object_name+"/"+seq_name+"/"+subseq_name;;
            vVideoAnno[i].jnt_rgb_un_file_path = par->jnt_path+"/"+subseq_name+".txt";

//            vVideoAnno[i].pos_3d_file_path     = par->data_path+"/table-object"+object_name;

//            cout<<vVideoAnno[i].bmf_file_path<<endl;
//            cout<<vVideoAnno[i].spx_file_path<<endl;
//            cout<<vVideoAnno[i].jnt_rgb_un_file_path<<endl;

            if(object_name.find("Subject1") != string::npos)
                vVideoAnno[i].actor_idx=0;
            else if(object_name.find("Subject3") != string::npos)
                vVideoAnno[i].actor_idx=1;
            else if(object_name.find("Subject4") != string::npos)
                vVideoAnno[i].actor_idx=2;
            else
                vVideoAnno[i].actor_idx=3;
        }
    }

#elif MPII_COOKING_DATASET
    // set paths for each action sequence
    try{

        const Setting &action_files = configFile.getRoot()["action_files"];
        vVideoAnno.resize(action_files.getLength());
        for(unsigned int i=0;i<vVideoAnno.size();++i) {

            string object_name = (const char*)action_files[i]["object_name"];
            string seq_name    = (const char*)action_files[i]["seq_name"];

            vVideoAnno[i].object_name          = object_name;
            vVideoAnno[i].seq_name             = seq_name;
            vVideoAnno[i].param                = par;
            vVideoAnno[i].start_frame          = action_files[i]["start_frame"];
            vVideoAnno[i].end_frame            = action_files[i]["end_frame"];
            vVideoAnno[i].end_frame           -= 1;

            if(vVideoAnno[i].end_frame-vVideoAnno[i].start_frame > 600){
                vVideoAnno[i].end_frame = vVideoAnno[i].start_frame+599;
            }
//            vVideoAnno[i].end_frame           -= 2;

            vVideoAnno[i].bmf_file_path        = par->data_path+"processed_data/bmf_files/"+seq_name+".bmf";
            vVideoAnno[i].rgb_file_path        = par->data_path+"images/"+seq_name+"/cam-002";
            vVideoAnno[i].dep_file_path        = par->data_path+"processed_data/table-object"+object_name+"/for_superpixel";
            vVideoAnno[i].opt_file_path        = par->data_path+"opticalFlow/"+seq_name+"/cam-002";
            vVideoAnno[i].spx_file_path        = par->data_path+"processed_data/superpixels/"+seq_name+"/cam-002";
            vVideoAnno[i].jnt_rgb_un_file_path = par->data_path+"processed_data/jnt_loc/"+seq_name+".txt";
            vVideoAnno[i].jnt_dep_file_path    = par->data_path+"processed_data/table-object"+object_name+"/joints_rgb/"+par->joints_dep_idn+"_0"+seq_name+".txt";
            vVideoAnno[i].pos_3d_file_path     = par->data_path+"/table-object"+object_name;

//            cout<<vVideoAnno[i].start_frame<<endl;
//            cout<<vVideoAnno[i].spx_file_path<<endl;
//            cout<<vVideoAnno[i].jnt_rgb_un_file_path<<endl;
        }
    }
#endif

    catch(const SettingNotFoundException &nfex) {
        cerr << "Not found in configuration file!" << endl;
        return -1;
    }
    return 0;
}

int ObjectAnno::size(){
    return vVideoAnno.size();
}


int loadTubes(const StructParam* param, const VideoAnno& anno, vector<vector<Rect> >& vvRect){

    string fileName;
    stringstream sSFr, sEFr;
    int num_tubes, num_frames;

    // get tube file name
    sSFr << anno.start_frame;
    sEFr << anno.end_frame;

#if ETHZ_ACTION_DATASET
    fileName =  param->tube_file_path+"_"+anno.object_name;
    fileName += "_"+anno.seq_name+"_"+sSFr.str()+"_"+sEFr.str()+".txt";
#elif CAD_120_DATASET
    fileName =  param->tube_file_path+"/"+anno.subseq_name+".txt";
#elif MPII_COOKING_DATASET
    fileName =  param->tube_file_path+"_"+anno.seq_name+"_"+sSFr.str()+"_"+sEFr.str()+".txt";
#endif

    vvRect.clear();
    ifstream iStream(fileName.c_str());
    iStream >> num_tubes;
    iStream >> num_frames;

    cout<<"tubes file name: "<<fileName << " " << num_tubes << " " << num_frames <<endl;

    vvRect.resize(num_tubes, vector<Rect>(num_frames));
    for(int tubeidx=0; tubeidx<num_tubes; ++tubeidx){
        for(int frameidx=0; frameidx<num_frames; ++frameidx){
            iStream >> vvRect[tubeidx][frameidx].x;
            iStream >> vvRect[tubeidx][frameidx].y;
            iStream >> vvRect[tubeidx][frameidx].width;
            iStream >> vvRect[tubeidx][frameidx].height;
        }
    }
    iStream.close();

    return 0;
}


int loadGtTubes(const StructParam* param, const VideoAnno& anno, vector<vector<Rect> >& vvRect){
    vvRect.resize(1);
    vvRect[0].resize(anno.end_frame-anno.start_frame+1);

    string fileName;
    stringstream sSFr, sEFr;
    ifstream iStream;

    // get ground truth tube file name
    sSFr << anno.start_frame;
    sEFr << anno.end_frame;

#if ETHZ_ACTION_DATASET
    fileName =  param->data_path+"processed_data/table-object"+anno.object_name+"/groundtruth/gt_0";
    fileName += anno.seq_name+"_"+sSFr.str()+"_"+sEFr.str()+".txt";
    iStream.open(fileName.c_str());
#elif CAD_120_DATASET
    fileName =  param->gt_path+"/"+anno.subseq_name+".txt";
    iStream.open(fileName.c_str());
#elif MPII_COOKING_DATASET
    fileName = param->gt_path+anno.seq_name+"/cam-002/gt_"+sSFr.str()+"_"+sEFr.str()+".txt";
    iStream.open(fileName.c_str());
    if(!iStream){
        sSFr.str(""); sEFr.str("");
        sSFr << anno.start_frame;
        sEFr << anno.end_frame-1;
        fileName =  param->gt_path+anno.seq_name+"/cam-002/gt_"+sSFr.str()+"_"+sEFr.str()+".txt";
        iStream.open(fileName.c_str());
    }
    int dummy;
    iStream >> dummy;
    iStream >> dummy;
    vvRect[0].resize(dummy);
#endif

    cout<<"gtrut file name: "<<fileName<< " size: " << vvRect[0].size() <<endl;

    for(unsigned int frameidx=0; frameidx<vvRect[0].size(); ++frameidx){
        iStream >> vvRect[0][frameidx].x;
        iStream >> vvRect[0][frameidx].y;
        iStream >> vvRect[0][frameidx].width;
        iStream >> vvRect[0][frameidx].height;
    }
    iStream.close();
    return 0;
}


#if MPII_COOKING_DATASET
int loadTubeBounds(const StructParam* param, const VideoAnno& anno, vector<pair<int,int> >& vBounds){
    string fileName;
    int num_tubes, dummy, startIdx, endIdx;
    stringstream sSFr, sEFr;
    int vFileListSize = anno.end_frame-anno.start_frame+1;

    // get tube file name
    sSFr << anno.start_frame;
    sEFr << anno.end_frame;
    fileName =  param->tube_file_path+"eeds_"+anno.seq_name+"_"+sSFr.str()+"_"+sEFr.str()+".txt";
    sSFr.str(""); sEFr.str("");

    vBounds.clear();
    ifstream iStream;

    iStream.open(fileName.c_str());
    if(!iStream){
        sSFr << anno.start_frame;
        sEFr << anno.end_frame-1;
        fileName =  param->tube_file_path+"eeds_"+anno.seq_name+"_"+sSFr.str()+"_"+sEFr.str()+".txt";
        sSFr.str(""); sEFr.str("");
        iStream.open(fileName.c_str());
    }
    cout<<fileName<<endl;

    if(iStream){
        iStream >> num_tubes;

        vBounds.resize(num_tubes);
        for(int tubeidx=0; tubeidx<num_tubes; ++tubeidx){
            iStream >> dummy;

            startIdx = dummy - param->tube_len/2;
            endIdx   = startIdx + param->tube_len;

            iStream >> dummy;
            iStream >> dummy;
            iStream >> dummy;
            iStream >> dummy;


            if(endIdx >= vFileListSize){
                endIdx   = vFileListSize-1;
                startIdx = endIdx - param->tube_len;
                if(startIdx < 0){
                    cout<<"warning: out of bounds during loading tubes; startIdx = "<<startIdx<<endl;
                    startIdx=0;
                }
            }
            else if(startIdx<0){
                startIdx = 0;
                endIdx   = startIdx + param->tube_len;
                if(endIdx >= vFileListSize){
                    cout<<"warning: out of bounds during loading tubes; endIdx = "<<endIdx<<endl;
                    endIdx = vFileListSize -1;
                }
            }

            vBounds[tubeidx].first = startIdx;
            vBounds[tubeidx].second = endIdx;
        }
        iStream.close();
    }
    else{
        cout<<"cannot open tubeseeds file: "<<fileName<<endl;
        exit(-1);
    }

    return 0;
}
#endif
