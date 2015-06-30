#include "customFunctions.h"

// Example: minimizing an energy function with Potts terms.
// See type*.h files for other types of terms.

void testGeneral(char* configFileName, int numBeliefIter, int paramIter)
{
	MRFEnergy<TypeGeneral>* mrf;
	MRFEnergy<TypeGeneral>::NodeId* nodes;
	MRFEnergy<TypeGeneral>::Options options;
	TypeGeneral::REAL energy, lowerBound;

	int numNodes; // number of nodes
    int numLabels;
	vector<TypeGeneral::REAL*>          unary_pot;
	vector<TypeGeneral::REAL*>          handpath_unary;
	vector<TypeGeneral::REAL*>          skin_unary;
	vector<TypeGeneral::REAL*>          saliency_unary;
	vector<TypeGeneral::REAL*>          size_unary;
	vector<vector<TypeGeneral::REAL*> > binary_pot;
	vector<vector<TypeGeneral::REAL*> > loc_binary;
	vector<vector<TypeGeneral::REAL*> > phog_binary;

	vector<vector<vector<Rect_<int> > > > gtTubes, detTubes;
    vector<vector<pair<int,int> > > tubeBounds;
	string writePath;

    loadPotentials(configFileName, paramIter, unary_pot, handpath_unary, skin_unary, saliency_unary, size_unary,
                   binary_pot, loc_binary, phog_binary, numNodes, numLabels, writePath, gtTubes, detTubes, tubeBounds);

//    float sal_un=0.5, skin_un=0.5, size_un=0.5, phog_bn=0.5;

    float wgt[] = {0.0,0.05,0.25,0.50,0.75,1.00};
    double overlap = 0;
    double precision, recall;

    // grid search for best parameters
    for(int sal_un=0; sal_un<6; sal_un++){
        for(int skin_un=0; skin_un<6; skin_un++){
            for(int size_un=0; size_un<6; size_un++){
                for(int handpath_un=0; handpath_un<6; handpath_un++){
                    for(int loc_bn=0; loc_bn<6; loc_bn++){
                        for(int phog_bn=0; phog_bn<6; phog_bn++){
                            // init the potentials
                            for(int nodeIdx=0; nodeIdx<numNodes; ++nodeIdx){
                                for(int labIdx=0; labIdx<numLabels; ++labIdx){
                                    unary_pot[nodeIdx][labIdx] = wgt[sal_un]*saliency_unary[nodeIdx][labIdx] + wgt[skin_un]*skin_unary[nodeIdx][labIdx] + wgt[size_un]*size_unary[nodeIdx][labIdx] + wgt[handpath_un]*handpath_unary[nodeIdx][labIdx];
                                }
                            }

                            for(int refIdx=0; refIdx<numNodes; ++refIdx){
                                for(int tstIdx=0; tstIdx<numNodes; ++tstIdx){
                                    for(int labIdx=0; labIdx<numLabels*numLabels; ++labIdx){
                                        binary_pot[refIdx][tstIdx][labIdx] = wgt[loc_bn]*loc_binary[refIdx][tstIdx][labIdx] + wgt[phog_bn]*phog_binary[refIdx][tstIdx][labIdx];
                                    }
                                }
                            }

                            // now initiate and run the mrf
                            mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
                            nodes = new MRFEnergy<TypeGeneral>::NodeId[numNodes];

                            // construct energy
                            for(int refidx=0; refidx<numNodes; ++refidx){
                                nodes[refidx] = mrf->AddNode(TypeGeneral::LocalSize(numLabels), TypeGeneral::NodeData(unary_pot[refidx]));
                            }
                            for(int refidx=0; refidx<numNodes; ++refidx){
                                for(int tstidx=0; tstidx<numNodes; ++tstidx){
                                    if(refidx<tstidx){
                                        mrf->AddEdge(nodes[refidx], nodes[tstidx], TypeGeneral::EdgeData(TypeGeneral::GENERAL, binary_pot[refidx][tstidx]));
                                    }
                                }
                            }

                            /////////////////////// TRW-S algorithm //////////////////////
                            options.m_iterMax = numBeliefIter; // maximum number of iterations
                            mrf->Minimize_TRW_S(options, lowerBound, energy);

                            //////////////// rerun TRW-S algorithm N times ///////////////
                            for(int idx=0; idx<NUM_RERUN_TRWS; ++idx){

                                // update unaries being kicked out
                                double max_val_unary=1e-10;
                                for(int nodeIdx=0; nodeIdx<numNodes; ++nodeIdx){
                                    int solIdx = mrf->GetSolution(nodes[nodeIdx]);

                                    // add binaries being kicked out to unaries
                                    for(int tstNodeIdx=0; tstNodeIdx<numNodes; ++tstNodeIdx){
//                                        if(nodeIdx == tstNodeIdx) continue;
                                        for(int tstSolIdx=0; tstSolIdx<numLabels; ++tstSolIdx){
                                            unary_pot[tstNodeIdx][tstSolIdx] += binary_pot[nodeIdx][tstNodeIdx][solIdx*numLabels+tstSolIdx];
                                            if(unary_pot[tstNodeIdx][tstSolIdx] > max_val_unary){
                                                max_val_unary = unary_pot[tstNodeIdx][tstSolIdx];
                                            }
                                        }
                                    }
                                }
                                assert(max_val_unary>0);

                                // kick out present solution
                                for(int nodeIdx=0; nodeIdx<numNodes; ++nodeIdx){
                                    int solIdx = mrf->GetSolution(nodes[nodeIdx]);
                                    if(solIdx<numLabels-1){
                                        unary_pot[nodeIdx][solIdx]=10*max_val_unary;
                                    }
                                    else{
                                        unary_pot[nodeIdx][solIdx];
                                    }
                                }

                                // reconstruct mrf
                                delete mrf;
                                delete nodes;
                                mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
                                nodes = new MRFEnergy<TypeGeneral>::NodeId[numNodes];
                                for(int refidx=0; refidx<numNodes; ++refidx){
                                    nodes[refidx] = mrf->AddNode(TypeGeneral::LocalSize(numLabels), TypeGeneral::NodeData(unary_pot[refidx]));
                                }
                                for(int refidx=0; refidx<numNodes; ++refidx){
                                    for(int tstidx=0; tstidx<numNodes; ++tstidx){
                                        if(refidx<tstidx){
                                            mrf->AddEdge(nodes[refidx], nodes[tstidx], TypeGeneral::EdgeData(TypeGeneral::GENERAL, binary_pot[refidx][tstidx]));
                                        }
                                    }
                                }

                                // resolve mrf
                                options.m_iterMax = numBeliefIter;
                                mrf->Minimize_TRW_S(options, lowerBound, energy);
                            }

                            ////////////////// solution-gt overlap //////////////////////
                            vector<int> solution_tube_ids(numNodes,0);
                            for(int idx=0; idx<numNodes; ++idx){
                                solution_tube_ids[idx] = mrf->GetSolution(nodes[idx]);
                                if(solution_tube_ids[idx] == numLabels-1)
                                    solution_tube_ids[idx] = -1;
                            }
                            overlap = calculate_solution_gt_overlap(gtTubes,detTubes,tubeBounds,solution_tube_ids);
                            calculate_precision_recall(gtTubes, detTubes, solution_tube_ids, precision, recall);


            #if WRITE_FLAG
                            // save the solution configuration
                            ofstream oStream(writePath.c_str(),ios::app);
                            oStream<<wgt[sal_un]<<" "<<wgt[skin_un]<<" "<<wgt[size_un]<<" "<<wgt[handpath_un]<<" "<<wgt[loc_bn]<<" "<<wgt[phog_bn]<<"\t";
                            for(int idx=0; idx<numNodes; ++idx){
                                oStream<<mrf->GetSolution(nodes[idx])<<" ";
                                cout<<mrf->GetSolution(nodes[idx])<<" ";
                            }
                            cout<< overlap <<" "<< precision <<" "<< recall <<endl;
                            oStream<<overlap <<" "<< precision <<" "<< recall <<endl;
                            oStream.close();
            #else
                            for(int idx=0; idx<numNodes; ++idx){
                                cout<<mrf->GetSolution(nodes[idx])+1<<" ";
                            }
                            cout<<overlap <<" "<< precision <<" "<< recall <<endl;
            #endif

                            // done
                            delete nodes;
                            delete mrf;
                        }
                    }
                }
            }
        }
    }

	// finally done
	destroyPotentials(unary_pot, handpath_unary, skin_unary, saliency_unary, size_unary, binary_pot, loc_binary, phog_binary, numNodes);
}

int main(int argc, char* argv[])
{
    if(argc!=3){
        cout<<"usage: <configFilePath> <numIter>"<<endl;
        return -1;
    }

#if REPLACE_LAST_TUBE_AS_NULL
    for(int iterIdx=0; iterIdx<=NUM_ITER; ++iterIdx){
	    testGeneral(argv[1],atoi(argv[2]),iterIdx);
	}
#else
    testGeneral(argv[1],atoi(argv[2]),NUM_ITER+1);
#endif
	return 0;
}


///** TESTING PREST MODEL **/
////#include <stdio.h>
////#include "MRFEnergy.h"
////#include <vector>
////#include <string>
////#include <iostream>
////#include <fstream>
////#include <stdlib.h>
////#include <math.h>
////
////using namespace std;
////#define WRITE_FLAG 0
////
////// Example: minimizing an energy function with Potts terms.
////// See type*.h files for other types of terms.
////
////void loadPotentials(char* configFile, vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& saliency_unary,
////                    vector<TypeGeneral::REAL*>& phog_unary, vector<TypeGeneral::REAL*>& surf_unary,
////                    vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& phog_binary, vector<vector<TypeGeneral::REAL*> >& surf_binary,
////                    int& nodeNum, int& numLabels, string& writePath){
////
////    string saliency_unary_path, phog_unary_path, surf_unary_path, phog_binary_path, surf_binary_path, filtered_tube_ind_path;
////    TypeGeneral::REAL temp;
////    int row, col;
////    float maxWeight=1000;
////
////    // read the potential file names
////    ifstream iStream(configFile);
////    iStream >> saliency_unary_path;
////    iStream >> phog_unary_path;
////    iStream >> surf_unary_path;
////    iStream >> phog_binary_path;
////    iStream >> surf_binary_path;
////    iStream >> writePath;
////    iStream.close();
////
////    // read saliency_unary potentials
////    iStream.open(saliency_unary_path.c_str());
////    iStream >> nodeNum;
////    iStream >> numLabels;
////    saliency_unary.resize(nodeNum);
////    unary.resize(nodeNum);
////    for(unsigned int rowidx=0; rowidx<saliency_unary.size(); ++rowidx){
////        saliency_unary[rowidx] = new TypeGeneral::REAL [numLabels];
////        unary[rowidx] = new TypeGeneral::REAL [numLabels];
////        for(int idx=0; idx<numLabels; ++idx){
////            iStream >> saliency_unary[rowidx][idx];
////            saliency_unary[rowidx][idx] = 1-saliency_unary[rowidx][idx];
////        }
////    }
////    iStream.close();
////
////    // read phog_unary potentials
////    iStream.open(phog_unary_path.c_str());
////    iStream >> nodeNum;
////    iStream >> numLabels;
////    phog_unary.resize(nodeNum);
////    for(unsigned int rowidx=0; rowidx<phog_unary.size(); ++rowidx){
////        phog_unary[rowidx] = new TypeGeneral::REAL [numLabels];
////        for(int idx=0; idx<numLabels; ++idx){
////            iStream >> phog_unary[rowidx][idx];
////            phog_unary[rowidx][idx] = phog_unary[rowidx][idx];
////        }
////    }
////    iStream.close();
////
////    // read surf_unary potentials
////    iStream.open(surf_unary_path.c_str());
////    iStream >> nodeNum;
////    iStream >> numLabels;
////    surf_unary.resize(nodeNum);
////    for(unsigned int rowidx=0; rowidx<surf_unary.size(); ++rowidx){
////        surf_unary[rowidx] = new TypeGeneral::REAL [numLabels];
////        for(int idx=0; idx<numLabels; ++idx){
////            iStream >> surf_unary[rowidx][idx];
////            surf_unary[rowidx][idx] = surf_unary[rowidx][idx];
////        }
////    }
////    iStream.close();
////
////    // read phog_binary potentials
////    iStream.open(phog_binary_path.c_str());
////    iStream >> nodeNum;
////    iStream >> numLabels;
////    phog_binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
////    binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
////    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
////        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
////            iStream >> row;
////            iStream >> col;
////            phog_binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
////            binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
////            for(int idx=0; idx<numLabels*numLabels; ++idx){
////                iStream >> phog_binary[row][col][idx];
////            }
////        }
////    }
////    iStream.close();
////
////    // read phog_binary potentials
////    iStream.open(surf_binary_path.c_str());
////    iStream >> nodeNum;
////    iStream >> numLabels;
////    surf_binary.resize(nodeNum,vector<TypeGeneral::REAL*>(nodeNum));
////    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
////        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
////            iStream >> row;
////            iStream >> col;
////            surf_binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
////            binary[row][col] = new TypeGeneral::REAL [numLabels*numLabels];
////            for(int idx=0; idx<numLabels*numLabels; ++idx){
////                iStream >> surf_binary[row][col][idx];
////            }
////        }
////    }
////    iStream.close();
////
////#if WRITE_FLAG
////    // clear the write file file
////    ofstream oStream(writePath.c_str());
////    oStream.close();
////#endif
////
////}
////
////void destroyPotentials(vector<TypeGeneral::REAL*>& unary, vector<TypeGeneral::REAL*>& saliency_unary, vector<TypeGeneral::REAL*>& phog_unary,
////                       vector<TypeGeneral::REAL*>& surf_unary,
////                       vector<vector<TypeGeneral::REAL*> >& binary, vector<vector<TypeGeneral::REAL*> >& phog_binary,
////                       vector<vector<TypeGeneral::REAL*> >& surf_binary, const int& nodeNum){
////
////    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
////        delete [] unary[rowidx];
////        delete [] saliency_unary[rowidx];
////        delete [] phog_unary[rowidx];
////        delete [] surf_unary[rowidx];
////    }
////
////    for(unsigned int rowidx=0; rowidx<nodeNum; ++rowidx){
////        for(unsigned int colidx=0; colidx<nodeNum; ++colidx){
////            delete [] binary[rowidx][colidx];
////            delete [] phog_binary[rowidx][colidx];
////            delete [] surf_binary[rowidx][colidx];
////        }
////    }
////}
////
////void testGeneral(char* configFileName, int numIter)
////{
////	MRFEnergy<TypeGeneral>* mrf;
////	MRFEnergy<TypeGeneral>::NodeId* nodes;
////	MRFEnergy<TypeGeneral>::Options options;
////	TypeGeneral::REAL energy, lowerBound;
////
////	int numNodes; // number of nodes
////    int numLabels;
////	vector<TypeGeneral::REAL*>          unary_pot;
////	vector<TypeGeneral::REAL*>          saliency_unary;
////	vector<TypeGeneral::REAL*>          phog_unary;
////	vector<TypeGeneral::REAL*>          surf_unary;
////	vector<vector<TypeGeneral::REAL*> > binary_pot;
////	vector<vector<TypeGeneral::REAL*> > surf_binary;
////	vector<vector<TypeGeneral::REAL*> > phog_binary;
////	string writePath;
////
////    loadPotentials(configFileName, unary_pot, saliency_unary, phog_unary, surf_unary,
////                   binary_pot, phog_binary, surf_binary, numNodes, numLabels, writePath);
////
//////    float sal_un=0.5, skin_un=0.5, size_un=0.5, phog_bn=0.5;
////
////    float wgt[] = {0.0,0.05,0.25,0.50,0.75,1.00};
////
////    // grid search for best parameters
////    for(int sal_un=0; sal_un<6; sal_un++){
////        for(int phog_un=0; phog_un<6; phog_un++){
////            for(int surf_un=0; surf_un<6; surf_un++){
////
////                for(int phog_bn=0; phog_bn<6; phog_bn++){
////                    for(int surf_bn=0; surf_bn<6; surf_bn++){
////                        // init the potentials
////                        for(int nodeIdx=0; nodeIdx<numNodes; ++nodeIdx){
////                            for(int labIdx=0; labIdx<numLabels; ++labIdx){
////                                unary_pot[nodeIdx][labIdx] =  wgt[sal_un]*saliency_unary[nodeIdx][labIdx];
////                                unary_pot[nodeIdx][labIdx] += wgt[phog_un]*phog_unary[nodeIdx][labIdx];
////                                unary_pot[nodeIdx][labIdx] += wgt[surf_un]*surf_unary[nodeIdx][labIdx];
////                            }
////                        }
////
////                        for(int refIdx=0; refIdx<numNodes; ++refIdx){
////                            for(int tstIdx=0; tstIdx<numNodes; ++tstIdx){
////                                for(int labIdx=0; labIdx<numLabels*numLabels; ++labIdx){
////                                    binary_pot[refIdx][tstIdx][labIdx] = wgt[surf_bn]*surf_binary[refIdx][tstIdx][labIdx] + wgt[phog_bn]*phog_binary[refIdx][tstIdx][labIdx];
////                                }
////                            }
////                        }
////
////                        // now initiate and run the mrf
////                        mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
////                        nodes = new MRFEnergy<TypeGeneral>::NodeId[numNodes];
////
////                        // construct energy
////                        for(int refidx=0; refidx<numNodes; ++refidx){
////                            nodes[refidx] = mrf->AddNode(TypeGeneral::LocalSize(numLabels), TypeGeneral::NodeData(unary_pot[refidx]));
////                        }
////                        for(int refidx=0; refidx<numNodes; ++refidx){
////                            for(int tstidx=0; tstidx<numNodes; ++tstidx){
////                                if(refidx<tstidx){
////                                    mrf->AddEdge(nodes[refidx], nodes[tstidx], TypeGeneral::EdgeData(TypeGeneral::GENERAL, binary_pot[refidx][tstidx]));
////                                }
////                            }
////                        }
////
////                        /////////////////////// TRW-S algorithm //////////////////////
////                        options.m_iterMax = numIter; // maximum number of iterations
////                        mrf->Minimize_TRW_S(options, lowerBound, energy);
////
////        #if WRITE_FLAG
////                        // save the solution configuration
////                        ofstream oStream(writePath.c_str(),ios::app);
////                        oStream<<wgt[sal_un]<<" "<<wgt[phog_un]<<" "<<wgt[surf_un]<<" "<<wgt[phog_bn]<<" "<<wgt[surf_bn]<<"\t";
////                        for(int idx=0; idx<numNodes; ++idx){
////                            oStream<<mrf->GetSolution(nodes[idx])<<" ";
////                            cout<<mrf->GetSolution(nodes[idx])<<" ";
////                        }
////                        cout<<endl;
////                        oStream<<endl;
////                        oStream.close();
////        #else
////                        for(int idx=0; idx<numNodes; ++idx){
////                            cout<<mrf->GetSolution(nodes[idx])+1<<" ";
////                        }
////                        cout<<endl;
////        #endif
////
////                        // done
////                        delete nodes;
////                        delete mrf;
////                    }
////                }
////            }
////        }
////    }
////
////	// finally done
////	destroyPotentials(unary_pot, saliency_unary, phog_unary, surf_unary, binary_pot, phog_binary, surf_binary, numNodes);
////}
////
////int main(int argc, char* argv[])
////{
////    if(argc!=3){
////        cout<<"usage: <configFilePath> <numIter>"<<endl;
////        return -1;
////    }
////
////	testGeneral(argv[1],atoi(argv[2]));
////	return 0;
////}
