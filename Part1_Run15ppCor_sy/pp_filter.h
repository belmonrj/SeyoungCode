#ifndef __PP_FILTER_H__
#define __PP_FILTER_H__

#include <SubsysReco.h>
#include "RpConst.h"
#include <vector>
#include <string>

class TFile;
class TTree;
class TF1;

class TFvtxCompactTrkMap;
class PHCentralTrack;

class BbcRaw;
class BbcCalib;
class BbcGeo;
class DoubleInteractionUtil;

class pp_filter: public SubsysReco 
{

 public:
  //pp_filter();
  pp_filter(const char* output = "filtered_tree.root");
  virtual ~pp_filter();
    
  virtual int Init(PHCompositeNode *topNode);
  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode);

	//int GetFVTXTrack(TFvtxCompactTrkMap* trkfvtx_map);
	//int GetCNTPC3_TOF(PHCentralTrack* central);
	//void GetBBCPMT(BbcRaw* bbcraw);

	double calcsdphi(double dphi, int arm, int ch, double mom);
	double calcsdz(double dphi, int arm, int ch, double mom);
	double getBBCX(float zvtx);
	double getBBCY(float zvtx);
 
 private:

	std::string _filename;

	TTree *T;
	TF1 *fx;
	TF1 *fy;

  BbcCalib* bbccalib;
  BbcGeo*   bbcgeo; 
	DoubleInteractionUtil *doubleint_util;

	int runnumber;

	/*

	*/

};

#endif
