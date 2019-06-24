#ifndef __PP_FILTER_EP_H__
#define __PP_FILTER_EP_H__

#include <SubsysReco.h>
#include "RpConst.h"
#include <vector>
#include <string>

class TFile;
class TTree;
class TF1;
class TH1;

class Fun4AllHistoManager;
class TFvtxCompactTrkMap;
class PHCentralTrack;

class BbcRaw;
class BbcCalib;
class BbcGeo;
class DoubleInteractionUtil;

class pp_filter_ep: public SubsysReco 
{
  public:
  pp_filter_ep(const char* output = "filtered_tree.root");
  virtual ~pp_filter_ep();
    
  virtual int Init(PHCompositeNode *topNode);
  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode);

  //int GetFVTXTrack(TFvtxCompactTrkMap* trkfvtx_map);
  //int GetCNTPC3_TOF(PHCentralTrack* central);
  //void GetBBCPMT(BbcRaw* bbcraw);

  double calcsdphi(double dphi, int arm, int ch, double mom);
  double calcsdz(double dphi, int arm, int ch, double mom);
  
  private:
  std::string _filename;
  TF1 *fx;
  TF1 *fy;

  Fun4AllHistoManager*  HistoManager;

  BbcCalib* bbccalib;
  BbcGeo*   bbcgeo; 
  DoubleInteractionUtil *doubleint_util;

  void  initHisto();

  int runnumber;

  TH1* h_Fvtx_phi[10];
  TH1* h_Bbc_qx[6][10];
  TH1* h_Bbc_qy[6][10];
  TH1* h_Bbc_psi[6][10];
	TH1* h_Bbc_psi_cos[6][10];
	TH1* h_Bbc_psi_sin[6][10];
  TH1* h_Cnt_qx[6][10];
  TH1* h_Cnt_qy[6][10];
  TH1* h_Cnt_psi[6][10];
};

#endif
