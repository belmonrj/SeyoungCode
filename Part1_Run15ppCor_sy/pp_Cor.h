#ifndef __PP_COR_H__
#define __PP_COR_H__

#include <SubsysReco.h>
#include "RpConst.h"
#include <vector>
#include <string>

class PreviousEvent;
class Fun4AllHistoManager;

class BbcCalib;
class BbcGeo; 
class BbcRaw;

class VtxOut;

class TFvtxCompactTrkMap;

class mpcTowerContainer;
class mpcRawContainer;
class TRxnpRawScintMap;

class PHCentralTrack;

class TFile;
class TProfile;
class TProfile2D;
class TH1;
class TH2;
class TH3;

class PHGlobal;

class pp_Cor: public SubsysReco 
{

 public:
  pp_Cor(const char* output = "pp_Cor");
  virtual ~pp_Cor();
    
  virtual int Init(PHCompositeNode *topNode);
  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode);
  
  virtual int Reset(PHCompositeNode *topNode) { return 0; }
  virtual int ResetEvent(PHCompositeNode *topNode) { return 0; }
  
 protected:
	int CreateNodeTree(PHCompositeNode *topNode) { return 0; }
  PHCompositeNode *dst_node;
  
 private:
  static const int NCENTANA = 4; //Number of centrality bins
  static const int NZANA    =10; //Number of z-vertex bins
  static const int NEVE     = 4; //Number of events in buffer
  static const int NTRIG    = 3; //Number of triggers (MB, FVTX-1, FVTX-2)
  static const int NPTBIN   = 5; //Number of pt bins


  Fun4AllHistoManager*  HistoManager;
  
  //For BBC RP
  BbcCalib* bbccalib;
  BbcGeo*   bbcgeo; 
  
  void  initHisto();

	//Get track information
  bool GetCNTPC3_TOF(const int ic_ana, const int iz_ana, const int ieve, PHCentralTrack* trk,PHGlobal* global);
  bool GetBBCPMT    (const int ic_ana, const int iz_ana, const int ieve, BbcRaw* bbcraw);
  bool GetFVTXTrack(const int ic_ana,const int iz_ana, const int ieve, PHGlobal* global,   TFvtxCompactTrkMap* trkfvtx_map);
  int GetFVTXCentralityBin (TFvtxCompactTrkMap* trkfvtx_map);

	//Reset buffer
  void do_Buff_Cnt_Fvt_Trk_Reset  (const int ic_ana, const int iz_ana);

	//Fill correlation histograms
  void do_Real_Mix_CNT_BBC_Pmt_Cor(const int ic_ana, const int iz_ana);
  void do_Real_Mix_BBC_S_N_Pmt_Cor(const int ic_ana, const int iz_ana);
  void do_Real_Mix_CNT_CNT_Cor    (const int ic_ana, const int iz_ana);
  void do_Real_Mix_CNT_FVT_Trk_Cor(const int ic_ana, const int iz_ana);
  void do_Real_Mix_FVT_S_N_Trk_Cor(const int ic_ana, const int iz_ana);

  std::string OutputFileName;
  int         m_RunNumber;
  int         m_EventNumber;
  
  //2PC
  int   m_Eve             [NCENTANA][NZANA];
  
  float m_Fvt_Phi         [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_Eta         [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_DcaX        [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_DcaY        [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_DcaX_pre    [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_DcaY_pre    [NCENTANA][NZANA][NEVE][1000];
  float m_Fvt_Ntrk        [NCENTANA][NZANA][NEVE]      ;

  float m_Cnt_Phi         [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Pt          [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Ntrk        [NCENTANA][NZANA][NEVE]      ;
  float m_Cnt_Pt_Bin      [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Eta         [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Phi1        [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Phi3        [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Zpc1        [NCENTANA][NZANA][NEVE][1000];
  float m_Cnt_Zpc3        [NCENTANA][NZANA][NEVE][1000];
  
  float m_Buf_BbcPmt_Phi  [NCENTANA][NZANA][NEVE][128];
  float m_Buf_BbcPmt_Eta  [NCENTANA][NZANA][NEVE][128];
  float m_Buf_BbcPmt_Cha  [NCENTANA][NZANA][NEVE][128];

  float fvtx_mult	  [3]; // South, North, Total
  float fvtx_x, fvtx_y, fvtx_z; 
  // event
  TH2* m_centzv[3];
  TH1* hFvtxTrk[3]; //SOUTH-NORTH-BOTH
  TH2* h_Cnt_Bbc_dPhi_Real[NCENTANA][2]; //[NCENTANA][ETA-BIN]
  TH2* h_Cnt_Bbc_dPhi_Mix [NCENTANA][2];
  TH1* h_Bbc_S_N_dPhi_Real[NCENTANA];  
  TH1* h_Bbc_S_N_dPhi_Mix [NCENTANA];
  TH2* h_Cnt_Bbc_dPhi_dEta_Real[NCENTANA][NPTBIN]; //[NCENTANA][PT-BIN]
  TH2* h_Cnt_Bbc_dPhi_dEta_Mix [NCENTANA][NPTBIN];;
  TH2* h_Cnt_Bbc_dPhi_dEta_Real_sqr[NCENTANA][NPTBIN]; // [NCENTANA][PT-BIN]
  TH2* h_Cnt_Bbc_dPhi_dEta_Mix_sqr [NCENTANA][NPTBIN];
  TH2* h_BbV_Fvt_Dca	[NCENTANA][2];
  TH2* h_BbV_Fvt_z	[NCENTANA][2];

  TH2* h_Bbc_S_N_dPhi_dEta_Real    [NCENTANA];
  TH2* h_Bbc_S_N_dPhi_dEta_Mix     [NCENTANA];
  TH2* h_Bbc_S_N_dPhi_dEta_Real_sqr[NCENTANA];
  TH2* h_Bbc_S_N_dPhi_dEta_Mix_sqr [NCENTANA];
    
  TH2* h_Cnt_Cnt_dPhi_dEta_Real[NCENTANA][NPTBIN][NPTBIN]; //[NCENTANA][PT-TRIG][PT-ASSOC]
  TH2* h_Cnt_Cnt_dPhi_dEta_Mix [NCENTANA][NPTBIN][NPTBIN];

  TH1* h_Fvt_S_N_dPhi_Real     [NCENTANA][NZANA][4]; //[NCENTANA][Z-BIN][SYS?]
  TH1* h_Fvt_S_N_dPhi_Mix      [NCENTANA][NZANA][4];
  TH2* h_Fvt_S_N_dPhi_dEta_Real[NCENTANA][NZANA][4];
  TH2* h_Fvt_S_N_dPhi_dEta_Mix [NCENTANA][NZANA][4];
  TH1* h_Cnt_Fvt_dPhi_Real     [NCENTANA][NZANA][2][4]; //[NCENTANA][Z-BIN][ETA-BIN][SYS?]
  TH1* h_Cnt_Fvt_dPhi_Mix      [NCENTANA][NZANA][2][4];
  TH2* h_Cnt_Fvt_dPhi_dEta_Real[NCENTANA][NZANA][NPTBIN][4]; //[NCENTANA][Z-BIN][PT-BIN][SYS]
  TH2* h_Cnt_Fvt_dPhi_dEta_Mix [NCENTANA][NZANA][NPTBIN][4];
  TH2* h_Fvt_phi_eta[NZANA][2];
  TH1* h_Fvt_Dca  [2][2][3];//[ARM][XY][MODIFIED]
  TH2* h_Bbc_S_N_correlation_Real[NCENTANA];
  TH2* h_Bbc_S_N_correlation_Mix[NCENTANA];

  TH2* h_dPc1_vs_dZpc1[2][2];
  TH2* h_dPc3_vs_dZpc3[2][2];
  TH2* h_NCnt_vs_NDch;
  TH1* h_Cnt_pt      ;
  TH2* h_Cnt_phi_vs_eta;
  TH2* m_FVTXtrk_BBCcha[3];

  TProfile2D* h_BbcS_x_vs_y[NCENTANA];
  TProfile2D* h_BbcN_x_vs_y[NCENTANA];
  TProfile* h_BbcS_ipmt[NCENTANA];
  TProfile* h_BbcN_ipmt[NCENTANA];

};

#endif
