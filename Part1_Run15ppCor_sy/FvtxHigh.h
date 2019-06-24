#ifndef __FVTXHIGH_H__
#define __FVTXHIGH_H__

#include <SubsysReco.h>
#include "RpConst.h"
#include <vector>
#include <string>

class PreviousEvent;
class RpSumXYObject;
class ReactionPlaneObject;

class Fun4AllHistoManager;
class ReactionPlaneCalibv1;

class TFvtxCompactTrkMap;

class SvxTrackMap;
class SvxSegmentList;
class SvxCentralTrackList;
class SvxCentralTrac;
class SvxCentralTrack;

class VtxOut;

class BbcCalib;
class BbcGeo; 
class BbcRaw;

class mpcTowerContainer;
//class mpcRawContainer;

class PHCentralTrack;

class TFile;
class TProfile;
class TProfile2D;
class TH1;
class TH2;
class TH3;
class TNtuple;

class PHGlobal;

class FvtxHigh: public SubsysReco 
{
 public:
  enum  ANA_FLAG{FLAG_EP_METHOD = 0, FLAG_2PC_METHOD=1, FLAG_CHECK_FVTX_TRK=2, FLAG_CNT_CHECK=3, FLAG_FVTX_2PC = 4};
 public:
  FvtxHigh(const char* output = "FvtxHigh");
  virtual ~FvtxHigh();
  
	/*
  virtual int  Init(   PHCompositeNode *topNode);
  virtual int  InitRun(PHCompositeNode *topNode);
  virtual int  process_event(PHCompositeNode *topNode);
  virtual int  End(PHCompositeNode *topNode);
  
  virtual int  Reset(     PHCompositeNode *topNode) { return 0; }
  virtual int  ResetEvent(PHCompositeNode *topNode) { return 0; }
  //--    virtual void Print(const char *what) const { return; }
  
  virtual void setAnaFlag  (int flag) { m_Ana_Flag = flag; }
  
 protected:
  int    CreateNodeTree(PHCompositeNode *topNode) { return 0; }
  
  int m_Ana_Flag;    // 0 : EP method
  // 1 : 2PC method
  
 private:
  static const unsigned int NHAR   = RP::NHAR4;
  static const unsigned int NHAR2D = 2;
  static const unsigned int NDET   = 7; // SVX, SEG, MPC, BBC, SMD, CNT, FVTX
  static const unsigned int NKIND  =70; // SVX, SEG, MPC, BBC, SMD, CNT, FVTX
  static const int NCENTANA =10;
  static const int NZANA    =12;
  static const int NKINDANA =70;
  static const int NDETANA  = 7;//BBC,SMD,CNT,MPC
  static const int NETAANA  = 8;
  static const int NETA     = 5;
  static const int NEVE     = 4;// real=1, mix =3;
  
  ReactionPlaneCalibv1* rpcalibv1;
  Fun4AllHistoManager*  HistoManager;
  //TFile* file;
  
  //For BBC RP
  BbcCalib* bbccalib;
  BbcGeo*   bbcgeo; 
  
  void  initHisto();
  void  Write_Histo();
  bool  isBadTick(PreviousEvent *peve);
  
  bool GetCNTPC3_TOF(PHCentralTrack* central);
  bool GetFVTXTrack (PHGlobal* global, TFvtxCompactTrkMap* trkfvtx_map, const int iz_ana, const int itrg);
  
  bool Get_Enable_Det_Self(const int detid);
  bool Get_Enable_Det_Self_vnRP(const int detid);

  void do_Buff_Cnt_Fvt_Trk_Reset();
    
  std::string OutputFileName;
  //char        *outNam;
  int         m_RunNumber;
  int         m_EventNumber;
  float       m_bgain[128];
  int         m_eve;
  float       m_ringgainS[4];
  float       m_ringgainN[4];
  float       m_bbcgainflag;
    
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
    float m_Cnt_Pc3sdZ      [NCENTANA][NZANA][NEVE][1000];
    float m_Cnt_Pc3sdPhi    [NCENTANA][NZANA][NEVE][1000];
    float m_Cnt_Ntrk        [NCENTANA][NZANA][NEVE]      ;

    float m_Buf_BbcPmt_Phi  [NCENTANA][NZANA][NEVE][128];
    float m_Buf_BbcPmt_Eta  [NCENTANA][NZANA][NEVE][128];
    float m_Buf_BbcPmt_Cha  [NCENTANA][NZANA][NEVE][128];
    
    float m_Buf_Mpc_Et      [NCENTANA][NZANA][NEVE][400];
    float m_Buf_Mpc_Phi     [NCENTANA][NZANA][NEVE][400];
    float m_Buf_Mpc_Eta     [NCENTANA][NZANA][NEVE][400];
    float m_Buf_MpcNTower   [NCENTANA][NZANA][NEVE];
    
    float m_Buf_FVTRP[NCENTANA][NZANA][NEVE][46][RP::NHAR4];
    float m_Buf_SMDRP[NCENTANA][NZANA][NEVE][3];
    float m_Buf_CNTRP[NCENTANA][NZANA][NEVE][RP::NHAR4];
    //--    int   EventSkip();
    
    // event
    TH2* m_centzv;
    TH2* m_centzv_svx;
    
    // recentering & flattening
    // ndet= Svx,     Seg, MPC, BBC, SMD, CNT
    // kind= 12+1+12   9    3    3    3    5
    TProfile** m_rawq  [RP::NMUL3][RP::NZPS2][RP::NHAR4][NDET]; // Raw Q vector to get recentering parameters
    TH3**      m_rawqxy[RP::NMUL3][NHAR2D][NDET];    // Raw Qx vs Qy vs iz

    // flattening
    TProfile** m_flatpar[RP::NMUL3][RP::NZPS2][RP::NHAR4][NDET]; // Raw RP to get flattening parameters
    TH1**      m_rawrp[RP::NMUL3/2][RP::NZPS2][RP::NHAR4][NDET];   // Raw RP to get flattening parameters
    TH1**      m_CalibRp[RP::NMUL3/2][RP::NHAR4][NDET];   // Raw RP to get flattening parameters    
    
    TProfile* bbcgain;
    TProfile* bbcgainn_ring;
    TProfile* bbcgains_ring;
    TH2* bbcmapp;
    TH2* bbcimap[2][4];
    TH2* bbcrad[2][4];
    TProfile2D* bbcpmt2dmaps[5];
    TProfile2D* bbcpmt2dmapn[5];
    
    TProfile* RPCORSN_COS [NDETANA][2][RP::NHAR4];
    TProfile* RPCORSN_SIN [NDETANA][2][RP::NHAR4];
    TProfile* RPCORCNT_COS[NDETANA][NKINDANA][RP::NHAR4];
    TProfile* RPCORCNT_SIN[NDETANA][NKINDANA][RP::NHAR4];
    TProfile* RPCORSMD_COS[3][2][NKINDANA][RP::NHAR4];
    TProfile* RPCORSMD_SIN[3][2][NKINDANA][RP::NHAR4];

    TH2* m_svx_nhits_vs_bbcq;
    //TH3* m_svx_dca2d_dca3d[NCENTANA][NETAANA][2];
    TH3* m_svx_dca2d_dca3d[NCENTANA][NETAANA];
    TH3* m_svx_nhits_cn_quality[NCENTANA][2][8];
    TH3* m_svx_eta_pt_phi[2][NCENTANA][2];
    TH3* m_SvxCntTrk[NCENTANA];
    TH2* m_CNT_Trk[NCENTANA];
    TH2* m_SvxCntTrk_DCA;
    TH2* m_SvxCntTrk_ptphi[NCENTANA];
    
    TProfile* vn_svxcnt_tpro_cos[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][3];
    TProfile* vn_svxcnt_tpro_sin[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][3];
    TProfile** vn_seg_tpro_cos_hit3[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][NETAANA];
    TProfile** vn_seg_tpro_sin_hit3[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][NETAANA];
    TProfile** vn_seg_tpro_cos_hit4[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][NETAANA];
    TProfile** vn_seg_tpro_sin_hit4[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][NETAANA];
    TProfile* RPCORBBC_COS[3][NDETANA][NKINDANA][RP::NHAR4];
    TProfile* RPCORBBC_SIN[3][NDETANA][NKINDANA][RP::NHAR4];

    
    //TTree* m_Tree;    
    float m_BBCRP[NKINDANA][RP::NHAR4];//N
    float m_FVTRP[NKINDANA+3][RP::NHAR4];
    float m_CNTRP[RP::NHAR4];//N
    float m_SMDRP[NKINDANA];//N
    float m_MPCRP[NKINDANA][RP::NHAR4];//N
    float m_phi[1000];
    float m_pt[1000];
    float m_eta[1000];
    int   m_cnt_ntrk;
    int   m_charge[1000];
    int   m_PID[1000];
    int   m_dcarm[1000];

    //Svx SRG Variables
    float m_seg_pt[1000];
    float m_seg_phi[1000];
    float m_seg_eta[1000];
    int   m_seg_hits[1000];
    float m_seg_dca2d[1000];
    float m_seg_dca3d[1000];
    float m_seg_Chisq_Ndf[1000];
    float m_seg_Quality[1000];
    int   m_seg_ntrk;
    
    //Svx CNT Trak
    float m_svxcnt_pt[1000];    
    float m_svxcnt_eta[1000];    
    float m_svxcnt_phi[1000];
    float m_svxcnt_dcaz[1000];
    float m_svxcnt_dca2dr[1000];
    float m_svxcnt_sdca2dr[1000];
    int   m_svxcnt_ntrk;

    float   m_Mean_qx_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4];
    float  m_Sigma_qx_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4];
    float   m_Mean_qy_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4];
    float  m_Sigma_qy_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4];
    float  m_Flat_Cos_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4][RP::NORD3];
    float  m_Flat_Sin_tmp[RP::NMUL3][RP::NZPS2][75][RP::NHAR4][RP::NORD3];

    float m_Dca_Mean[ 7]; 
    float m_Dca_Sigma[7];
    float m_pc3sdphi[1000];
    float m_pc3sdz[  1000];

    //Test
    TProfile* Re_QSMD[  10][3];
    TH2*      Re_Q2DSMD[10][3];
    TProfile* Re_QBBC[  10][3];
    TH2*      Re_Q2DBBC[10][3];
    TProfile* Re_QCNT[  10];
    TH2*      Re_Q2DCNT[10];
    
    TProfile*  vn_tpro_cos_alleta[NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    TProfile*  vn_tpro_sin_alleta[NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    
    int m_NFvtx_Trk;
    float m_fvtx_eta[2000];
    float m_fvtx_phi[2000];
    TH2* Fvtx_Phi_Eta[10];
    TH2* fvtx_xy[10][2];
    TH2* fvtx_xy_DCACUT_x2_y2    [10][12][3][2];
    TH2* fvtx_xy_DCACUT_xcos_ysin[10][12][3][2];
    TH2* m_Fvtx_Eta[3];
    TH2* m_Fvtx_Phi[3][2];
    TProfile* vn_tpro_cos_fvtx[NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    TProfile* vn_tpro_sin_fvtx[NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    TH1* BbcCha[3];
    TH2* BbcCha2D;
    TH2* seg_dphi_dPsi[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][NETA];
    TH2* m_m2tof_e;
    TH2* m_m2tof_w;
    TH2* m_m2tof_e_PID;
    TH2* m_m2tof_w_PID;
    TH2* QA_TOF_dPhi;
    TH2* QA_TOF_dZ;
    TH2* QA_PC3_dPhi;
    TH2* QA_PC3_dZ;
    TH2* QA_TOF_SdPhi;
    TH2* QA_TOF_SdZ;
    TH2* QA_PC3_SdPhi;
    TH2* QA_PC3_SdZ;
    TProfile* vn_tpro_cos_allpt[3][4][10];
    TProfile* vn_tpro_sin_allpt[3][4][10];
    TH2* Fvtx_DCA[10][12][2];
    TH2* Fvtx_DCA_X2_Y2    [NCENTANA][NZANA][2];
    TH2* Fvtx_DCA_XCos_YSin[NCENTANA][NZANA][2];
    TH2*  h_Cnt_Fvt_dPhi_Real[NCENTANA][NZANA][3][5];
    TH2*  h_Cnt_Fvt_dPhi_Mix [NCENTANA][NZANA][3][5];

    TH2* Fvtx_X_eta[NCENTANA][NZANA][2];
    TH2* Fvtx_Y_eta[NCENTANA][NZANA][2];
    TH2* Fvtx_R_eta[NCENTANA][NZANA][2];

    TH1* h_Fvt_S_N_dPhi_Real[NCENTANA][NZANA][5];
    TH1* h_Fvt_S_N_dPhi_Mix [NCENTANA][NZANA][5];

    TH2* dPhi_Real_trk_EP[NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    TH2* dPhi_Mix_trk_EP [NDETANA][NKINDANA][RP::NHAR4][NCENTANA];
    TH2* Real_Re_S_N_COR [NDETANA][2][RP::NHAR4];
    TH2* Mix_Re_S_N_COR  [NDETANA][2][RP::NHAR4];

    TH2* Real_CNTReRP_ReRP_COR[NDETANA][NKINDANA][RP::NHAR4];
    TH2* Mix_CNTReRP_ReRP_COR [NDETANA][NKINDANA][RP::NHAR4];
    TH2* Fvtx_Phi_DcaR_X2_Y2    [NCENTANA][NZANA][3];
    TH2* Fvtx_Phi_DcaR_Xcos_Ysin[NCENTANA][NZANA][3];
    TNtuple* m_Ntple_Fvt_Trk_Info;
    TNtuple* m_Ntple_Fvt_Trk_DCA_Info;
    TH2* h_Cnt_Bbc_dPhi_Real[NCENTANA][2];
    TH2* h_Cnt_Bbc_dPhi_Mix [NCENTANA][2];
    TH1* h_Bbc_S_N_dPhi_Real[NCENTANA];  
    TH1* h_Bbc_S_N_dPhi_Mix [NCENTANA];
    TH2* h_Cnt_Mpc_dPhi_Real[NCENTANA][2];
    TH2* h_Cnt_Mpc_dPhi_Mix [NCENTANA][2];
    TH1* h_Mpc_S_N_dPhi_Real[NCENTANA];  
    TH1* h_Mpc_S_N_dPhi_Mix [NCENTANA];
    TH1* dPhi_Real_FvtTrk_ReRP[NDETANA][NKINDANA][RP::NHAR4][NCENTANA][3];
    TH1* dPhi_Mix_FvtTrk_ReRP [NDETANA][NKINDANA][RP::NHAR4][NCENTANA][3];
    TH2* Nfvtx_S_N[3][10];
    TH2* Nfvtx_BbcCha[3][3];
    TH1* m_Cnt_Pc3dZ  [10];
    TH1* m_Cnt_Pc3dPhi[10];
    TH2* m_bbcq;

  TH2* h_Mpc_E_vs_Adc[2];
  TH2* h_Tower_E;
  TH2* h_Tower_Val;
  //TH2* h_Tower_Adc;
  TH2* h_BbcV_vs_FvtV[10];
  TH2* h_FvtV_Z_X_Y  [10][3];
  TH2* h_Chi2_R      [10][4][6][2];
  TH2* h_Phi_R       [10][4][6][2];
  TH2* h_Phi_Eta     [10][4][2];
  TH2* h_Hit_Hitp    [10][4][2];
  TH2* h_Hit_R       [10][4][6][2];
  TH2* h_Hit_Phi     [10][4][3];
  TH3* Nfvtx_Eta_Phi_Vertex[10];
  TH2* Nfvtx_Hits_Vertex[10][2];

  float m_NFvtx_Trk_S;
  float m_NFvtx_Trk_N;;
  float m_cnt_ntrk_good;
  float m_cnt_ntrk_good_Q;
  float m_cnt_nemc_good;
  float m_cnt_nemc_good_Q;
  TH2* h_NBbc_vs_NFvt   [3][3];
  // S vs N ( Bbc, Fvtx )
  TH2* h_NBbc_S_vs_N    [3];
  TH2* h_NFvt_S_vs_N    [3];
  // Cnt vs Fvtx 
  TH2* h_NCntTrk_vs_NFvt [3][3];
  TH2* h_NCntTrkQ_vs_NFvt[3][3];
  TH2* h_NCntEmc_vs_NFvt [3][3];
  TH2* h_NCntEmcQ_vs_NFvt[3][3];
  TH2* h_Fvt_phi_eta;
	*/

};

#endif
