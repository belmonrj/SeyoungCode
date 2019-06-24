#define M_PI 3.14159265358979323846
#include <stdio.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <math.h>
#include <TSpectrum.h>
#include "def_Const.h"
#include <TProfile2D.h>

int	nBbcPmt;
short	BbcPmt_arm[128];
float	BbcPmt_time[128];
float	BbcPmt_phi[128];
float	BbcPmt_val[128];
float	BbcPmt_x[128];
float	BbcPmt_y[128];
float	BbcPmt_z[128];
float	BbcPmt_eta[128];

int	nFvtxTrack;
short	FvtxTrack_nhits[1024];
float	FvtxTrack_phi[1024];
float	FvtxTrack_eta[1024];
float	FvtxTrack_chi2[1024];
float	FvtxTrack_dca_x[1024];
float   FvtxTrack_dca_y[1024];
float	FvtxTrack_dca_r[1024];
short	FvtxTrack_has_svxhit[1024];

int	nCntTrack;
float	CntTrack_phi[1024];
float	CntTrack_eta[1024];
float	CntTrack_pt[1024];
//float	CntTrack_phi1[1024];
//float	CntTrack_phi3[1024];
//float	CntTrack_zpc1[1024];
//float	CntTrack_zpc3[1024];
float	CntTrack_mom[1024];
short	CntTrack_ch[1024];
short	CntTrack_arm[1024];
float	CntTrack_pc3sdz[1024];
float	CntTrack_pc3sdphi[1024];

int	runnumber, trigbit_scaled, eventnumber;
short	centrality;
float	bbc_z;
float	fvtx_x, fvtx_y, fvtx_z;
float	beam_x, beam_y;
float	doubleint_frac;
int	doubleint_good;

typedef struct GetFVTXtrack
{
  int rawTrack;
  int goodTrack[NARM+1];
  int hits[1000];
  int arm[1000];
  float phi[1000];
  float eta[1000];
  float chi2[1000];
  float dca_x[1000];
  float dca_y[1000];
  float dca_r[1000];
  int svxhit[1000];

  GetFVTXtrack(){
	rawTrack = 0;
	goodTrack[0] = goodTrack[1] = goodTrack[2] = 0;
	for (int ii=0; ii<1000; ii++){
		hits[ii] = arm[ii] = 0;
		phi[ii] = eta[ii] = chi2[ii] = 0;
		dca_x[ii] = dca_y[ii] = dca_r[ii] = 0.0;
		svxhit[ii] = 0;
	}
  }//

} FVTX_INFO;

typedef struct GetBBCpmt
{
  int rawPmt;
  int arm[128];
  float time[128];
  float phi[128];
  float cha[128];
  float eta[128];

  GetBBCpmt(){
	rawPmt = 0;
	for (int ii=0; ii<128; ii++){
		arm[ii] = 0;
		time[ii] = phi[ii] = cha[ii] = eta[ii] = 0;
	}
  }//
} BBC_INFO;

typedef struct GetCNTtrack
{
  int track;
  int ipt_bin[500];
  float phi[500];
  float eta[500];
  float pt[500];
  //float phi1[500];
  //float phi3[500];
  //float zpc1[500];
  //float zpc3[500];
  float	mom[500];
  int	ch[500];
  int	arm[500];
  float	pc3sdz[500];
  float	pc3sdphi[500];

  GetCNTtrack(){
	for (int ii=0; ii<500; ii++){
		ipt_bin[ii] = ch[ii] = arm[ii] = 0;
		//phi[ii] = eta[ii] = pt[ii] = phi1[ii] = phi3[ii] = 0.0;
		//zpc1[ii] = zpc3[ii] = mom[ii] = pc3dz[ii] = pc3dphi[ii] = 0.0;
		phi[ii] = eta[ii] = 0.0;
		pt[ii] = mom[ii] = pc3sdz[ii] = pc3sdphi[ii] = 0.0;
	}
  }//

} CNT_INFO;

typedef struct GetRunInfo
{
  int Run;
  float bbc_vertex_z;
  float fvtx_vertex_x;
  float fvtx_vertex_y;
  float fvtx_vertex_z;
  float beam_offset_x;
  float beam_offset_y;

  GetRunInfo(){
	Run = -9999;
	bbc_vertex_z = -9999;
	fvtx_vertex_x = fvtx_vertex_y = fvtx_vertex_z = -9999;
  }
} RUN_INFO;


struct EVNT
{
  FVTX_INFO	FVTX;
  BBC_INFO	BBC;
  CNT_INFO	CNT;
  RUN_INFO	RUN;
 
  float vertex;
  int cent;  //centrality bin
  int zvtx;  //z vertex bin
  int npeaks[NARM];

  EVNT(){
	vertex = -9999;
	cent = zvtx = -9999;
	npeaks[0] = npeaks[1] = 0;
  }
};


// Set histograms
TH1*  h_Bbc_S_N_dPhi_Real         [NCENTANA];
TH1*  h_Bbc_S_N_dPhi_Mix          [NCENTANA];
TH1*  h_Fvt_S_N_dPhi_Real         [NCENTANA][NZANA];
TH1*  h_Fvt_S_N_dPhi_Mix          [NCENTANA][NZANA];

//TH2*  h_Cnt_Cnt_dPhi_dEta_Real    [NCENTANA][NPTBIN][NPTBIN];
//TH2*  h_Cnt_Cnt_dPhi_dEta_Mix     [NCENTANA][NPTBIN][NPTBIN];
TH2*  h_Cnt_Bbc_dPhi_Real         [NCENTANA][NARM];
TH2*  h_Cnt_Bbc_dPhi_Mix          [NCENTANA][NARM];
TH2*  h_Cnt_Bbc_dPhi_dEta_Real    [NCENTANA][NPTBIN];
TH2*  h_Cnt_Bbc_dPhi_dEta_Real_sqr[NCENTANA][NPTBIN];
TH2*  h_Cnt_Bbc_dPhi_dEta_Mix     [NCENTANA][NPTBIN];
TH2*  h_Cnt_Bbc_dPhi_dEta_Mix_sqr [NCENTANA][NPTBIN];
TH2*  h_Bbc_S_N_dPhi_dEta_Real    [NCENTANA];
TH2*  h_Bbc_S_N_dPhi_dEta_Real_sqr[NCENTANA];
TH2*  h_Bbc_S_N_dPhi_dEta_Mix     [NCENTANA];
TH2*  h_Bbc_S_N_dPhi_dEta_Mix_sqr [NCENTANA];
TH2*  h_Cnt_Fvt_dPhi_Real         [NCENTANA][NZANA][NARM];
TH2*  h_Cnt_Fvt_dPhi_Mix          [NCENTANA][NZANA][NARM];
TH2*  h_Fvt_S_N_dPhi_dEta_Real    [NCENTANA][NZANA];
TH2*  h_Fvt_S_N_dPhi_dEta_Mix     [NCENTANA][NZANA];
TH2*  h_Cnt_Fvt_dPhi_dEta_Real    [NCENTANA][NZANA][NPTBIN];
TH2*  h_Cnt_Fvt_dPhi_dEta_Mix     [NCENTANA][NZANA][NPTBIN];
TH2*  h_Bbc_Fvt_S_dPhi_dEta_Real  [NCENTANA][NZANA];
TH2*  h_Bbc_Fvt_S_dPhi_dEta_Mix	  [NCENTANA][NZANA];
TH2*  h_Bbc_Fvt_N_dPhi_dEta_Real  [NCENTANA][NZANA];
TH2*  h_Bbc_Fvt_N_dPhi_dEta_Mix	  [NCENTANA][NZANA];
TH2*  h_BbcS_FvtN_dPhi_dEta_Real  [NCENTANA][NZANA];
TH2*  h_BbcS_FvtN_dPhi_dEta_Mix   [NCENTANA][NZANA];
TH2*  h_BbcN_FvtS_dPhi_dEta_Real  [NCENTANA][NZANA];
TH2*  h_BbcN_FvtS_dPhi_dEta_Mix   [NCENTANA][NZANA];

//QA histos
TH1*  h_Fvt_trk	[NCENTANA];
TH1*  h_Trk_eta [10][2];
TH1*  h_Bbc_time[NARM];
TH1*  h_Cnt_pt	[NCENTANA];
TH1*  h_Bbc_cha	[NCENTANA];
TH2*  h_Fvt_x_z;
TH2*  h_Fvt_y_z;
TSpectrum*	spec[NARM];
TH2*  h_Fvt_Bbc	[NCENTANA];
TH1*  h_Fvt_dca	[NCENTANA][NARM];
TH2*  h_Fvt_eta_dca			[NCENTANA][NARM];
TSpectrum*	h_Fvt_trk_phi	[NARM];
TProfile2D* h_Bbc_cha_x_y	[NARM];
TH1*  h_Bbc_cha_dist		[128];
TH2*  h_Bbc_index			[NARM];

// Set event buffer variable
int 	n_Eve_buff	[NCENTANA][NZANA];
int	b_fvtx_trk	[NCENTANA][NZANA][NEVE];
int	b_fvtx_arm	[NCENTANA][NZANA][NEVE][1024];
float	b_fvtx_phi	[NCENTANA][NZANA][NEVE][1024];
float	b_fvtx_eta	[NCENTANA][NZANA][NEVE][1024];
int	b_fvtx_svxhit	[NCENTANA][NZANA][NEVE][1024];

int	b_cnt_trk	[NCENTANA][NZANA][NEVE];
int	b_cnt_pt_bin	[NCENTANA][NZANA][NEVE][1024];
float	b_cnt_phi	[NCENTANA][NZANA][NEVE][1024];
float	b_cnt_pt	[NCENTANA][NZANA][NEVE][1024];
float	b_cnt_eta	[NCENTANA][NZANA][NEVE][1024];
//float	b_cnt_phi1	[NCENTANA][NZANA][NEVE][1024];
//float	b_cnt_phi3	[NCENTANA][NZANA][NEVE][1024];
//float	b_cnt_zpc1	[NCENTANA][NZANA][NEVE][1024];
//float	b_cnt_zpc3	[NCENTANA][NZANA][NEVE][1024];

int	b_bbc_pmt	[NCENTANA][NZANA][NEVE];
int	b_bbc_arm	[NCENTANA][NZANA][NEVE][1024];
float	b_bbc_phi	[NCENTANA][NZANA][NEVE][1024];
float	b_bbc_eta	[NCENTANA][NZANA][NEVE][1024];
float	b_bbc_cha	[NCENTANA][NZANA][NEVE][1024];

