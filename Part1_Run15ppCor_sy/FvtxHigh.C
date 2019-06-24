// Quality assurance 
// CNT  : pt, phi, eta ,, N tracks,, 
// BBC  : Total charge, each pmt charge ....eta, phi, 
// FVTX : phi, eta, N track .......
//////////////////////////////////////////////////////
#include "FvtxHigh.h"

#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>
#include <RunNumberRanges.h>
#include <getClass.h>
#include <recoConsts.h>

#include <PHGlobal.h>
#include <PreviousEvent.h>
#include <EventHeader.h>
#include <PHCentralTrack.h>
#include <PHSnglCentralTrack.h>

#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "MpcMap.h"
#include "mpcTowerContainer.h"
#include "mpcTowerContent.h"
//#include "mpcRawContainer.h"
//#include "mpcRawContent.h"

#include "TrigLvl1.h"
//Reaction Plane Object 
#include "RpSumXYObject.h"
#include "RpSnglSumXY.h"
#include "ReactionPlaneObject.h"
#include "ReactionPlaneSngl.h"
#include "ReactionPlaneCalibv1.h"

//BBC PMT Info.
#include <Bbc.hh>
#include <BbcRaw.h>
#include <BbcCalib.hh>
#include <BbcGeo.hh>
#include <RunToTime.hh>

#include <TOAD.h>

#include <TFvtxCompactTrkMap.h>
#include <PHPoint.h>

//SVX Track 
#include <compactCNT/SvxCentralTrackMap.h>
#include <compactCNT/SvxCentralTrackMapEntry.h>
#include "SvxClusterList.h"
#include "SvxCluster.h"
#include <compactCNT/SvxTrackMapEntry.h>
#include <compactCNT/SvxTrackMap.h>
#include "SvxSegmentList.h"
#include "SvxSegment.h"
#include <SvxClusterInfo.h>
#include <SvxCentralTrackList.h>
#include <SvxCentralTrack.h>
#include <SvxSegmentList.h>

#include <VtxOut.h>

//ROOT Hists
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TNtuple.h>

//#include "CuAu200_TOF_matching.h"
//#include "CuAu200_PC3_matching.h"
//#include "EmcPC3Matching.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

using namespace std;

const float MomfactorWest = 0.9742;
const float MomfactorEast = 0.9727;

const float XOffsetE = 1.19326e-01;
const float YOffsetE = 1.30366e-01;

const float XOffsetW = 1.81180e-01;
const float YOffsetW = 1.48351e-01;
//==============================================================
FvtxHigh::FvtxHigh(const char* output) : 
  SubsysReco("RPCALIBRATOR")
  //,m_Ana_Flag(0)
  //,rpcalibv1(NULL)  
  //,HistoManager(new Fun4AllHistoManager("RPCalibHisto"))
  //,file(NULL)
  //,bbccalib(NULL)  
  //,bbcgeo(NULL)
  //,OutputFileName(output)
{
	/*
  //outNam = (char*)output; 
  m_RunNumber   = 0;
  m_EventNumber = 0;
  // histo NULL
  m_centzv = NULL;
  m_centzv_svx = NULL;
  m_SvxCntTrk_DCA = NULL;
  m_bbcgainflag = 0;
  
  //CNT
  fill_n(m_pt,      sizeof(m_pt)/sizeof(m_pt[0]), -9999.);
  fill_n(m_phi,     sizeof(m_phi)/sizeof(m_phi[0]), -9999.);
  fill_n(m_eta,     sizeof(m_eta)/sizeof(m_eta[0]), -9999.);
  fill_n(m_pc3sdz,  sizeof(m_pc3sdz)/sizeof(m_pc3sdz[0]), -9999);
  fill_n(m_pc3sdphi,sizeof(m_pc3sdphi)/sizeof(m_pc3sdphi[0]),-9999);
  m_cnt_ntrk = 0.;
  
  //Svx CNT
  fill_n(m_pt,           sizeof(m_pt)/sizeof(m_pt[0]), -9999.);
  fill_n(m_svxcnt_pt    ,sizeof(m_svxcnt_pt)    /sizeof(    m_svxcnt_pt[0]),-9999.);
  fill_n(m_svxcnt_phi   ,sizeof(m_svxcnt_phi)   /sizeof(   m_svxcnt_phi[0]),-9999.);
  fill_n(m_svxcnt_eta   ,sizeof(m_svxcnt_eta)   /sizeof(   m_svxcnt_eta[0]),-9999.);
  fill_n(m_svxcnt_dcaz  ,sizeof(m_svxcnt_dcaz)  /sizeof(  m_svxcnt_dcaz[0]),-9999.);
  fill_n(m_svxcnt_dca2dr,sizeof(m_svxcnt_dca2dr)/sizeof(m_svxcnt_dca2dr[0]),-9999.);
  fill_n(m_svxcnt_sdca2dr,sizeof(m_svxcnt_sdca2dr)/sizeof(m_svxcnt_sdca2dr[0]),-9999.);
  m_svxcnt_ntrk=0;
	*/
  
//  for(int ic_ana=0;ic_ana<10;ic_ana++)
//    for(int iz_ana=0;iz_ana<12;iz_ana++)
//      do_Buff_Cnt_Fvt_Trk_Reset(ic_ana,iz_ana);
}
//==============================================================
FvtxHigh::~FvtxHigh()
{
	/*
  delete bbccalib;
  delete bbcgeo;
  delete HistoManager;
  //delete file;
	*/
}
/*
//==============================================================
int FvtxHigh::Init(PHCompositeNode *topNode)
{
  cout << "FvtxHigh::Init started..." << endl;
  
  bbccalib = new BbcCalib();
  bbcgeo   = new BbcGeo();
  //file     = new TFile(outNam,"recreate");  

  cout << "FvtxHigh::Init ended." << endl;
  return 0;
}
//==============================================================
int FvtxHigh::InitRun(PHCompositeNode *topNode)
{
  cout << "FvtxHigh::InitRun started..." << endl;
  
  // Get Run Number
  recoConsts *rc = recoConsts::instance();
  m_RunNumber = rc->get_IntFlag("RUNNUMBER");  
  
  ///Calib bbc pmt
  int icalibversion = 4002;// Calibrated PMTs are used  
  RunToTime* runTime = RunToTime::instance();
  PHTimeStamp* ts( runTime->getBeginTime(m_RunNumber) );
  PHTimeStamp tstart = *ts;
  delete ts;
  int bbccalib_version = 4002;
  cout << "FvtxHigh::InitRun, run number= " << m_RunNumber
       << " Version of BbcCalib = " << bbccalib_version << endl;
  bbccalib->restore(tstart, bbccalib_version);
  rc->set_IntFlag("BBCCALIBVERSION", icalibversion);  
  
  cout << "FvtxHigh::InitRun   run number      = " << m_RunNumber << endl;
  
  //histogram Initialization
  initHisto();
  
  
  // this must be after getting rpcalib
  cout << "FvtxHigh::InitRun ended." << endl;
  
  return 0;    
}
//==============================================================
int FvtxHigh::process_event(PHCompositeNode *topNode)
{
  //cout << "//-------------------------- process-------------------------------- " << m_EventNumber << endl;
  if(m_EventNumber%10000==0) cout << m_EventNumber<< "//------------------------------------" << endl;
  m_EventNumber++;

  PHGlobal                *global = findNode::getClass<PHGlobal>(     topNode, "PHGlobal");
  PreviousEvent *peve             = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent"); // Tick
  //RpSumXYObject *sumxy  = findNode::getClass<RpSumXYObject>(topNode, "RpSumXYObject");
  //ReactionPlaneObject        *rp  = findNode::getClass<ReactionPlaneObject>(topNode, "ReactionPlaneObject"); 
  PHCentralTrack         *central = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  BbcRaw                 *bbcraw  = findNode::getClass<BbcRaw>(topNode,"BbcRaw");  
  TFvtxCompactTrkMap* trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode,"TFvtxCompactTrkMap");
  //VtxOut                  *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut");
  TrigLvl1 *triglvl1              = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  //mpcTowerContainer     *mpctower = findNode::getClass<mpcTowerContainer>(topNode,"mpcTowerContainer");
  //mpcRawContainer   *mpc_Rawtower = findNode::getClass<mpcRawContainer>(topNode,"MpcRaw");
  //if (!mpc_Rawtower) {
  //cout << PHWHERE << "could not get mpcRawContainer, is Node missing?" << endl;
  //return EVENT_OK;
  //}
  
  // find SvxSegment
  //SvxSegmentList *svxseglist = findNode::getClass<SvxSegmentList>(topNode, "SvxSegmentList");
  //Svx Central Track
  //SvxCentralTrackList *d_svxcnt = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
  //if(d_svxcnt == NULL) { cerr << PHWHERE << " SvxCentralTrackList node not found." << endl; return DISCARDEVENT;}
  // Svx Cluster
  //SvxClusterList *d_svxhit = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  //if(d_svxhit==NULL) { cerr << PHWHERE << " SvxClusterList node not found." << endl; return DISCARDEVENT; }    
  
  //SvxTrackMap* d_svxtrk = findNode::getClass<SvxTrackMap>(topNode, "SvxTrack_comp");  
  if(!global){ cerr << PHWHERE << " No PHGloval object !"      << endl; return EVENT_OK; }
  if(!peve)  { cerr << PHWHERE << " No PreviousEvent object !" << endl; return EVENT_OK; }
  //if(!sumxy) { cerr << PHWHERE << " No RpSumXYObject object !" << endl; return EVENT_OK; }
  if(!central){cerr << PHWHERE << " CNT Event skip "        << endl;  return EVENT_OK; }    
  //else { rp = NULL; }
  //if(!vtxout) {cerr << PHWHERE << " VTXOUT Event skip "     << endl;  return EVENT_OK; }    
  //if (!mpctower) { cout << PHWHERE << "could not get mpcTowerContainer, is Node missing?" << endl;    return EVENT_OK; }  
  if (!trkfvtx_map) {cerr << PHWHERE << " No TFvtxCompactTrkMap object !" << endl;  return EVENT_OK; }
  if (!bbcraw){cerr << PHWHERE << "Could not find BbcRaw !" << endl;  return EVENT_OK; }
  
  //  //Trigger selection
  if(!triglvl1)
    {
      cout << PHWHERE << "Could not find TrigLvl1 !" << endl;
      return 1;
    }
  


  int Scaled = triglvl1->get_lvl1_trigscaled();
  int itrig=0;                       // MB
  if(Scaled & 0x00000400 ) itrig = 1;// Fvtx High Mul And
  if(Scaled & 0x00000800 ) itrig = 2;// Fvtx high Mul Or 
  
  //if ((Scaled & 0x00000002) <= 0) return 1;//BBCLL1(>1 tubes) narrow vertex
  //int bbc0Scale = 0;
  //if ((Scaled & 0x00000004) > 0 || //BBCLL1(>1 tubes) narrowvtx central
  //    (Scaled & 0x00000008) > 0) //BBCLL1(>1 tubes) narrowvtx central
  //  //(Scaled & 0x00000010) > 0) //BBCLL1(>1 tubes) narrowvtx central
  //  bbc0Scale = 2;
  //if(bbc0Scale!=2) return 1;
  
  // Z_Vertex Cut //
  float vertex = -9999;
  int  iz_ana  = -9999;
  vertex = global->getBbcZVertex(); 
  if(fabs(vertex)>10) return EVENT_OK;
  if     (vertex<-5) iz_ana=0;
  else if(vertex<0 ) iz_ana=1;
  else if(vertex<5 ) iz_ana=2;
  else if(vertex<10) iz_ana=3;
  
  //float vertex_fvt = vtxout->get_Vertex("FVTX").getZ();
  float bbcq = global->getBbcChargeN()+global->getBbcChargeS();
  float bbcs = global->getBbcChargeS();
  float bbcn = global->getBbcChargeN();
  if(bbcq<0) return EVENT_OK;
  //int ic_ana=-9999;
  //int ic_ana=0;
  
  // event histo
  m_bbcq  ->Fill(bbcq,vertex);
  //m_centzv->Fill(centrality, vertex);
  
  //int ieve = m_Eve[ic_ana][iz_ana];
  
  // Reset 
  do_Buff_Cnt_Fvt_Trk_Reset();  
  
  // Track Check
  GetCNTPC3_TOF(central);                           // Get CNT  information
  GetFVTXTrack (global, trkfvtx_map, iz_ana, itrig);// Get FVTX information
  
  // Bbc vs Fvtx 
  h_NBbc_vs_NFvt   [0][itrig]->Fill(bbcs, m_NFvtx_Trk_S              );
  h_NBbc_vs_NFvt   [1][itrig]->Fill(bbcn, m_NFvtx_Trk_N              );
  h_NBbc_vs_NFvt   [2][itrig]->Fill(bbcq, m_NFvtx_Trk_S+m_NFvtx_Trk_N);
  // S vs N ( Bbc, Fvtx )
  h_NBbc_S_vs_N       [itrig]->Fill(bbcs, bbcn                       );
  h_NFvt_S_vs_N       [itrig]->Fill(m_NFvtx_Trk_S,m_NFvtx_Trk_N      );
  // Cnt vs Fvtx 
  h_NCntTrk_vs_NFvt [0][itrig]->Fill(m_cnt_ntrk_good  , m_NFvtx_Trk_S              );
  h_NCntTrk_vs_NFvt [1][itrig]->Fill(m_cnt_ntrk_good  , m_NFvtx_Trk_N              );
  h_NCntTrk_vs_NFvt [2][itrig]->Fill(m_cnt_ntrk_good  , m_NFvtx_Trk_S+m_NFvtx_Trk_N);
  h_NCntTrkQ_vs_NFvt[0][itrig]->Fill(m_cnt_ntrk_good_Q, m_NFvtx_Trk_S              );
  h_NCntTrkQ_vs_NFvt[1][itrig]->Fill(m_cnt_ntrk_good_Q, m_NFvtx_Trk_N              );
  h_NCntTrkQ_vs_NFvt[2][itrig]->Fill(m_cnt_ntrk_good_Q, m_NFvtx_Trk_S+m_NFvtx_Trk_N);
  h_NCntEmc_vs_NFvt [0][itrig]->Fill(m_cnt_nemc_good  , m_NFvtx_Trk_S              );
  h_NCntEmc_vs_NFvt [1][itrig]->Fill(m_cnt_nemc_good  , m_NFvtx_Trk_N              );
  h_NCntEmc_vs_NFvt [2][itrig]->Fill(m_cnt_nemc_good  , m_NFvtx_Trk_S+m_NFvtx_Trk_N);
  h_NCntEmcQ_vs_NFvt[0][itrig]->Fill(m_cnt_nemc_good_Q, m_NFvtx_Trk_S              );
  h_NCntEmcQ_vs_NFvt[1][itrig]->Fill(m_cnt_nemc_good_Q, m_NFvtx_Trk_N              );
  h_NCntEmcQ_vs_NFvt[2][itrig]->Fill(m_cnt_nemc_good_Q, m_NFvtx_Trk_S+m_NFvtx_Trk_N);
  
//  // Get Tracks, pmts, towers, for 2PC method
//  if( m_Ana_Flag == FLAG_FVTX_2PC)
//    {
//      GetFVTXTrack (ic_ana, vertex, ieve, global, trkfvtx_map, vtxout);
//      GetCNTPC3_TOF(ic_ana, vertex, ieve, central);
//    }
//  if(m_Ana_Flag == FLAG_2PC_METHOD)
//    {
//      GetCNTPC3_TOF(ic_ana, vertex, ieve, central);
//      GetBBCPMT    (ic_ana, iz_ana, ieve, bbcraw);
//      GetMPCTower  (ic_ana, iz_ana, vertex, ieve, mpctower);
//    }
  
//  // Get Q-Vector from RpObject 
//  if(m_Ana_Flag == FLAG_EP_METHOD)
//    {
//      if (rp) 
//	{
//	  fill_n(&m_BBCRP[0][0],    sizeof(m_BBCRP)/sizeof(float), -9999.);
//	  fill_n(&m_MPCRP[0][0],    sizeof(m_MPCRP)/sizeof(float), -9999.);
//	  fill_n(&m_CNTRP[0],       sizeof(m_CNTRP)/sizeof(float), -9999.);
//
//	  GetAllRP(RP::ID_BBC, icent, iz_ana, rp);
//	  GetAllRP(RP::ID_MPC, icent, iz_ana, rp);
//	  GetAllRP(RP::ID_CNT, icent, iz_ana, rp);
//	  
//	  //fill vn
//	  do_vn(RP::ID_BBC, ic_ana, iz_ana, 0);
//	  do_vn(RP::ID_MPC, ic_ana, iz_ana, 0);
//	  do_vn(RP::ID_CNT, ic_ana, iz_ana, 0);
//	  //if( m_Vn == FLAG_PID)
//	  //{
//	  //do_vn(RP::ID_BBC, ic_ana, 1);
//	  //do_vn(RP::ID_SMD, ic_ana, 1);
//	  //}
//	  
//	  //fill rp correlations
//	  do_RP_SN_correlation ( RP::ID_BBC, 0, ic_ana, iz_ana);
//	  do_RP_CNT_correlation( RP::ID_BBC,    ic_ana, iz_ana);
//	  do_RP_SN_correlation ( RP::ID_MPC, 0, ic_ana, iz_ana);
//	  do_RP_CNT_correlation( RP::ID_MPC,    ic_ana, iz_ana);
//	  //do_RP_SMD_correlation( RP::ID_BBC,    ic_ana, iz_ana);
//	  
//	  //if(3<iz_ana && iz_ana<8) do_RP_SN_correlation ( RP::ID_FVT, 0, ic_ana, iz_ana);
//	  //if(3<iz_ana && iz_ana<8) do_RP_SN_correlation ( RP::ID_FVT, 1, ic_ana, iz_ana);
//	  //if(3<iz_ana && iz_ana<8) do_RP_CNT_correlation( RP::ID_FVT,    ic_ana, iz_ana);
//	  //if(3<iz_ana && iz_ana<8) do_RP_SMD_correlation( RP::ID_FVT,    ic_ana, iz_ana);
//	}
//    } //if(rp)
//  
//  if(m_Ana_Flag == FLAG_2PC_METHOD)
//    {
//      //if(rp)
//      //{
//      //GetAllRP(RP::ID_SMD, icent, iz_ana, rp);
//      //GetAllRP(RP::ID_CNT, icent, iz_ana, rp);
//      //GetAllRP(RP::ID_FVT, icent, iz_ana, rp);
//      
//      //2PC CNT-FVTX
//m_Eve[ic_ana][iz_ana]++;
//      if(m_Eve[ic_ana][iz_ana]==NEVE)
//	    {
//	      //2PC method CNT vs BBC
//	      do_Real_Mix_CNT_BBC_Pmt_Cor(ic_ana,iz_ana);
//	      do_Real_Mix_BBC_S_N_Pmt_Cor(ic_ana,iz_ana);	      
//	      
//	      //2PC method CNT vs MPC
//	      do_Real_Mix_CNT_MPC_Tow_Cor(ic_ana,iz_ana);
//	      do_Real_Mix_MPC_S_N_Tow_Cor(ic_ana,iz_ana);
//	      
//	      //2PC method FVT vs SMD
//	      //do_Real_Mix_FVT_Trk_ReRP_Cor(RP::ID_SMD, ic_ana, iz_ana);	      
//
//	      //EP method CNT vs SMD (FVT)
//	      //do_Real_Mix_CNT_Trk_ReRP_Cor( RP::ID_SMD, ic_ana, iz_ana);
//	      //do_Real_Mix_CNT_Trk_ReRP_Cor( RP::ID_FVT, ic_ana, iz_ana);
//	      
//	      //EP method S vs N cor (FVT, SMD)
//	      //do_Real_Mix_ReRP_SN_Cor( RP::ID_SMD, 0, ic_ana, iz_ana);
//	      //do_Real_Mix_ReRP_SN_Cor( RP::ID_FVT, 0, ic_ana, iz_ana);
//	      //do_Real_Mix_ReRP_SN_Cor( RP::ID_FVT, 1, ic_ana, iz_ana);
//	      
//	      //do_Real_Mix_CNTReRP_ReRPCor( RP::ID_FVT, ic_ana, iz_ana);
//	      
//	      do_Buff_Cnt_Fvt_Trk_Reset(ic_ana,iz_ana);
//	    }     
//	  //}
//    }
//  
//  if(m_Ana_Flag == FLAG_FVTX_2PC)
//    {
//      m_Eve[ic_ana][iz_ana]++;
//      if(m_Eve[ic_ana][iz_ana]==NEVE)
//	{
//	  //2PC method CNT vs FVT
//	  do_Real_Mix_CNT_FVT_Trk_Cor(ic_ana,iz_ana);
//	  do_Real_Mix_FVT_S_N_Trk_Cor(ic_ana,iz_ana);
//	  
//	  do_Buff_Cnt_Fvt_Trk_Reset(ic_ana,iz_ana);
//	}
//    }
  
  return 0;
}
//==============================================================
int FvtxHigh::End(PHCompositeNode *topNode)
{
  cout << "FvtxHigh::End" << endl;
  
  //Write_Histo();
  //file->Write();
  //file->Close();
  HistoManager->dumpHistos(OutputFileName);
  return 0; 
}

bool FvtxHigh::GetFVTXTrack(PHGlobal* global, TFvtxCompactTrkMap *trkfvtx_map, const int iz_ana, const int itrig)
{
  
  //Initialize
  for(int i_trk=0;i_trk<2000;i_trk++)
    {
      m_fvtx_eta[i_trk] =-9999.;
      m_fvtx_phi[i_trk] =-9999.;
    }
  m_NFvtx_Trk=0;

  int itrk_S=0;
  int itrk_N=0;
  TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
  while (TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next()) {
    TFvtxCompactTrk* fvtx_trk = trk_ptr->get();
    float the = fvtx_trk->get_fvtx_theta();
    float eta = fvtx_trk->get_fvtx_eta();
    float phi = fvtx_trk->get_fvtx_phi();
    //float hit = (float)fvtx_trk->get_nhits();
    int   arm = (int)fvtx_trk->get_arm();
    float fvx = fvtx_trk->get_fvtx_vtx().getX();
    float fvy = fvtx_trk->get_fvtx_vtx().getY();
    float fvz = fvtx_trk->get_fvtx_vtx().getZ();
    //short short_chi2 = fvtx_trk->get_short_chi2_ndf();
    //float chi2 = fvtx_trk->get_chi2_ndf();
    //int   hit_p = fvtx_trk->get_hit_pattern();
    if( the==0 && phi==0 && fvx==0 && fvy==0 && fvz==0           ) continue;
    if( (arm==0 && eta>0) || (arm==1 && eta<0) || arm<0 || arm>1 ) continue;
    
    //float DCA_x_BbV                = fvx + tan(the)*cos(phi)*(vertex        - fvz)                ;
    //float DCA_y_BbV                = fvy + tan(the)*sin(phi)*(vertex        - fvz)                ;
    //float DCA_x_BbV_FvX            = fvx + tan(the)*cos(phi)*(vertex        - fvz) - Zvertex_Fvt_X;
    //float DCA_y_BbV_FvY            = fvy + tan(the)*sin(phi)*(vertex        - fvz) - Zvertex_Fvt_Y;
    //float DCA_x_FvV_FvX            = fvx + tan(the)*cos(phi)*(Zvertex_Fvt_Z - fvz) - Zvertex_Fvt_X;
    //float DCA_y_FvV_FvY            = fvy + tan(the)*sin(phi)*(Zvertex_Fvt_Z - fvz) - Zvertex_Fvt_Y;
    //float DCA_R_BbV_x2_y2          = sqrt((DCA_x_BbV*DCA_x_BbV) + (DCA_y_BbV*DCA_y_BbV));
    //float DCA_R_BbV_xcos_ysin      = DCA_x_BbV*cos(phi)+DCA_y_BbV*sin(phi);
    //float DCA_R_BbV_FvXY_x2_y2     = sqrt((DCA_x_BbV_FvX*DCA_x_BbV_FvX) + (DCA_y_BbV_FvY*DCA_y_BbV_FvY));
    //float DCA_R_BbV_FvXY_xcos_ysin = DCA_x_BbV_FvX*cos(phi)+DCA_y_BbV_FvY*sin(phi);
    //float DCA_R_FvV_FvXY_x2_y2     = sqrt((DCA_x_FvV_FvX*DCA_x_FvV_FvX) + (DCA_y_FvV_FvY*DCA_y_FvV_FvY));
    //float DCA_R_FvV_FvXY_xcos_ysin = DCA_x_FvV_FvX*cos(phi)+DCA_y_FvV_FvY*sin(phi);
    
    //m_fvtx_phi[itrk] = phi;	       
    //m_fvtx_eta[itrk] = eta;		   
    //m_Fvt_Phi [ic_ana][ieve][itrk] = phi;
    //m_Fvt_Eta [ic_ana][ieve][itrk] = eta;
    //m_Fvt_DcaX[ic_ana][ieve][itrk] = DCA_x_BbV;
    //m_Fvt_DcaY[ic_ana][ieve][itrk] = DCA_y_BbV;
    //m_Fvt_DcaX_pre[ic_ana][ieve][itrk] = DCA_x_pre;
    //m_Fvt_DcaY_pre[ic_ana][ieve][itrk] = DCA_y_pre;	
    
    if(eta<0) ++itrk_S;
    else      ++itrk_N;
    h_Fvt_phi_eta->Fill(phi,eta);
  } // while
  m_NFvtx_Trk_S = itrk_S;
  m_NFvtx_Trk_N = itrk_N;
  //m_Fvt_Ntrk[ic_ana][ieve] = itrk;  
  
  bool status = true;
  return status;
}

bool FvtxHigh::GetCNTPC3_TOF(PHCentralTrack* central)
{
  
  fill_n(m_pt, sizeof(m_pt)/sizeof(m_pt[0]), -9999.);
  fill_n(m_phi,sizeof(m_phi)/sizeof(m_phi[0]), -9999.);
  fill_n(m_eta,sizeof(m_eta)/sizeof(m_eta[0]), -9999.);
  fill_n(m_pc3sdz, sizeof(m_pc3sdz)/sizeof(m_pc3sdz[0]), -9999);
  fill_n(m_pc3sdphi,sizeof(m_pc3sdphi)/sizeof(m_pc3sdphi[0]),-9999);
  fill_n(m_PID, sizeof(m_PID)/sizeof(m_PID[0]),-9999);
  fill_n(m_dcarm, sizeof(m_dcarm)/sizeof(m_dcarm[0]),-9999);
  m_cnt_ntrk = 0.;
  
  if(!central) return false;  
  int i_trk_Good         = 0;
  int i_trk_Good_Quality = 0;
  float  i_emc_Good         = 0;
  float  i_emc_Good_Quality = 0;
  for(unsigned int itrk=0;itrk<central->get_npart();itrk++)
    {
      PHSnglCentralTrack *d_cnt = central->get_track(itrk);
      
      float zed    = d_cnt->get_zed();
      float mom    = d_cnt->get_mom();
      //short dcarm  = d_cnt->get_dcarm();
      //float the0   = d_cnt->get_the0();
      //float charge = d_cnt->get_charge();
      //float eta0   = -log(tan(the0/2.0)); 
      //float phi0   = d_cnt->get_phi0();      
      //phi0         = atan2(sin(phi0), cos(phi0)); 	  
      //float pt     = mom*sin(the0);
      
      float pc3dphi  = d_cnt->get_pc3dphi();
      float pc3dz    = d_cnt->get_pc3dz();
      
      float ecorr    = d_cnt->get_ecore();
      
      bool Trk_good         = mom<15.0 && mom>0.2 && fabs(pc3dphi)<10.0 && fabs(pc3dz)<10.0 && fabs(zed)<70; 
      bool Trk_good_Quality = d_cnt->get_quality() == 31 || d_cnt->get_quality() == 63                     ;
      if (Trk_good        )	  { i_trk_Good++;           if(ecorr>0) i_emc_Good         += ecorr; }
      if (Trk_good_Quality)	  { i_trk_Good_Quality++;   if(ecorr>0) i_emc_Good_Quality += ecorr; }
    }
  
  m_cnt_ntrk_good   = i_trk_Good        ;
  m_cnt_ntrk_good_Q = i_trk_Good_Quality;
  m_cnt_nemc_good   = i_emc_Good        ; 
  m_cnt_nemc_good_Q = i_emc_Good_Quality; 
  
  return true;
}                // end of PHCentralTrack



void FvtxHigh::do_Buff_Cnt_Fvt_Trk_Reset()
{
  m_NFvtx_Trk_S     = -9999;
  m_NFvtx_Trk_N     = -9999;
  m_cnt_ntrk_good   = -9999;
  m_cnt_ntrk_good_Q = -9999;
  m_cnt_nemc_good   = -9999;
  m_cnt_nemc_good_Q = -9999;
}

//void FvtxHigh::do_Buff_Cnt_Fvt_Trk_Reset(const int ic_ana, const int iz_ana)
//{
//  
//  m_Eve  [ic_ana][iz_ana] = 0;
//  for(int ieve=0;ieve<NEVE;ieve++)
//    {
//      for(int ikind=0;ikind<3;ikind++)
//	m_Buf_SMDRP[ic_ana][iz_ana][ieve][ikind]       = -9999.;
//      
//      for(int ihar=0;ihar<4;ihar++)
//	{
//	  m_Buf_CNTRP[ic_ana][iz_ana][ieve][ihar] = -9999.;
//	  for(int ikind_fvt=0;ikind_fvt<6;ikind_fvt++)
//	    m_Buf_FVTRP[ic_ana][iz_ana][ieve][ikind_fvt][ihar] = -9999.;
//	}
//      
//      m_Fvt_Ntrk[ic_ana][iz_ana][ieve] = 0;
//      m_Cnt_Ntrk[ic_ana][iz_ana][ieve] = 0;
//      for(int itrk=0;itrk<1000;itrk++)
//	{
//	  m_Fvt_Phi     [ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  m_Fvt_Eta     [ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  m_Fvt_DcaX    [ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  m_Fvt_DcaY    [ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  m_Fvt_DcaX_pre[ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  m_Fvt_DcaY_pre[ic_ana][iz_ana][ieve][itrk] = -9999.;
//	  
//	  m_Cnt_Phi     [ic_ana][iz_ana][ieve][itrk] = -9999;
//	  m_Cnt_Pt      [ic_ana][iz_ana][ieve][itrk] = -9999;
//	  m_Cnt_Pc3sdZ  [ic_ana][iz_ana][ieve][itrk] = -9999;
//	  m_Cnt_Pc3sdPhi[ic_ana][iz_ana][ieve][itrk] = -9999;
//	}
//
//      for(int ipmt=0;ipmt<128;ipmt++)
//	{
//	  m_Buf_BbcPmt_Phi[ic_ana][iz_ana][ieve][ipmt] = -9999;
//	  m_Buf_BbcPmt_Eta[ic_ana][iz_ana][ieve][ipmt] = -9999;
//	  m_Buf_BbcPmt_Cha[ic_ana][iz_ana][ieve][ipmt] = -9999;
//	}
//
//
//      m_Buf_MpcNTower[ic_ana][iz_ana][ieve] = 0;
//      for(int itower=0;itower<400;itower++)
//	{
//	  m_Buf_Mpc_Et [ic_ana][iz_ana][ieve][itower] = -9999;
//	  m_Buf_Mpc_Phi[ic_ana][iz_ana][ieve][itower] = -9999;
//	  m_Buf_Mpc_Eta[ic_ana][iz_ana][ieve][itower] = -9999;	  
//	}
//
//
//    }  
//  
//}


//==============================================================
void FvtxHigh::initHisto()
{
  cout << " Start Init Hist " << endl;
  
  char name[200];
  //Hist Event check
  //unsigned int nzbin = rpcalibv1->GetDetNz(RP::ID_BBC);
  unsigned int nzbin = 12;
  
  sprintf(name, "centzv");
  HistoManager->registerHisto(new TH2F(name, name, 100, 0.,800., nzbin, -30.,30.) );
  m_bbcq = static_cast<TH2*>(HistoManager->getHisto(name));
  //m_bbcq   = new TH2F(name, name, 100, 0.,800., nzbin, -30.,30.);
  
  sprintf(name, "centzv");
  HistoManager->registerHisto(new TH2F(name, name, RP::NMUL3, 0.,100., nzbin, -30.,30.) );
  m_centzv = static_cast<TH2*>(HistoManager->getHisto(name));
  //m_centzv = new TH2F(name, name, RP::NMUL3, 0.,100., nzbin, -30.,30.);  

  //nzbin = rpcalibv1->GetDetNz(RP::ID_SVX);
  nzbin = 10;
  sprintf(name, "centzv_svx");
  //m_centzv_svx  = new TH2F(name, name, RP::NMUL3, 0.,100., nzbin, -10.,10.); 
  HistoManager->registerHisto(new TH2F(name, name, RP::NMUL3, 0.,100., nzbin, -10.,10.) ); 
  m_centzv_svx = static_cast<TH2*>(HistoManager->getHisto(name));
  
  sprintf(name,"NFvtx_phi_vs_eta");  
  HistoManager->registerHisto(new TH2F(name,name,100,-M_PI,M_PI,100,-5,5));
  h_Fvt_phi_eta = static_cast<TH2*>(HistoManager->getHisto(name));

  // Bbc vs Fvtx 
  for(int itrig=0;itrig<3;itrig++)
    {
      // S vs N ( Bbc, Fvtx )
      sprintf(name, "NBbc_S_vs_N_Tr%d", itrig);
      HistoManager->registerHisto(new TH2F(name, name, 100,0,200, 100,0,200));
      h_NBbc_S_vs_N       [itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
      
      sprintf(name, "NFvtx_S_vs_N_Tr%d", itrig);
      HistoManager->registerHisto(new TH2F(name, name, 100,0,200, 100,0,200));
      h_NFvt_S_vs_N       [itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
      
      for(int isn=0;isn<3;isn++)
	{
	  // Bbc vs Fvtx
	  sprintf(name, "NBbc_vs_NFvt_S%d_Tr%d", isn, itrig);
	  HistoManager->registerHisto(new TH2F(name, name, 100,0,200, 100,0,200));
	  h_NBbc_vs_NFvt   [isn][itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
	  
	  // Cnt vs Fvtx 
	  sprintf(name, "NCnt_vs_NFvt_S%d_Tr%d", isn, itrig);
	  HistoManager->registerHisto(new TH2F(name, name, 100,0,100, 100,0,200));
	  h_NCntTrk_vs_NFvt[isn][itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
	  
	  // Cnt Q vs Fvtx
	  sprintf(name, "NCntQ_vs_NFvt_S%d_Tr%d", isn, itrig);
	  HistoManager->registerHisto(new TH2F(name, name, 100,0,100, 100,0,200));
	  h_NCntTrkQ_vs_NFvt[isn][itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
	  
	  // CntEmc vs Fvtx
	  sprintf(name, "NCntEmc_vs_NFvt_S%d_Tr%d", isn, itrig);
	  HistoManager->registerHisto(new TH2F(name, name, 100,0,100, 100,0,200));
	  h_NCntEmc_vs_NFvt [isn][itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
	  
	  // CntEmcQ vs Fvtx
	  sprintf(name, "NCntEmcQ_vs_NFvt_S%d_Tr%d", isn, itrig);
	  HistoManager->registerHisto(new TH2F(name, name, 100,0,100, 100,0,200));
	  h_NCntEmcQ_vs_NFvt[isn][itrig] = static_cast<TH2*>(HistoManager->getHisto(name));
	}//isn
      cout << " aa " << itrig << endl;
    }//itrig

  cout << " End Init Hist " << endl;
  
  return;
}


bool FvtxHigh::Get_Enable_Det_Self(const int detid)
{
  bool status = false;
  //if(detid==RP::ID_SVX) return true;
  //if(detid==RP::ID_SEG) return true;
  if(detid==RP::ID_MPC) return true;
  if(detid==RP::ID_BBC) return true;
  //if(detid==RP::ID_SMD) return true;
  if(detid==RP::ID_CNT) return true;
  //if(detid==RP::ID_FVT) return true;
  
  return status;
}


bool FvtxHigh::Get_Enable_Det_Self_vnRP(const int detid)
{
  bool status = false;
  //if(detid==RP::ID_SVX) return true;
  //if(detid==RP::ID_SEG) return true;
  if(detid==RP::ID_MPC) return true;
  if(detid==RP::ID_BBC) return true;
  //if(detid==RP::ID_SMD) return true;
  if(detid==RP::ID_CNT) return true;
  //if(detid==RP::ID_FVT) return true;
  
  return status;
}
*/
