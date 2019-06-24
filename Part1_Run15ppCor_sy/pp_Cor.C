#include "pp_Cor.h"

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
#include "mpcRawContainer.h"
#include "mpcRawContent.h"

#include <TRxnpRawScintMap.h>
#include <TRxnpScintMap.h>
#include <TRxnpRawXangMap.h>

#include "TrigLvl1.h"
//Reaction Plane Object 

#include "VtxOut.h"

//BBC PMT Info.
#include <Bbc.hh>
#include <BbcRaw.h>
#include <BbcCalib.hh>
#include <BbcGeo.hh>
#include <RunToTime.hh>

#include <TFvtxCompactTrkMap.h>

#include "lpcRawv2.h"
#include "lpcRaw.h"
#include "lpcRawHitv2.h"
#include "lpcRawHit.h"

//ROOT Hists
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TNtuple.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

#include "Run15pppc3dphidzcalibsmoothpass1.h"
#include <PHPoint.h>
#include <TFvtxGlobalParCntrl.h>


using namespace std;
//==============================================================
pp_Cor::pp_Cor(const char* output) : 
  SubsysReco("RPCALIBRATOR")
  ,HistoManager(new Fun4AllHistoManager("RPCalibHisto"))
  ,bbccalib(NULL)  
  ,OutputFileName(output)
  ,m_RunNumber(0)
  ,m_EventNumber(0)
{
	//Initialization
	m_centzv[0] = NULL;
	m_centzv[1] = NULL;
	m_centzv[2] = NULL;

	fvtx_mult[0] = NULL;
	fvtx_mult[1] = NULL;
	fvtx_mult[2] = NULL;
	  
	for(int ic_ana=0;ic_ana<NCENTANA;ic_ana++){ // centrality //for pp 1/5/10/20/40/100 
		for(int iz_ana=0;iz_ana<NZANA;iz_ana++){ //z-vtx bin, 10 bins from -10 to 10 cm
			do_Buff_Cnt_Fvt_Trk_Reset(ic_ana,iz_ana);
		}
	}

}
//==============================================================
pp_Cor::~pp_Cor()
{
  delete bbccalib;
  delete bbcgeo;
  delete HistoManager;
}
//==============================================================
int pp_Cor::Init(PHCompositeNode *topNode)
{
  cout << "pp_Cor::Init started..." << endl;
  
  bbccalib = new BbcCalib();
  bbcgeo   = new BbcGeo();
  
  cout << "pp_Cor::Init ended." << endl;
  return 0;
}
//==============================================================
int pp_Cor::InitRun(PHCompositeNode *topNode)
{
  cout << "pp_Cor::InitRun started..." << endl;
  
  //Get Run number 
  recoConsts *rc = recoConsts::instance();
  m_RunNumber = rc->get_IntFlag("RUNNUMBER");  
  
  ///Calibration BBC PMT
  int icalibversion = 4002;		// Calibration ID(calibrated PMTs are used)
  RunToTime* runTime = RunToTime::instance();
  PHTimeStamp* ts( runTime->getBeginTime(m_RunNumber) );
  PHTimeStamp tstart = *ts;
  delete ts;
  int bbccalib_version = 4002;
	cout << "pp_Cor::InitRun, run number= " << m_RunNumber
		   << " Version of BbcCalib = " << bbccalib_version << endl;
  bbccalib->restore(tstart, bbccalib_version);
  rc->set_IntFlag("BBCCALIBVERSION", icalibversion);  
  
  cout << "pp_Cor::InitRun   run number      = " << m_RunNumber << endl;
  
	//Initialization of histograms
	initHisto();
  
  cout << "pp_Cor::InitRun ended." << endl;
  
  return 0;    
}
//==============================================================
int pp_Cor::process_event(PHCompositeNode *topNode)
{
	//Check event number
	if (m_EventNumber%10000==0) cout << m_EventNumber<< "//------------------------------------" << endl;
	m_EventNumber++;

	//Get nodes and Check
	PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal"); //Global info
	PreviousEvent *peve = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent"); //Tick
	PHCentralTrack *central = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack"); //Central arm track
	BbcRaw *bbcraw  = findNode::getClass<BbcRaw>(topNode,"BbcRaw"); //BBC info
	TFvtxCompactTrkMap* trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode,"TFvtxCompactTrkMap"); //FVTX track
	TrigLvl1 *triglvl1 = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1"); //Trigger info
	VtxOut *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut"); //Trigger info

	if(!global      ) { cerr << PHWHERE << " No PHGloval object !"            << endl; return EVENT_OK; }
	if(!peve        ) { cerr << PHWHERE << " No PreviousEvent object !"       << endl; return EVENT_OK; }
	if(!bbcraw      ) { cerr << PHWHERE << " Could not find BbcRaw !"         << endl; return EVENT_OK; }
	if(!central     ) { cerr << PHWHERE << " CNT Event skip "                 << endl; return EVENT_OK; }    
	//if(!bbcraw      ) {  }
	//if(!central     ) {  }    
	if(!trkfvtx_map ) { cerr << PHWHERE << " No TFvtxCompactTrkMap object !"  << endl; return EVENT_OK; }
	if(!triglvl1    ) { cout << PHWHERE << "Could not find TrigLvl1 !" << endl; return 1; } 
	if(!vtxout          ) { cout << PHWHERE << "Could not find VtxOut !" << endl; return 1; }
  
	//Get scaled trigger bit
	//int Scaled = triglvl1->get_lvl1_trigscaled();

	//Check trigger bit : Trigger bit needed to be confirm in different run
  	//int itrig=-9999; 

	//comment out when FVTX triggered dataset use
	//if ( (Scaled&0x00000001) || (Scaled&0x00000002) || (Scaled&0x00000010) ) itrig = 0; //MB BBCLL1 +-30, novertex, narrowvtx 
	//if ( itrig<0 ) return ABORTEVENT;

	//To see only FVTX and triggered events
	//if ( Scaled & 0x00000400 ) itrig = 1; //FVTX And triggered 
	//if ( itrig<1 ) return ABORTEVENT;
  
	//BBC Charge
	float bbccha_s = -9999;	//BBC charge, South
	float bbccha_n = -9999;	//BBC charge, North
	float bbccha_t = -9999;	//BBC charge, Sum (S+N)

	//BBC charge, How to weight (1GeV=1ptl, 3GeV=3ptls)
	bbccha_s = global->getBbcChargeS();
	bbccha_n = global->getBbcChargeN();
	bbccha_t = bbccha_s + bbccha_n;   

	//BBC Centrality binning for pp
	int ic_ana = -9999;

	//Get FVTX trk vs. BBC charge histogram
        m_FVTXtrk_BBCcha[0]->Fill(bbccha_s, fvtx_mult[0]);
        m_FVTXtrk_BBCcha[1]->Fill(bbccha_n, fvtx_mult[1]);
        m_FVTXtrk_BBCcha[2]->Fill(bbccha_t, fvtx_mult[2]);
	
	//BBC Centrality binning for pp MB
	//if	( bbccha_t>26.5 ) ic_ana = 0;	//0-1%
	//else if ( bbccha_t>19.5 ) ic_ana = 1;	//1-5%
	//else if ( bbccha_t>13.5 ) ic_ana = 2;	//5-20%
	//else			  ic_ana = 3;	//20-100%

	//BBC Centrality binning for pp FVTX
	if	( bbccha_t>40.5 ) ic_ana = 0;	//0-1%
	else if ( bbccha_t>30.5 ) ic_ana = 1;	//1-5%
	else if ( bbccha_t>20.5 ) ic_ana = 2;	//5-20%
	else			  ic_ana = 3;	//20-100%

	//FVTX Track Centrality
	//ic_ana = GetFVTXCentralityBin(trkfvtx_map);
	GetFVTXCentralityBin(trkfvtx_map);

	if ( ic_ana<0 ) return EVENT_OK;

	//This could be redundant due to the trigger bit check above
	bool MinBias = bbccha_s>0 && bbccha_n>0;
	if(!MinBias) return EVENT_OK;
   
	//Get BBC-Z vertex
	float vertex = global->getBbcZVertex();

        //Get Fvtx vertex info
        PHPoint fvtx_vtx = vtxout->get_Vertex("FVTX");
        fvtx_x = fvtx_vtx.getX();
        fvtx_y = fvtx_vtx.getY();
        fvtx_z = fvtx_vtx.getZ();
        
        //beam average for DCA calculation
	if (TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy")){
		//cout << "using beam center" << endl;
		fvtx_x = TFvtxGlobalParCntrl::get_float_par("beam_x_seed");
		fvtx_y = TFvtxGlobalParCntrl::get_float_par("beam_y_seed");
        }

	if (fabs(vertex)>10.)	return EVENT_OK; //BBC narrow vertex && FVTX trig
	if (isnan(fvtx_z))	return EVENT_OK; //skip nan-fvtxz event

	//Z-vertex binning
	int iz_ana =-9999;
	iz_ana = (int)(vertex+10.0)/2.;    //--10~-8, -8~-6, -6~-4, -4~-2, -2~0, 0-2, 2-4, 4-6, 6-8, 8-10
	if(iz_ana<0 || iz_ana>9) { cout << " Out of range " << iz_ana << " " << vertex << endl; return EVENT_OK;}
   
	//Fill histogram (BBC Charge sum, BBC-Z vertex)
	m_centzv[2]->Fill(bbccha_t, vertex);
	m_centzv[1]->Fill(bbccha_n, vertex);
	m_centzv[0]->Fill(bbccha_s, vertex);

	int ieve   = m_Eve[ic_ana][iz_ana]; //Number of events[Trigger][Z-vertex bin]
	//Get BBC information
	GetBBCPMT (ic_ana, iz_ana, ieve, bbcraw);
	//Get FVTX track information
	GetFVTXTrack (ic_ana, iz_ana, ieve, global, trkfvtx_map);
	//Get central arm track information
	GetCNTPC3_TOF (ic_ana, iz_ana, ieve, central,global);
   
	//Count number of buffered events
	m_Eve[ic_ana][iz_ana]++;
	//If the number of buffered events in [Trigger][Z-vtx bin] is 6, Get correlation between buffered events
	if(m_Eve[ic_ana][iz_ana]==NEVE) // Real=0, Mix=1,2,3 : NEVE=4
	{
		//2PC method CNT vs BBC
		do_Real_Mix_CNT_BBC_Pmt_Cor(ic_ana,iz_ana);
		do_Real_Mix_BBC_S_N_Pmt_Cor(ic_ana,iz_ana);	      

		//2PC method CNT vs FVT
		do_Real_Mix_CNT_FVT_Trk_Cor(ic_ana,iz_ana);
		do_Real_Mix_FVT_S_N_Trk_Cor(ic_ana,iz_ana);

		//2PC method CNT vs CNT
		do_Real_Mix_CNT_CNT_Cor(ic_ana,iz_ana);

		//Initialize buffer
		do_Buff_Cnt_Fvt_Trk_Reset(ic_ana,iz_ana);
	}     
   
   return 0;
 }

//==============================================================
int pp_Cor::End(PHCompositeNode *topNode)
{
  cout << "pp_Cor::End" << endl;
  
  HistoManager->dumpHistos(OutputFileName);
  return 0; 
}

//==============================================================
bool pp_Cor::GetBBCPMT(int ic_ana, int iz_ana, int ieve, BbcRaw* bbcraw)
{
  if(!bbcraw) return false;

	// beam beam counter (bbc) r.p. // -----------------------------------
	for (int ipmt=0; ipmt<128; ipmt++) {
		short adc = bbcraw->get_Adc(ipmt);
		short tdc = bbcraw->get_Tdc0(ipmt);
		float time0  = bbccalib->getHitTime0(ipmt, tdc, adc);
		float charge = bbccalib->getCharge(ipmt, adc);
		float bbcx = bbcgeo->getX(ipmt);
		float bbcy = bbcgeo->getY(ipmt);
		float bbcz = bbcgeo->getZ(ipmt);
		if (time0>0 && charge>0) {
			//int iarm = 0;          //S
			//if (bbcz > 0) iarm = 1;//N
			float phi = atan2(bbcy,bbcx);
			float val = charge;
			float rad = sqrt(bbcx*bbcx+bbcy*bbcy);
			float the = atan2(rad,bbcz); // bbcz-bbcv*10.0);
			float eta = -log(tan(0.5*the));
			//float fab_eta = fabs(eta);
			//Symmetry
			//if(ipmt==10||ipmt==21||ipmt==42||ipmt==53||ipmt==74||ipmt==85||ipmt==106||ipmt==117) continue;

			//Remove hot pmts
			//if(ipmt==118 || ipmt==72 || ipmt==122 || ipmt==86 || ipmt==108) continue;
                        //Fill track buffer [Trigger][Z-vtx bin][Event Number][PMT number]
			m_Buf_BbcPmt_Phi[ic_ana][iz_ana][ieve][ipmt] = phi;
			m_Buf_BbcPmt_Eta[ic_ana][iz_ana][ieve][ipmt] = eta;
			m_Buf_BbcPmt_Cha[ic_ana][iz_ana][ieve][ipmt] = val;

			//Fill QA histograms (BBC MAP)
			if(eta<0) { h_BbcS_x_vs_y[ic_ana]->Fill(bbcx,bbcy,val);
	                            h_BbcS_ipmt[ic_ana]->Fill(ipmt, val);		}
			else      { h_BbcN_x_vs_y[ic_ana]->Fill(bbcx,bbcy,val);
				    h_BbcN_ipmt[ic_ana]->Fill(ipmt, val);		}
		}//time , charge
	}//ipmt
  
  return true;
}

//==============================================================
bool pp_Cor::GetCNTPC3_TOF(const int ic_ana, const int iz_ana, const int ieve,PHCentralTrack* central, PHGlobal* global)
{
  
	if(!central) return false;  
	int i_trk = 0;
	//float NEcorr = 0;

	//Loop over Central arm tracks
	for(unsigned int itrk=0;itrk<central->get_npart();itrk++)
	{
		PHSnglCentralTrack *d_cnt = central->get_track(itrk);

		float zed    = d_cnt->get_zed();
		float mom    = d_cnt->get_mom();
		float the0   = d_cnt->get_the0();
		float ch     = d_cnt->get_charge();
		float eta0   = -log(tan(the0/2.0)); 
		if(fabs(eta0)>0.35) continue;
		float phi0   = d_cnt->get_phi0();      
		phi0         = atan2(sin(phi0), cos(phi0)); 	  
		float xpc1 = d_cnt->get_ppc1x();
		float ypc1 = d_cnt->get_ppc1y();
		float xpc3 = d_cnt->get_ppc3x();
		float ypc3 = d_cnt->get_ppc3y();

		float phi1 = atan2(ypc1,xpc1);
		float phi3 = atan2(ypc3,xpc3);

		float zpc1 = d_cnt->get_ppc1z();
		float zpc3 = d_cnt->get_ppc3z();
		float pt   = mom*sin(the0);
		float pc3dphi  = d_cnt->get_pc3dphi();
		float pc3dz    = d_cnt->get_pc3dz();
		//float ecorr = d_cnt->get_ecore();
		float sdphi = -9999;
		float sdz   = -9999;
		int   arm   = -100;
		arm	    = d_cnt->get_dcarm();

		//Central arm track quality cut
		bool good_trk_1 = mom<15.0 && mom>0.2  &&  fabs(pc3dphi)<10.0  &&  fabs(pc3dz)<10.0  
		  &&  (d_cnt->get_quality()  == 31 || d_cnt->get_quality()  == 63)  &&  fabs(zed)<70;

		sdphi = calcsdphi(pc3dphi, arm, int(ch), mom);
		sdz   = calcsdz  (pc3dz,   arm, int(ch), mom);
		bool good_trk = fabs(sdphi)<3. && fabs(sdz)<3.;
		bool Pt_max = 0.25<pt && pt<3.;

		//Fill Central arm track buffer
		if (good_trk && Pt_max && good_trk_1)
		{                                               
			if(i_trk<1000 )
			{
				//Fill QA histograms
				h_Cnt_phi_vs_eta->Fill(phi0,eta0);
				h_Cnt_pt        ->Fill(pt       );
				//Fill buffer
				m_Cnt_Phi [ic_ana][iz_ana][ieve][i_trk] = phi0;
				m_Cnt_Pt  [ic_ana][iz_ana][ieve][i_trk] = pt  ;
				m_Cnt_Eta [ic_ana][iz_ana][ieve][i_trk] = eta0;  //Not change
				m_Cnt_Phi1[ic_ana][iz_ana][ieve][i_trk] = phi1;  //pc1
				m_Cnt_Phi3[ic_ana][iz_ana][ieve][i_trk] = phi3;  //pc3
				m_Cnt_Zpc1[ic_ana][iz_ana][ieve][i_trk] = zpc1;
				m_Cnt_Zpc3[ic_ana][iz_ana][ieve][i_trk] = zpc3;

				//pT bin
				int i_pt_bin=0;
				if     ( 0.25< pt && pt<0.35) i_pt_bin = 0;	//30%
				else if( 0.35<=pt && pt<0.45) i_pt_bin = 1;	//50%
				else if( 0.45<=pt && pt<0.65) i_pt_bin = 2;	//70%
				else if( 0.65<=pt && pt<1.00) i_pt_bin = 3;	//90%
				else 			      i_pt_bin = 4;	//100%
				m_Cnt_Pt_Bin[ic_ana][iz_ana][ieve][i_trk] = i_pt_bin;

				/*
				if(ecorr>0)
					NEcorr += ecorr;
				*/

				i_trk++;
			}
		}
	} //CUTS on momentum, pc3sdz && pc3sdphi < 10, quality, zed, n0

	//Number of Central arm tracks in the buffer
	m_Cnt_Ntrk[ic_ana][iz_ana][ieve] = i_trk;  

return true;
}

//==============================================================
void pp_Cor::do_Real_Mix_CNT_BBC_Pmt_Cor(const int ic_ana, const int iz_ana)
{

	//Loop over buffered events
	for(int i_eve=0;i_eve<NEVE;i_eve++)
	{
		//CNT
		int Ntrk_Cnt = m_Cnt_Ntrk[ic_ana][iz_ana][i_eve];
		for(int i_Cnt_trk=0;i_Cnt_trk<Ntrk_Cnt;i_Cnt_trk++)
		{
			float cnt_Phi  = m_Cnt_Phi   [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			float cnt_Pt   = m_Cnt_Pt    [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			float cnt_Eta  = m_Cnt_Eta   [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			int   ipt_bin  = m_Cnt_Pt_Bin[ic_ana][iz_ana][i_eve][i_Cnt_trk];

			for(int j_eve=0;j_eve<NEVE;j_eve++)
			{
				//BBC
				for(int ipmt=0;ipmt<128;ipmt++)
				{
					float bbc_Phi  = m_Buf_BbcPmt_Phi[ic_ana][iz_ana][j_eve][ipmt];
					float bbc_Eta  = m_Buf_BbcPmt_Eta[ic_ana][iz_ana][j_eve][ipmt];
					float bbc_Cha  = m_Buf_BbcPmt_Cha[ic_ana][iz_ana][j_eve][ipmt];
					float bbc_Cha2 = bbc_Cha*bbc_Cha;

					float dPhi = cnt_Phi - bbc_Phi;
					dPhi = atan2(sin(dPhi), cos(dPhi));
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;

					float dEta = cnt_Eta - bbc_Eta;//-0.35 ~ 0.35 - ( -3.8 ~ 3.8 )
					if(cnt_Phi<-9000) continue;
					if(bbc_Phi<-9000) continue;

					if(i_eve==j_eve)
					{
						if(bbc_Eta<0)     h_Cnt_Bbc_dPhi_Real[ic_ana][   0   ]->Fill(cnt_Pt, dPhi, bbc_Cha );
						else              h_Cnt_Bbc_dPhi_Real[ic_ana][   1   ]->Fill(cnt_Pt, dPhi, bbc_Cha );
						h_Cnt_Bbc_dPhi_dEta_Real             [ic_ana][ipt_bin]->Fill(dPhi  , dEta, bbc_Cha );
						h_Cnt_Bbc_dPhi_dEta_Real_sqr	     [ic_ana][ipt_bin]->Fill(dPhi  , dEta, bbc_Cha2);
					}
					else
					{
						if(bbc_Eta<0)     h_Cnt_Bbc_dPhi_Mix [ic_ana][   0   ]->Fill(cnt_Pt, dPhi, bbc_Cha );
						else              h_Cnt_Bbc_dPhi_Mix [ic_ana][   1   ]->Fill(cnt_Pt, dPhi, bbc_Cha );
						h_Cnt_Bbc_dPhi_dEta_Mix              [ic_ana][ipt_bin]->Fill(dPhi  , dEta, bbc_Cha );
						h_Cnt_Bbc_dPhi_dEta_Mix_sqr	     [ic_ana][ipt_bin]->Fill(dPhi  , dEta, bbc_Cha2);
					}
				}//ipmt
			}//cnt trk
		}//j_eve
	}//i_eve
}

//==============================================================
void pp_Cor::do_Real_Mix_BBC_S_N_Pmt_Cor(const int ic_ana, const int iz_ana)
{
  
	//Loop over buffered events
	for(int i_eve=0;i_eve<NEVE;i_eve++)
	{
		//BBC-S
		for(int i_Pmt_S=0;i_Pmt_S<128;i_Pmt_S++)
		{
			float bbc_Phi_S = m_Buf_BbcPmt_Phi[ic_ana][iz_ana][i_eve][i_Pmt_S];
			float bbc_Eta_S = m_Buf_BbcPmt_Eta[ic_ana][iz_ana][i_eve][i_Pmt_S];
			float bbc_Cha_S = m_Buf_BbcPmt_Cha[ic_ana][iz_ana][i_eve][i_Pmt_S];
			float bbc_Cha_S_sqr = bbc_Cha_S*bbc_Cha_S;

			if(bbc_Eta_S>0) continue;
			for(int j_eve=0;j_eve<NEVE;j_eve++)
			{
				//BBC-N
				for(int i_Pmt_N=0;i_Pmt_N<128;i_Pmt_N++)
				{
					float bbc_Phi_N = m_Buf_BbcPmt_Phi[ic_ana][iz_ana][j_eve][i_Pmt_N];
					float bbc_Eta_N = m_Buf_BbcPmt_Eta[ic_ana][iz_ana][j_eve][i_Pmt_N];
					float bbc_Cha_N = m_Buf_BbcPmt_Cha[ic_ana][iz_ana][j_eve][i_Pmt_N];
					float bbc_Cha_N_sqr = bbc_Cha_N*bbc_Cha_N;

					float bbc_Cha_S_N_sqr = bbc_Cha_S_sqr+bbc_Cha_N_sqr;

					if(bbc_Eta_N<0) continue;		  

					if(bbc_Phi_S<-9000) continue;
					if(bbc_Phi_N<-9000) continue;

					float dPhi = bbc_Phi_S - bbc_Phi_N;// -2pi(-pi - (pi)) ~ 2pi(pi - (-pi)) 
					float dEta = bbc_Eta_S - bbc_Eta_N;//-3.8~-3.0  ~ 3.0 ~ 3.8 -> -7.6 ~ - -6.0
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;

					float bbc_Cha_S_N = bbc_Cha_S * bbc_Cha_N;
					if(i_eve==j_eve) 
					{
						h_Bbc_S_N_dPhi_Real         [ic_ana]->Fill(dPhi,        bbc_Cha_S_N);		   
						h_Bbc_S_N_dPhi_dEta_Real    [ic_ana]->Fill(dPhi, dEta,  bbc_Cha_S_N);
						h_Bbc_S_N_dPhi_dEta_Real_sqr[ic_ana]->Fill(dPhi, dEta,  bbc_Cha_S_N_sqr);
						h_Bbc_S_N_correlation_Real  [ic_ana]->Fill(bbc_Phi_S,   bbc_Phi_N, bbc_Cha_S_N);
					}
					else
					{
						h_Bbc_S_N_dPhi_Mix         [ic_ana]->Fill(dPhi,       bbc_Cha_S_N);		    
						h_Bbc_S_N_dPhi_dEta_Mix    [ic_ana]->Fill(dPhi, dEta, bbc_Cha_S_N);
						h_Bbc_S_N_dPhi_dEta_Mix_sqr[ic_ana]->Fill(dPhi, dEta, bbc_Cha_S_N_sqr);
						h_Bbc_S_N_correlation_Mix  [ic_ana]->Fill(bbc_Phi_S,  bbc_Phi_N, bbc_Cha_S_N);
					}
					//////////////////////////////////////////////////////////////////////////////////
				}//fvt trk
			}//cnt trk
		}//j_eve
	}//i_eve
  
}

//==============================================================
void pp_Cor::do_Real_Mix_CNT_CNT_Cor(const int ic_ana, const int iz_ana)
{

	//Loop over buffered event
	for(int i_eve=0;i_eve<NEVE;i_eve++)
	{
		int Ntrk_Cnt_T = m_Cnt_Ntrk[ic_ana][iz_ana][i_eve];
		//CNT Track
		for(int i_Cnt_trk_T=0;i_Cnt_trk_T<Ntrk_Cnt_T;i_Cnt_trk_T++)
		{
			float cnt_Phi_T  = m_Cnt_Phi   [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			//float cnt_Pt_T   = m_Cnt_Pt    [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			float cnt_Eta_T  = m_Cnt_Eta   [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			int   ipt_bin_T  = m_Cnt_Pt_Bin[ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			float cnt_Phi1_T = m_Cnt_Phi1  [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			float cnt_Phi3_T = m_Cnt_Phi3  [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			float cnt_Zpc1_T = m_Cnt_Zpc1  [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];
			float cnt_Zpc3_T = m_Cnt_Zpc3  [ic_ana][iz_ana][i_eve][i_Cnt_trk_T];

			for(int j_eve=0;j_eve<NEVE;j_eve++)
			{
				int Ntrk_Cnt_A = m_Cnt_Ntrk[ic_ana][iz_ana][j_eve];
				//CNT trigger
				for(int i_Cnt_trk_A=0;i_Cnt_trk_A<Ntrk_Cnt_A;i_Cnt_trk_A++)
				{
					float cnt_Phi_A  = m_Cnt_Phi   [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					//float cnt_Pt_A   = m_Cnt_Pt    [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					float cnt_Eta_A  = m_Cnt_Eta   [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					int   ipt_bin_A  = m_Cnt_Pt_Bin[ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					float cnt_Phi1_A = m_Cnt_Phi1  [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					float cnt_Phi3_A = m_Cnt_Phi3  [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					float cnt_Zpc1_A = m_Cnt_Zpc1  [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];
					float cnt_Zpc3_A = m_Cnt_Zpc3  [ic_ana][iz_ana][j_eve][i_Cnt_trk_A];

					if(i_eve==j_eve && i_Cnt_trk_T == i_Cnt_trk_A) continue;

					float dPc1 = cnt_Phi1_T - cnt_Phi1_A;
					dPc1 = atan2(sin(dPc1),cos(dPc1));
					dPc1 = fabs(dPc1);

					float dPc3 = cnt_Phi3_T - cnt_Phi3_A;
					dPc3 = atan2(sin(dPc3),cos(dPc3));
					dPc3 = fabs(dPc3);

					float dZpc1 = fabs(cnt_Zpc1_T - cnt_Zpc1_A);
					float dZpc3 = fabs(cnt_Zpc3_T - cnt_Zpc3_A);

					//dPc1 vs dZpc1
					if(i_eve==j_eve) h_dPc1_vs_dZpc1[0][0]->Fill(dPc1, dZpc1); 
					else             h_dPc1_vs_dZpc1[1][0]->Fill(dPc1, dZpc1);

					//dPc3 vs dZpc3
					if(i_eve==j_eve) h_dPc3_vs_dZpc3[0][0]->Fill(dPc3, dZpc3); 
					else             h_dPc3_vs_dZpc3[1][0]->Fill(dPc3, dZpc3);

					float Cut1 = sqrt( pow(fabs(dPc1)/0.040, 2.0)  +  pow(fabs( dZpc1 )/90.0, 2.0) );
					float Cut2 = sqrt( pow(fabs(dPc1)/0.080, 2.0)  +  pow(fabs( dZpc1 )/8.00, 2.0) );
					float Cut3 = sqrt( pow(fabs(dPc3)/0.070, 2.0)  +  pow(fabs( dZpc3 )/25.0, 2.0) );		  
					if( Cut1<1. || Cut2<1. || Cut3<1.) continue;

					//dPc1 vs dZpc1 Cut
					if(i_eve==j_eve) h_dPc1_vs_dZpc1[0][1]->Fill(dPc1, dZpc1); 
					else             h_dPc1_vs_dZpc1[1][1]->Fill(dPc1, dZpc1);
					//dPc3 vs dZpc3 Cut		      
					if(i_eve==j_eve) h_dPc3_vs_dZpc3[0][1]->Fill(dPc3, dZpc3); 
					else             h_dPc3_vs_dZpc3[1][1]->Fill(dPc3, dZpc3);


					float dPhi = cnt_Phi_T - cnt_Phi_A;
					dPhi = atan2(sin(dPhi), cos(dPhi));
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
					float dEta = cnt_Eta_T - cnt_Eta_A;//-0.35 ~ 0.35 - ( -3.8 ~ 3.8 )
					if(cnt_Phi_T<-9000) continue;
					if(cnt_Phi_A<-9000) continue;

					if(i_eve==j_eve)
					{
						h_Cnt_Cnt_dPhi_dEta_Real[ic_ana][ipt_bin_T][ipt_bin_A]->Fill(dPhi, dEta);
					}
					else
					{
						h_Cnt_Cnt_dPhi_dEta_Mix [ic_ana][ipt_bin_T][ipt_bin_A]->Fill(dPhi, dEta);
					}
				}//i_Cnt_trk_A
			}//j_eve
		}//i_Cnt_trk_T
	}//i_eve
}

//==============================================================
void pp_Cor::do_Real_Mix_CNT_FVT_Trk_Cor(const int ic_ana, const int iz_ana)
{

	//Loop over buffered events
	for(int i_eve=0;i_eve<NEVE;i_eve++)
	{
		int Ntrk_Cnt = m_Cnt_Ntrk[ic_ana][iz_ana][i_eve];
		//CNT Track
		for(int i_Cnt_trk=0;i_Cnt_trk<Ntrk_Cnt;i_Cnt_trk++)
		{
			float cnt_Phi  = m_Cnt_Phi   [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			float cnt_Pt   = m_Cnt_Pt    [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			float cnt_Eta  = m_Cnt_Eta   [ic_ana][iz_ana][i_eve][i_Cnt_trk];
			int   ipt_bin  = m_Cnt_Pt_Bin[ic_ana][iz_ana][i_eve][i_Cnt_trk];

			for(int j_eve=0;j_eve<NEVE;j_eve++)
			{
				int Ntrk_Fvt = m_Fvt_Ntrk[ic_ana][iz_ana][j_eve];
				//FVTX Track
				for(int i_Fvt_trk=0;i_Fvt_trk<Ntrk_Fvt;i_Fvt_trk++)
				{
					float fvt_Phi             = m_Fvt_Phi     [ic_ana][iz_ana][j_eve][i_Fvt_trk];
					float fvt_Eta             = m_Fvt_Eta     [ic_ana][iz_ana][j_eve][i_Fvt_trk];
					float fvt_dca_X           = m_Fvt_DcaX    [ic_ana][iz_ana][j_eve][i_Fvt_trk];
					float fvt_dca_Y           = m_Fvt_DcaY    [ic_ana][iz_ana][j_eve][i_Fvt_trk];
					//float fvt_dca_R_xcos_ysin = fvt_dca_X*cos(fvt_Phi) + fvt_dca_Y*sin(fvt_Phi);
					float fvt_dca_R = fvt_dca_X + fvt_dca_Y;

					bool   fvt_dca_R_2  = fabs( fvt_dca_R )<2  ;
					bool   fvt_dca_R_1  = fabs( fvt_dca_R )<1  ;
					bool   fvt_dca_R_05 = fabs( fvt_dca_R )<0.5;

					float dPhi = cnt_Phi - fvt_Phi;
					dPhi = atan2(sin(dPhi),cos(dPhi));
					float dEta = cnt_Eta - fvt_Eta;
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;

					if(cnt_Phi<-9000) continue;
					if(fvt_Phi<-9000) continue;

					if(i_eve==j_eve)
					{
						if(fvt_Eta<0)
						{
							h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][0][0]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][0][1]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][0][2]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][0][3]->Fill(cnt_Pt, dPhi);
						}
						else
						{
							h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][1][0]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][1][1]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][1][2]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][1][3]->Fill(cnt_Pt, dPhi);
						}
						
						//2D hist eta vs phi
						h_Cnt_Fvt_dPhi_dEta_Real[ic_ana][iz_ana][ipt_bin][0]->Fill(dPhi,dEta);
						if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_dEta_Real[ic_ana][iz_ana][ipt_bin][1]->Fill(dPhi,dEta);
						if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_dEta_Real[ic_ana][iz_ana][ipt_bin][2]->Fill(dPhi,dEta);
						if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_dEta_Real[ic_ana][iz_ana][ipt_bin][3]->Fill(dPhi,dEta);
					}
					else
					{
						if(fvt_Eta<0)
						{
							h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][0][0]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][0][1]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][0][2]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][0][3]->Fill(cnt_Pt, dPhi);
						}
						else
						{
							h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][1][0]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][1][1]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][1][2]->Fill(cnt_Pt, dPhi);
							if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][1][3]->Fill(cnt_Pt, dPhi);
						}

						h_Cnt_Fvt_dPhi_dEta_Mix [ic_ana][iz_ana][ipt_bin][0]->Fill(dPhi,dEta);
						if(fvt_dca_R_2 ) h_Cnt_Fvt_dPhi_dEta_Mix [ic_ana][iz_ana][ipt_bin][1]->Fill(dPhi,dEta);
						if(fvt_dca_R_1 ) h_Cnt_Fvt_dPhi_dEta_Mix [ic_ana][iz_ana][ipt_bin][2]->Fill(dPhi,dEta);
						if(fvt_dca_R_05) h_Cnt_Fvt_dPhi_dEta_Mix [ic_ana][iz_ana][ipt_bin][3]->Fill(dPhi,dEta);
					}
					//////////////////////////////////////////////////////////////////////////////////
				}//fvt trk
			}//cnt trk
		}//j_eve
	}//i_eve

}

//==============================================================
void pp_Cor::do_Real_Mix_FVT_S_N_Trk_Cor(const int ic_ana, const int iz_ana)
{

	//Loop over buffered event
	for(int i_eve=0;i_eve<NEVE;i_eve++)
	{
		for(int j_eve=0;j_eve<NEVE;j_eve++)
		{
			int Ntrk_Fvt_S = m_Fvt_Ntrk[ic_ana][iz_ana][i_eve];
			int Ntrk_Fvt_N = m_Fvt_Ntrk[ic_ana][iz_ana][j_eve];

			//FVT-S
			for(int i_Fvt_trk_S=0;i_Fvt_trk_S<Ntrk_Fvt_S;i_Fvt_trk_S++)
			{
				float fvt_Phi_S       = m_Fvt_Phi     [ic_ana][iz_ana][i_eve][i_Fvt_trk_S];
				float fvt_Eta_S       = m_Fvt_Eta     [ic_ana][iz_ana][i_eve][i_Fvt_trk_S];
				float fvt_dca_X_S = m_Fvt_DcaX[ic_ana][iz_ana][i_eve][i_Fvt_trk_S];
				float fvt_dca_Y_S = m_Fvt_DcaY[ic_ana][iz_ana][i_eve][i_Fvt_trk_S];
				//float fvt_dca_R_S_xcos_ysin = fvt_dca_X_S*cos(fvt_Phi_S) + fvt_dca_Y_S*sin(fvt_Phi_S);
				float fvt_dca_R_S = fvt_dca_X_S + fvt_dca_Y_S;

				if(fvt_Eta_S>0) continue;
				//Fvt-N
				for(int i_Fvt_trk_N=0;i_Fvt_trk_N<Ntrk_Fvt_N;i_Fvt_trk_N++)
				{
					float fvt_Phi_N       = m_Fvt_Phi     [ic_ana][iz_ana][j_eve][i_Fvt_trk_N];
					float fvt_Eta_N       = m_Fvt_Eta     [ic_ana][iz_ana][j_eve][i_Fvt_trk_N];
					float fvt_dca_X_N = m_Fvt_DcaX[ic_ana][iz_ana][j_eve][i_Fvt_trk_N];
					float fvt_dca_Y_N = m_Fvt_DcaY[ic_ana][iz_ana][j_eve][i_Fvt_trk_N];
					//float fvt_dca_R_N_xcos_ysin = fvt_dca_X_N*cos(fvt_Phi_N) + fvt_dca_Y_N*sin(fvt_Phi_N);
					float fvt_dca_R_N = fvt_dca_X_N + fvt_dca_Y_N;

					if(fvt_Eta_N<0) continue;		  

					bool fvt_dca_R_S_N_2  = fabs(fvt_dca_R_S)<2   && fabs(fvt_dca_R_N)<2  ;
					bool fvt_dca_R_S_N_1  = fabs(fvt_dca_R_S)<1   && fabs(fvt_dca_R_N)<1  ;
					bool fvt_dca_R_S_N_05 = fabs(fvt_dca_R_S)<0.5 && fabs(fvt_dca_R_N)<0.5;

					float dPhi = fvt_Phi_S - fvt_Phi_N;
					dPhi = atan2(sin(dPhi),cos(dPhi));
					float dEta = fvt_Eta_S - fvt_Eta_N;//-2.5 ~ -1.5  , 1.5 ~ 2.5 -> dEta = -5 ~ -3
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;

					if(fvt_Phi_S<-9000) continue;
					if(fvt_Phi_N<-9000) continue;

					if(i_eve==j_eve)
					{
						h_Fvt_S_N_dPhi_Real[ic_ana][iz_ana][0]->Fill(dPhi);
						if( fvt_dca_R_S_N_2 ) h_Fvt_S_N_dPhi_Real[ic_ana][iz_ana][1]->Fill(dPhi);
						if( fvt_dca_R_S_N_1 ) h_Fvt_S_N_dPhi_Real[ic_ana][iz_ana][2]->Fill(dPhi);
						if( fvt_dca_R_S_N_05) h_Fvt_S_N_dPhi_Real[ic_ana][iz_ana][3]->Fill(dPhi);

						h_Fvt_S_N_dPhi_dEta_Real[ic_ana][iz_ana][0]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_2 ) h_Fvt_S_N_dPhi_dEta_Real[ic_ana][iz_ana][1]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_1 ) h_Fvt_S_N_dPhi_dEta_Real[ic_ana][iz_ana][2]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_05) h_Fvt_S_N_dPhi_dEta_Real[ic_ana][iz_ana][3]->Fill(dPhi,dEta);
					}
					else
					{
						h_Fvt_S_N_dPhi_Mix [ic_ana][iz_ana][0]->Fill(dPhi);
						if( fvt_dca_R_S_N_2 ) h_Fvt_S_N_dPhi_Mix [ic_ana][iz_ana][1]->Fill(dPhi);
						if( fvt_dca_R_S_N_1 ) h_Fvt_S_N_dPhi_Mix [ic_ana][iz_ana][2]->Fill(dPhi);
						if( fvt_dca_R_S_N_05) h_Fvt_S_N_dPhi_Mix [ic_ana][iz_ana][3]->Fill(dPhi);

						h_Fvt_S_N_dPhi_dEta_Mix[ic_ana][iz_ana][0]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_2 ) h_Fvt_S_N_dPhi_dEta_Mix[ic_ana][iz_ana][1]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_1 ) h_Fvt_S_N_dPhi_dEta_Mix[ic_ana][iz_ana][2]->Fill(dPhi,dEta);
						if( fvt_dca_R_S_N_05) h_Fvt_S_N_dPhi_dEta_Mix[ic_ana][iz_ana][3]->Fill(dPhi,dEta);
					}
					//////////////////////////////////////////////////////////////////////////////////
				}//fvt trk
			}//cnt trk
		}//j_eve
	}//i_eve
  
}

//==============================================================
void pp_Cor::do_Buff_Cnt_Fvt_Trk_Reset(const int ic_ana, const int iz_ana)
{
  
  m_Eve  [ic_ana][iz_ana] = 0;
	for(int ieve=0;ieve<NEVE;ieve++)
	{
		//CNT
		m_Cnt_Ntrk[ic_ana][iz_ana][ieve] = 0;
		for(int itrk=0;itrk<1000;itrk++)
		{
			m_Cnt_Phi   [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Pt    [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Eta   [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Pt_Bin[ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Phi1  [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Phi3  [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Zpc1  [ic_ana][iz_ana][ieve][itrk] = -9999;
			m_Cnt_Zpc3  [ic_ana][iz_ana][ieve][itrk] = -9999;
		}
      
		//BBC
		for(int ipmt=0;ipmt<128;ipmt++)
		{
			m_Buf_BbcPmt_Phi[ic_ana][iz_ana][ieve][ipmt] = -9999;
			m_Buf_BbcPmt_Eta[ic_ana][iz_ana][ieve][ipmt] = -9999;
			m_Buf_BbcPmt_Cha[ic_ana][iz_ana][ieve][ipmt] = -9999;
		}
      
		//FVTX
		m_Fvt_Ntrk[ic_ana][iz_ana][ieve] = 0;
		for(int itrk=0;itrk<1000;itrk++)
		{
			m_Fvt_Phi     [ic_ana][iz_ana][ieve][itrk] = -9999.;
			m_Fvt_Eta     [ic_ana][iz_ana][ieve][itrk] = -9999.;
			m_Fvt_DcaX    [ic_ana][iz_ana][ieve][itrk] = -9999.;
			m_Fvt_DcaY    [ic_ana][iz_ana][ieve][itrk] = -9999.;
			m_Fvt_DcaX_pre[ic_ana][iz_ana][ieve][itrk] = -9999.;
			m_Fvt_DcaY_pre[ic_ana][iz_ana][ieve][itrk] = -9999.;
		}    
	}//ieve
  
}


//==============================================================
void pp_Cor::initHisto()
{
  cout << " Start Init Hist " << endl;

  char name[200];
  //Hist Event check
	for (int iarm=0; iarm<3; iarm++){
		sprintf(name, "QA_hFvtxTrk_arm%d",iarm);
		HistoManager->registerHisto(new TH1F(name,name,200,0,200));
		hFvtxTrk[iarm] = static_cast<TH1*>(HistoManager->getHisto(name));

		sprintf(name, "QA_centzv_arm%d", iarm);
       		HistoManager->registerHisto(new TH2F(name, name, 250, 0.,250., NZANA, -10.,10.) );
       		m_centzv[iarm] = static_cast<TH2*>(HistoManager->getHisto(name));

                sprintf(name, "QA_BBCvsFVTX_arm%d", iarm);
                HistoManager->registerHisto(new TH2F(name, name, 250, 0.,250., 250, 0.,250.) );
                m_FVTXtrk_BBCcha[iarm] = static_cast<TH2*>(HistoManager->getHisto(name));
	}

	sprintf(name, "QA_Cnt_phi_vs_eta");
	HistoManager->registerHisto(new TH2F(name, name, 100, -M_PI, M_PI, 100, -0.5,0.5));
	h_Cnt_phi_vs_eta = static_cast<TH2*>(HistoManager->getHisto(name));

	sprintf(name, "QA_Cnt_pT");
	HistoManager->registerHisto(new TH1F(name, name, 100, 0,10));
	h_Cnt_pt       = static_cast<TH1*>(HistoManager->getHisto(name));

	sprintf(name, "QA_NCnt_vs_NDch");
	HistoManager->registerHisto(new TH2F(name, name, 15, 0, 15, 100, 0, 100));
	h_NCnt_vs_NDch  = static_cast<TH2*>(HistoManager->getHisto(name));


	for(int ic_ana=0; ic_ana<NCENTANA; ic_ana++){
		sprintf(name,"QA_bbcpmt2dmapS_C%d", ic_ana); 
		HistoManager->registerHisto(new TProfile2D(name,name,75,-150,150,75,-150,150,0,1000));
		h_BbcS_x_vs_y[ic_ana] = static_cast<TProfile2D*>(HistoManager->getHisto(name));

		sprintf(name,"QA_bbcpmt2dmapN_C%d", ic_ana); 
		HistoManager->registerHisto(new TProfile2D(name,name,75,-150,150,75,-150,150,0,1000));
		h_BbcN_x_vs_y[ic_ana] = static_cast<TProfile2D*>(HistoManager->getHisto(name));

                sprintf(name,"QA_bbc_pmt1D_S_C%d", ic_ana);
                HistoManager->registerHisto(new TProfile(name,name, 128, 0, 128));
                h_BbcS_ipmt[ic_ana] = static_cast<TProfile*>(HistoManager->getHisto(name));

                sprintf(name,"QA_bbc_pmt1D_N_C%d", ic_ana);
                HistoManager->registerHisto(new TProfile(name,name, 128, 0, 128));
                h_BbcN_ipmt[ic_ana] = static_cast<TProfile*>(HistoManager->getHisto(name));
	}
	
	for (int iz_ana=0; iz_ana<NZANA; iz_ana++){
	sprintf(name,"QA_FVT_phi_eta_Z%d_Sys0",iz_ana);
        HistoManager->registerHisto(new TH2D(name,name, 100, -M_PI, M_PI, 100, -5, 5));
        h_Fvt_phi_eta[iz_ana][0] = static_cast<TH2D*>(HistoManager->getHisto(name));	

        sprintf(name,"QA_FVT_phi_eta_Z%d_Sys1",iz_ana);
        HistoManager->registerHisto(new TH2D(name,name, 100, -M_PI, M_PI, 100, -5, 5));
        h_Fvt_phi_eta[iz_ana][1] = static_cast<TH2D*>(HistoManager->getHisto(name));
	}

	for (int iarm=0; iarm<2; iarm++){
	for (int ixy=0 ; ixy<2 ; ixy++ ){ 
	for (int im_ana=0; im_ana<3; im_ana++){
	sprintf(name,"QA_FVTX_DCA_ARM%d_XY%d_M%d",iarm, ixy, im_ana);
        HistoManager->registerHisto(new TH1D(name,name, 300, -15, 15));
        h_Fvt_Dca[iarm][ixy][im_ana] = static_cast<TH1D*>(HistoManager->getHisto(name));
	}}}


  // BBC CNT 2PC
	for(int ic_ana=0;ic_ana<NCENTANA; ic_ana++)
	{
		// BBC vs CNT dPhi vs pt
		// MPC vs CNT dPhi vs pt
		for(int ieta=0;ieta<2;ieta++)
		{
			if(ieta==0) sprintf(name,"Real_dPhi_CNT_BBC_C%d_South", ic_ana );
			if(ieta==1) sprintf(name,"Real_dPhi_CNT_BBC_C%d_North", ic_ana );	      
			HistoManager->registerHisto(new TH2D(name, name, 5, 0, 5, 15,-0.5*M_PI,1.5*M_PI));
			h_Cnt_Bbc_dPhi_Real[ic_ana][ieta] = static_cast<TH2D*>(HistoManager->getHisto(name));

			if(ieta==0) sprintf(name,"Mix_dPhi_CNT_BBC_C%d_South", ic_ana);
			if(ieta==1) sprintf(name,"Mix_dPhi_CNT_BBC_C%d_North", ic_ana);
			HistoManager->registerHisto(new TH2D(name, name, 5, 0, 5, 15,-0.5*M_PI,1.5*M_PI));
			h_Cnt_Bbc_dPhi_Mix [ic_ana][ieta] = static_cast<TH2D*>(HistoManager->getHisto(name));

			sprintf(name,"QA_BBC_FVTX_DCA_C%d_ARM%d", ic_ana, ieta);
                        HistoManager->registerHisto(new TH2D(name, name, 100, -5, 5, 100, -5, 5));
                        h_BbV_Fvt_Dca	   [ic_ana][ieta] = static_cast<TH2D*>(HistoManager->getHisto(name));

			sprintf(name,"QA_BBC_FVTX_Z_C%d_ARM%d", ic_ana, ieta);
                        HistoManager->registerHisto(new TH2D(name, name, 600, -30, 30, 600, -30, 30));
                        h_BbV_Fvt_z        [ic_ana][ieta] = static_cast<TH2D*>(HistoManager->getHisto(name));
		}

		for(int ipt_bin=0;ipt_bin<NPTBIN;ipt_bin++)
		{
			//CNT vs BBC Real and Mix
			sprintf(name,"Real_dPhi_dEta_CNT_BBC_C%d_Pt%d", ic_ana, ipt_bin );
			HistoManager->registerHisto(new TH2D(name, name, 18,-0.5*M_PI,1.5*M_PI, 200,-5,5));
			h_Cnt_Bbc_dPhi_dEta_Real[ic_ana][ipt_bin] = static_cast<TH2D*>(HistoManager->getHisto(name));
			
			sprintf(name,"Real_dPhi_dEta_CNT_BBC_sqr_C%d_Pt%d", ic_ana, ipt_bin );	// To calculate error bar
                        HistoManager->registerHisto(new TH2D(name, name, 18,-0.5*M_PI,1.5*M_PI, 200,-5,5));
                        h_Cnt_Bbc_dPhi_dEta_Real_sqr[ic_ana][ipt_bin] = static_cast<TH2D*>(HistoManager->getHisto(name));

			sprintf(name,"Mix_dPhi_dEta_CNT_BBC_C%d_Pt%d", ic_ana, ipt_bin );
			HistoManager->registerHisto(new TH2D(name, name, 18,-0.5*M_PI,1.5*M_PI, 200,-5,5));
			h_Cnt_Bbc_dPhi_dEta_Mix [ic_ana][ipt_bin] = static_cast<TH2D*>(HistoManager->getHisto(name));

                        sprintf(name,"Mix_dPhi_dEta_CNT_BBC_sqr_C%d_Pt%d", ic_ana, ipt_bin );
                        HistoManager->registerHisto(new TH2D(name, name, 18,-0.5*M_PI,1.5*M_PI, 200,-5,5));  // To calculate error bar
                        h_Cnt_Bbc_dPhi_dEta_Mix_sqr[ic_ana][ipt_bin] = static_cast<TH2D*>(HistoManager->getHisto(name));
		}
      
		// dPhi hist
		// BBC S vs N dPhi
		sprintf(name,"Real_dPhi_BBC_S_N_C%d", ic_ana);
		HistoManager->registerHisto(new TH1D(name, name, 12,-0.5*M_PI,1.5*M_PI));
		h_Bbc_S_N_dPhi_Real[ic_ana] = static_cast<TH1D*>(HistoManager->getHisto(name));

		sprintf(name,"Mix_dPhi_BBC_S_N_C%d", ic_ana);
		HistoManager->registerHisto(new TH1D(name, name, 12,-0.5*M_PI,1.5*M_PI));
		h_Bbc_S_N_dPhi_Mix [ic_ana] = static_cast<TH1D*>(HistoManager->getHisto(name));

		/// dPhi vs dEta BBC S vs N
		sprintf(name,"Real_dPhi_dEta_BBC_S_N_C%d", ic_ana );
		HistoManager->registerHisto(new TH2D(name, name, 12,-0.5*M_PI,1.5*M_PI, 60,-8,-5));
		h_Bbc_S_N_dPhi_dEta_Real[ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));

                sprintf(name,"Real_dPhi_dEta_BBC_S_N_sqr_C%d", ic_ana );
                HistoManager->registerHisto(new TH2D(name, name, 12,-0.5*M_PI,1.5*M_PI, 60,-8,-5));
                h_Bbc_S_N_dPhi_dEta_Real_sqr[ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));

		sprintf(name,"Mix_dPhi_dEta_BBC_S_N_C%d", ic_ana );
		HistoManager->registerHisto(new TH2D(name, name, 12,-0.5*M_PI,1.5*M_PI, 60,-8,-5));
		h_Bbc_S_N_dPhi_dEta_Mix [ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));

                sprintf(name,"Mix_dPhi_dEta_BBC_S_N_sqr_C%d", ic_ana );
                HistoManager->registerHisto(new TH2D(name, name, 12,-0.5*M_PI,1.5*M_PI, 60,-8,-5));
                h_Bbc_S_N_dPhi_dEta_Mix_sqr [ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));

                sprintf(name,"QA_Real_Phi_BBC_S_N_Correlation_C%d", ic_ana );
                HistoManager->registerHisto(new TH2D(name, name, 12,-1*M_PI,1*M_PI, 20,-1*M_PI,1*M_PI));
		h_Bbc_S_N_correlation_Real  [ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));

		sprintf(name,"QA_Mix_Phi_BBC_S_N_Correlation_C%d", ic_ana );
                HistoManager->registerHisto(new TH2D(name, name, 12,-1*M_PI,1*M_PI, 20,-1*M_PI,1*M_PI));
                h_Bbc_S_N_correlation_Mix   [ic_ana] = static_cast<TH2D*>(HistoManager->getHisto(name));
	}
  
	for(int ic_ana=0;ic_ana<NCENTANA;ic_ana++)
	{
		// Cnt dPhi vs dEta
		for(int ipt_bin_T=0;ipt_bin_T<NPTBIN;ipt_bin_T++)
		{
			for(int ipt_bin_A=0;ipt_bin_A<NPTBIN;ipt_bin_A++)
			{
				sprintf(name,"Real_dPhi_dEta_CNT_CNT_C%d_TpT%d_ApT%d", ic_ana, ipt_bin_T, ipt_bin_A);
				HistoManager->registerHisto(new TH2D(name, name, 20,-0.5*M_PI,1.5*M_PI, 20,-0.75,0.75));
				h_Cnt_Cnt_dPhi_dEta_Real[ic_ana][ipt_bin_T][ipt_bin_A] = static_cast<TH2D*>(HistoManager->getHisto(name));

				sprintf(name,"Mix_dPhi_dEta_CNT_CNT_C%d_TpT%d_ApT%d", ic_ana, ipt_bin_T, ipt_bin_A);
				HistoManager->registerHisto(new TH2D(name, name, 20,-0.5*M_PI,1.5*M_PI, 20,-0.75,0.75));
				h_Cnt_Cnt_dPhi_dEta_Mix [ic_ana][ipt_bin_T][ipt_bin_A] = static_cast<TH2D*>(HistoManager->getHisto(name));
			}//ipt_bin_A
		}//ipt_bin_A
	}//Itrig
  
	for(int icut=0;icut<2;icut++)
	{
		sprintf(name,"Real_dPc1_dZpc1_%d",icut);
		HistoManager->registerHisto(new TH2D(name, name, 50, 0, 0.5, 100, 0,100));
		h_dPc1_vs_dZpc1[0][icut] = static_cast<TH2D*>(HistoManager->getHisto(name));

		sprintf(name,"Mix_dPc1_dZpc1_%d",icut);
		HistoManager->registerHisto(new TH2D(name, name, 50, 0, 0.5, 100, 0,100));
		h_dPc1_vs_dZpc1[1][icut] = static_cast<TH2D*>(HistoManager->getHisto(name));

		sprintf(name,"Real_dPc3_dZpc3_%d",icut);
		HistoManager->registerHisto(new TH2D(name, name, 50, 0, 0.5, 100, 0,100));
		h_dPc3_vs_dZpc3[0][icut] = static_cast<TH2D*>(HistoManager->getHisto(name));

		sprintf(name,"Mix_dPc3_dZpc3_%d",icut);
		HistoManager->registerHisto(new TH2D(name, name, 50, 0, 0.5, 100, 0,100));
		h_dPc3_vs_dZpc3[1][icut] = static_cast<TH2D*>(HistoManager->getHisto(name));
	}//icut


	for(int ic_ana=0;ic_ana<NCENTANA;ic_ana++)
	{	  
		for(int iz_ana=0;iz_ana<NZANA;iz_ana++)
		{
			for(int isys=0;isys<4;isys++)
			{
				for(int ieta=0;ieta<2;ieta++)
				{
					sprintf(name,"Real_dPhi_CNT_FVT_C%d_Z%d_Eta%d_Sys%d", ic_ana, iz_ana, ieta,isys);
					HistoManager->registerHisto(new TH2D(name, name, 5,0,5, 48,-0.5*M_PI,1.5*M_PI));
					h_Cnt_Fvt_dPhi_Real[ic_ana][iz_ana][ieta][isys] = static_cast<TH2*>(HistoManager->getHisto(name)); 

					sprintf(name,"Mix_dPhi_CNT_FVT_C%d_Z%d_Eta%d_Sys%d", ic_ana, iz_ana, ieta,isys);
					HistoManager->registerHisto(new TH2D(name, name, 5,0,5, 48,-0.5*M_PI,1.5*M_PI));
					h_Cnt_Fvt_dPhi_Mix [ic_ana][iz_ana][ieta][isys] = static_cast<TH2*>(HistoManager->getHisto(name));
				}

				sprintf(name,"Real_dPhi_FVT_S_N_C%d_Z%d_Sys%d", ic_ana, iz_ana, isys);
				HistoManager->registerHisto(new TH1D(name, name,48,-0.5*M_PI,1.5*M_PI));
				h_Fvt_S_N_dPhi_Real[ic_ana][iz_ana][isys] = static_cast<TH1*>(HistoManager->getHisto(name));

				sprintf(name,"Mix_dPhi_FVT_S_N_C%d_Z%d_Sys%d", ic_ana, iz_ana, isys);
				HistoManager->registerHisto(new TH1D(name, name,48,-0.5*M_PI,1.5*M_PI));
				h_Fvt_S_N_dPhi_Mix [ic_ana][iz_ana][isys] = static_cast<TH1*>(HistoManager->getHisto(name));

				sprintf(name,"Real_dPhi_dEta_FVT_S_N_C%d_Z%d_Sys%d", ic_ana, iz_ana, isys);
				HistoManager->registerHisto(new TH2D(name, name,48,-0.5*M_PI,1.5*M_PI, 100, -6, 0));
				h_Fvt_S_N_dPhi_dEta_Real[ic_ana][iz_ana][isys] = static_cast<TH2*>(HistoManager->getHisto(name));

				sprintf(name,"Mix_dPhi_dEta_FVT_S_N_C%d_Z%d_Sys%d", ic_ana, iz_ana, isys);
				HistoManager->registerHisto(new TH2D(name, name,48,-0.5*M_PI,1.5*M_PI, 100, -6, 0));
				h_Fvt_S_N_dPhi_dEta_Mix [ic_ana][iz_ana][isys] = static_cast<TH2*>(HistoManager->getHisto(name));
				for(int ipt_bin=0;ipt_bin<NPTBIN;ipt_bin++)
				{
					sprintf(name,"Real_dPhi_dEta_CNT_FVT_C%d_Z%d_Sys%d_Pt%d", ic_ana, iz_ana, isys,ipt_bin);
					HistoManager->registerHisto(new TH2D(name, name,48,-0.5*M_PI,1.5*M_PI, 100, -3,3));
					h_Cnt_Fvt_dPhi_dEta_Real[ic_ana][iz_ana][ipt_bin][isys]= static_cast<TH2*>(HistoManager->getHisto(name));

					sprintf(name,"Mix_dPhi_dEta_CNT_FVT_C%d_Z%d_Sys%d_Pt%d", ic_ana, iz_ana, isys,ipt_bin);
					HistoManager->registerHisto(new TH2D(name, name,48,-0.5*M_PI,1.5*M_PI, 100, -3,3));
					h_Cnt_Fvt_dPhi_dEta_Mix[ic_ana][iz_ana][ipt_bin][isys]= static_cast<TH2*>(HistoManager->getHisto(name));
				}// ipt_bin
			}//isys
		}//iz_ana
	}//ic_ana
 
  return;
}


bool pp_Cor::GetFVTXTrack(const int ic_ana, const int iz_ana, const int ieve, PHGlobal* global, TFvtxCompactTrkMap *trkfvtx_map)
{
  
	float vertex = global->getBbcZVertex();
	float bbcs = global->getBbcChargeS();
	float bbcn = global->getBbcChargeN();
	//Why doing this again?
	if(bbcs<0) return false;
	if(bbcn<0) return false;

	int ntrk[2] = {0};
  
	int itrk=0;
	//Get FVTXCompactTrkMap
	TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
	//Loop over FVTX tracks in an event
	while (TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next()) {
		TFvtxCompactTrk* fvtx_trk = trk_ptr->get();
		float the = fvtx_trk->get_fvtx_theta();
		float eta = fvtx_trk->get_fvtx_eta();
		float phi = fvtx_trk->get_fvtx_phi();
		float hit = (float)fvtx_trk->get_nhits();
		int   arm = (int)fvtx_trk->get_arm();
		float fvx = fvtx_trk->get_fvtx_vtx().getX();
		float fvy = fvtx_trk->get_fvtx_vtx().getY();
		float fvz = fvtx_trk->get_fvtx_vtx().getZ();
		float chi2 = fvtx_trk->get_chi2_ndf();
		if(the==0 && phi==0 && fvx==0 && fvy==0 && fvz==0) continue;
		if((arm==0 && eta>0) || (arm==1 && eta<0) || arm<0 || arm>1) continue;

		//FVTX track condition
		if( fabs(eta)> 2.5 ) continue;
		if( fabs(eta)< 1.5 ) continue;	//1.2? 1.5?
		if( hit<3 ) continue;
		if( chi2>4) continue;

		ntrk[arm]++;

		float DCA_x_fvtx = (fvx-(tan(the)*cos(phi))*(fvz-fvtx_z)-fvtx_x)*cos(phi);	//Track matching at the Z-vtx, FVTX-z
		float DCA_y_fvtx = (fvy-(tan(the)*sin(phi))*(fvz-fvtx_z)-fvtx_y)*sin(phi);
		//float DCA_r_fvtx = DCA_x_fvtx + DCA_y_fvtx;
		float DCA_x_BeMean = DCA_x_fvtx;
                float DCA_y_BeMean = DCA_y_fvtx;

		float DCA_x_BbV = fvx + tan(the)*cos(phi)*(vertex - fvz); //Track matching at the Z-vtx
		float DCA_y_BbV = fvy + tan(the)*sin(phi)*(vertex - fvz);

		//float DCA_x_BeMean = DCA_x_BbV;
		//float DCA_y_BeMean = DCA_y_BbV;

		float mean_Dca_X=-9999;
		float mean_Dca_Y=-9999;
		int ieta=0;

		if(eta<0)	//FVTX south
		{
			mean_Dca_X = 0.0 ;
			mean_Dca_Y = 0.0 ;
			ieta=0;
		}
		else		//FVTX north
		{
			mean_Dca_X = 0.0 ;
			mean_Dca_Y = 0.0 ;
			ieta=1;
		}

                if(mean_Dca_X==-9999 || mean_Dca_Y==-9999) continue;
		DCA_x_BbV = DCA_x_BbV - mean_Dca_X; //Track matching at the Z-vtx
                DCA_y_BbV = DCA_y_BbV - mean_Dca_Y; //Track matching at the Z-vtx

		DCA_x_fvtx = DCA_x_fvtx - mean_Dca_X; //Track matching at the Z-vtx
                DCA_y_fvtx = DCA_y_fvtx - mean_Dca_Y; //Track matching at the Z-vtx

                float fvt_dca_R_xcos_ysin = DCA_x_BbV*cos(phi) + DCA_y_BbV*sin(phi);
		float fvt_dca_R = DCA_x_fvtx + DCA_y_fvtx;

		bool GoodChi = chi2<4;    

		//Fill FVTX track buffer
		if(itrk<1000 && GoodChi)
		{
			m_Fvt_Phi [ic_ana][iz_ana][ieve][itrk] = phi;
			m_Fvt_Eta [ic_ana][iz_ana][ieve][itrk] = eta;
			//m_Fvt_DcaX[ic_ana][iz_ana][ieve][itrk] = DCA_x_BbV;
			//m_Fvt_DcaY[ic_ana][iz_ana][ieve][itrk] = DCA_y_BbV;
			m_Fvt_DcaX[ic_ana][iz_ana][ieve][itrk] = DCA_x_fvtx;
                        m_Fvt_DcaY[ic_ana][iz_ana][ieve][itrk] = DCA_y_fvtx;
			//++itrk;

                        h_Fvt_Dca[ieta][0][0]->Fill(DCA_x_BeMean);
                        h_Fvt_Dca[ieta][1][0]->Fill(DCA_y_BeMean);
			//h_Fvt_Dca[ieta][0][1]->Fill(DCA_x_BbV);  //[ieta][Mean]  [S/N][Mean=0/Mean=#]
			//h_Fvt_Dca[ieta][1][1]->Fill(DCA_y_BbV);
                        h_Fvt_Dca[ieta][0][1]->Fill(DCA_x_fvtx);  //[ieta][Mean]  [S/N][Mean=0/Mean=#]
                        h_Fvt_Dca[ieta][1][1]->Fill(DCA_y_fvtx);

			if (fabs(fvt_dca_R) < 2)
			{
				//h_Fvt_Dca	[ieta][0][2]->Fill(DCA_x_BbV);
				//h_Fvt_Dca	[ieta][1][2]->Fill(DCA_y_BbV);
                                h_Fvt_Dca	[ieta][0][2]->Fill(DCA_x_fvtx);
                                h_Fvt_Dca	[ieta][1][2]->Fill(DCA_y_fvtx);
				++itrk;

				//compare DCA_r BBC & FVTX
				if (fabs(fvt_dca_R_xcos_ysin)<2) {
					h_BbV_Fvt_Dca	[ic_ana][ieta]->Fill(fvt_dca_R_xcos_ysin, fvt_dca_R);
					h_BbV_Fvt_z	[ic_ana][ieta]->Fill(vertex, fvtx_z);
				}
			}

			h_Fvt_phi_eta[iz_ana][0]->Fill(phi,eta);
			if (fabs(fvt_dca_R) < 2)  h_Fvt_phi_eta[iz_ana][1]->Fill(phi,eta);

			//?????????
			else if(fabs(fvt_dca_R) > 2) continue;

		}// itrk <1000

	} // while
	//Fill number of FVTX tracks
	m_Fvt_Ntrk[ic_ana][iz_ana][ieve] = itrk;  

	//Fvtx track distribution
	hFvtxTrk[0]->Fill(ntrk[0]);
	hFvtxTrk[1]->Fill(ntrk[1]);
	hFvtxTrk[2]->Fill(ntrk[0]+ntrk[1]);

	if(itrk>1000) cout << " N trk  "<< itrk << " "  << ic_ana << " " << iz_ana << endl;   
	bool status = true;
	return status;
}


//When use FVTX centrality bin
int pp_Cor::GetFVTXCentralityBin(TFvtxCompactTrkMap *trkfvtx_map)
{
  
	int centbin  = -9999;
	if ( !trkfvtx_map ) return centbin;

	int ntrk[2] = {0};
  
	//Get FVTXCompactTrkMap
	TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
	//Loop over FVTX tracks in an event
	while (TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next()) {
		TFvtxCompactTrk* fvtx_trk = trk_ptr->get();
		float the = fvtx_trk->get_fvtx_theta();
		float eta = fvtx_trk->get_fvtx_eta();
		float phi = fvtx_trk->get_fvtx_phi();
		float hit = (float)fvtx_trk->get_nhits();
		int   arm = (int)fvtx_trk->get_arm();
		float fvx = fvtx_trk->get_fvtx_vtx().getX();
		float fvy = fvtx_trk->get_fvtx_vtx().getY();
		float fvz = fvtx_trk->get_fvtx_vtx().getZ();
		float chi2 = fvtx_trk->get_chi2_ndf();
		if(the==0 && phi==0 && fvx==0 && fvy==0 && fvz==0) continue;
		if((arm==0 && eta>0) || (arm==1 && eta<0) || arm<0 || arm>1) continue;

		//FVTX track condition
		if( fabs(eta)> 2.5 ) continue;
		if( fabs(eta)< 1.5 ) continue; //1.5? 1.2?
		if( hit<3 ) continue; //Why?
                if( chi2>4) continue;

		//
		ntrk[arm]++;
	} // while

	//if ( (ntrk[0]+ntrk[1])==0 ) return centbin;		//why?
	
	//FVTX multiplicity binning for pp MB
	if ( (ntrk[0]+ntrk[1])>11.5 ) centbin = 0;		//0-1%
	else if ( (ntrk[0]+ntrk[1])>8.5 ) centbin = 1;		//1-5%
	else if ( (ntrk[0]+ntrk[1])>5.5 ) centbin = 2;		//5-20%
	else	 			  centbin = 3;		//20-100%

	//FVTX multiplicity binning for pp FVTX
	//if	( (ntrk[0]+ntrk[1])>19.5 ) centbin = 0;		//0-1%
	//else if ( (ntrk[0]+ntrk[1])>15.5 ) centbin = 1;		//1-5%
	//else if ( (ntrk[0]+ntrk[1])>11.5 ) centbin = 2;		//5-20%
	//else				   centbin = 3;		//20-100%

	//FVTX multiplicity binning for pp FVTX datasets
	//if      ( (ntrk[0]+ntrk[1])>22.5 ) centbin = 0;         //0-1%
        //else if ( (ntrk[0]+ntrk[1])>17.5 ) centbin = 1;         //1-5%
        //else if ( (ntrk[0]+ntrk[1])>12.5 ) centbin = 2;         //5-20%
        //else				   centbin = 3;		//20-100%

	fvtx_mult[0] = ntrk[0];
	fvtx_mult[1] = ntrk[1];
	fvtx_mult[2] = ntrk[0]+ntrk[1];

	return centbin;

 }
