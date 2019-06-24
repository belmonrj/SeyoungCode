#include "pp_filter.h"

#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>
#include <RunNumberRanges.h>
#include <getClass.h>
#include <recoConsts.h>

#include <PHGlobal.h>
#include <PHTFileServer.h>
#include <PreviousEvent.h>
#include <EventHeader.h>

//#include "emcClusterContainer.h"
//#include "emcClusterContent.h"

#include "TrigLvl1.h"
#include "VtxOut.h"
#include <RunToTime.hh>

#include <Bbc.hh>
#include <BbcRaw.h>
#include <BbcCalib.hh>
#include <BbcGeo.hh>
#include <RunToTime.hh>
//#include "BbcMultipleVtx.hh"

//#include "lpcRawv2.h"
//#include "lpcRaw.h"
//#include "lpcRawHitv2.h"
//#include "lpcRawHit.h"

//ROOT Hists
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <PHPoint.h>

//#include "Run15pppc3dphidzcalibsmoothpass1.h"
#include <TFvtxGlobalParCntrl.h>
#include <TFvtxCompactTrkMap.h>
#include <PHCentralTrack.h>
#include <PHSnglCentralTrack.h>
#include <DoubleInteractionUtil.h>

using namespace std;

//==============================================================
//pp_filter::pp_filter() 
pp_filter::pp_filter(const char* output):
	_filename(output),
	T(NULL),
	fx(NULL),
	fy(NULL),
	bbccalib(NULL),
	bbcgeo(NULL),
	doubleint_util(NULL),
	runnumber(0)
{


}

//==============================================================
pp_filter::~pp_filter()
{

}

//==============================================================
int pp_filter::Init(PHCompositeNode *topNode)
{
	cout << "pp_filter::Init started..." << endl;
 
	bbccalib = new BbcCalib();
	bbcgeo   = new BbcGeo();
	doubleint_util = new DoubleInteractionUtil();

	PHTFileServer::get().open(_filename,"RECREATE");
	cout << "Opening TFile : " << _filename.c_str() << endl;
	
	T = new TTree("T","Track tree");
	T->Branch("runnumber",0,"runnumber/I");
	T->Branch("eventnumber",0,"eventnumber/I");
	T->Branch("trigbit_scaled",0,"trigbit_scaled/I");
	T->Branch("centrality",0,"centrality/S");

	T->Branch("bbc_z",0,"bbc_z/F");
	T->Branch("beam_x", 0, "beam_x/F");
	T->Branch("beam_y", 0, "beam_y/F");

	T->Branch("fvtx_x",0,"fvtx_x/F");
	T->Branch("fvtx_y",0,"fvtx_y/F");
	T->Branch("fvtx_z",0,"fvtx_z/F");

	T->Branch("nFvtxTrack",0,"nFvtxTrack/I");
	//T->Branch("FvtxTrack_nhits",0,"FvtxTrack_nhits[nFvtxTrack]/S");
	//T->Branch("FvtxTrack_has_svxhit",0,"FvtxTrack_has_svxhit[nFvtxTrack]/S");
	T->Branch("FvtxTrack_phi",0,"FvtxTrack_phi[nFvtxTrack]/F");
	T->Branch("FvtxTrack_eta",0,"FvtxTrack_eta[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_x",0,"FvtxTrack_x[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_y",0,"FvtxTrack_y[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_z",0,"FvtxTrack_z[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_dca_r",0,"FvtxTrack_dca_r[nFvtxTrack]/F");
	T->Branch("FvtxTrack_dca_x",0,"FvtxTrack_dca_x[nFvtxTrack]/F");
	T->Branch("FvtxTrack_dca_y",0,"FvtxTrack_dca_y[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_pvalue",0,"FvtxTrack_pvalue[nFvtxTrack]/F");
	//T->Branch("FvtxTrack_chi2",0,"FvtxTrack_chi2[nFvtxTrack]/F");

	T->Branch("nCntTrack",0,"nCntTrack/I");
	T->Branch("CntTrack_phi",0,"CntTrack_phi[nCntTrack]/F");
	T->Branch("CntTrack_eta",0,"CntTrack_eta[nCntTrack]/F");
	T->Branch("CntTrack_pt",0,"CntTrack_pt[nCntTrack]/F");
	//T->Branch("CntTrack_phi1",0,"CntTrack_phi1[nCntTrack]/F");
	//T->Branch("CntTrack_phi3",0,"CntTrack_phi3[nCntTrack]/F");
	//T->Branch("CntTrack_zpc1",0,"CntTrack_zpc1[nCntTrack]/F");
	//T->Branch("CntTrack_zpc3",0,"CntTrack_zpc3[nCntTrack]/F");
	//T->Branch("CntTrack_mom",0,"CntTrack_mom[nCntTrack]/F");
	//T->Branch("CntTrack_ch",0,"CntTrack_ch[nCntTrack]/S");
	T->Branch("CntTrack_arm",0,"CntTrack_arm[nCntTrack]/S");
	//T->Branch("CntTrack_pc3dphi",0,"CntTrack_pc3dphi[nCntTrack]/F");
	//T->Branch("CntTrack_pc3dz",0,"CntTrack_pc3dz[nCntTrack]/F");
	T->Branch("CntTrack_pc3sdphi",0,"CntTrack_pc3sdphi[nCntTrack]/F");
	T->Branch("CntTrack_pc3sdz",0,"CntTrack_pc3sdz[nCntTrack]/F");

	T->Branch("nBbcPmt",0,"nBbcPmt/I");
	T->Branch("BbcPmt_arm",0,"BbcPmt_arm[nBbcPmt]/S");
	//T->Branch("BbcPmt_phi",0,"BbcPmt_phi[nBbcPmt]/F");
	T->Branch("BbcPmt_val",0,"BbcPmt_val[nBbcPmt]/F");
	T->Branch("BbcPmt_x",0,"BbcPmt_x[nBbcPmt]/F");
	T->Branch("BbcPmt_y",0,"BbcPmt_y[nBbcPmt]/F");
	T->Branch("BbcPmt_z",0,"BbcPmt_z[nBbcPmt]/F");
	T->Branch("BbcPmt_eta",0,"BbcPmt_eta[nBbcPmt]/F");

	T->Branch("doubleint_frac",0,"doubleint_frac/F");
	T->Branch("doubleint_good",0,"doubleint_good/I");


	//Beam center correction p+p
	//fx = new TF1("fx", "-0.00311056*x + 1.68313", -10, 10);
	//fy = new TF1("fy", "0.00175398*x + 0.663009", -10, 10);
	//Beam center correction p+Au
	////fx = new TF1("fx", "0.0330196*x + 2.01022", -10, 10);
	////fy = new TF1("fy", "9.68418e-04*x + 0.80678", -10, 10);
	//Beam center correction p+Al
	fx = new TF1("fx", "0*x + 0", -10, 10);	//0.0260934*x + 2.073
	fy = new TF1("fy", "0*x + 0", -10, 10);	//0.00208715*x + 0.726095
	cout << "pp_filter::Init ended." << endl;
	return 0;
}

//==============================================================
int pp_filter::InitRun(PHCompositeNode *topNode)
{
  cout << "pp_filter::InitRun started..." << endl;

	recoConsts *rc = recoConsts::instance();
	runnumber = rc->get_IntFlag("RUNNUMBER");

  ///Calibration BBC PMT
  int icalibversion = 4002;		// Calibration ID(calibrated PMTs are used)
  RunToTime* runTime = RunToTime::instance();
  PHTimeStamp* ts( runTime->getBeginTime(runnumber) );
  PHTimeStamp tstart = *ts;
  delete ts;
  int bbccalib_version = 4002;
	cout << "pp_Cor::InitRun, run number= " << runnumber
		   << " Version of BbcCalib = " << bbccalib_version << endl;
  bbccalib->restore(tstart, bbccalib_version);
  rc->set_IntFlag("BBCCALIBVERSION", icalibversion);  

	TFvtxGlobalParCntrl::set_pdb_run_number(runnumber);
	TFvtxGlobalParCntrl::init_run();
	TFvtxGlobalParCntrl::print();

	doubleint_util->setBbcCalib(topNode);
  
  return 0; 
}

//==============================================================
int pp_filter::process_event(PHCompositeNode *topNode)
{

	//Get nodes and Check
	VtxOut *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut"); //Vertex info
	TrigLvl1 *triglvl1 = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1"); //Trigger info
	PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal"); //Global info
	BbcRaw *bbcraw  = findNode::getClass<BbcRaw>(topNode,"BbcRaw"); //BBC info
	TFvtxCompactTrkMap *trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode, "TFvtxCompactTrkMap"); //FVTX track
	PHCentralTrack *central = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack"); //Central arm track
	EventHeader *evth = findNode::getClass<EventHeader>(topNode,"EventHeader"); //Event header

	if(!global      ) { cout << PHWHERE << "No PHGloval object !"		<< endl; return ABORTEVENT; }
	if(!vtxout      ) { cout << PHWHERE << "Could not find VtxOut !"	<< endl; return ABORTEVENT; }
	if(!triglvl1    ) { cout << PHWHERE << "Could not find TrigLvl1 !"	<< endl; return ABORTEVENT; } 
	if(!trkfvtx_map ) { cout << PHWHERE << "No TFvtxCompactTrkMap object !" << endl; return ABORTEVENT; }
	if(!central     ) { cout << PHWHERE << "CNT Event skip "			<< endl; return ABORTEVENT; }    
	if(!bbcraw      ) { cerr << PHWHERE << " Could not find BbcRaw !"	<< endl; return ABORTEVENT; }
	if(!evth		) { cerr << PHWHERE << " Could not find EventHeader !"<< endl; return ABORTEVENT; }


	//Check vertex
	float bbc_z = global->getBbcZVertex();
	if ( fabs(bbc_z)>10.0 ) return ABORTEVENT;

	//Check trigger
	int scaled = triglvl1->get_lvl1_trigscaled();
	int itrig = -9999; 
	
	//MB BBCLL1 +-30, novertex, narrowvtx 
	//if ( (scaled&0x00000001) || (scaled&0x00000002) || (scaled&0x00000010) || (scaled&0x00000008) ) itrig = 1; //p+Au
	if ( (scaled&0x00000010) || (scaled&0x00000008) ) itrig = 1;
	//if ( (scaled&0x00000001) || (scaled&0x00000002) || (scaled&0x00000010) || (scaled&0x00000400) || (scaled&0x00000800)) itrig = 1; //p+p
	//FVTX high mult.
	//if ( Scaled&0x00000400 ) itrig = 1; //FVTX And triggered 
	if ( itrig<0 ) return ABORTEVENT;
	short centrality = short(global->getCentrality()+0.01);
	if ( centrality<0 || centrality>99) return ABORTEVENT;
	//if ( (scaled&0x00000008) && centrality>5.5 ) return ABORTEVENT;


	PHPoint fvtx_vertex = vtxout->get_Vertex("FVTX");
	float fvtx_x = fvtx_vertex.getX();
	float fvtx_y = fvtx_vertex.getY();
	float fvtx_z = fvtx_vertex.getZ();

	PHPoint vertex(0,0,bbc_z);
	/*if ( !isnan(fvtx_z) && fabs(fvtx_z)<15.0 ){
		vertex.setX(fvtx_x);
		vertex.setY(fvtx_y);
		vertex.setZ(fvtx_z);
	}*/

	float beam_x = -9999;
	float beam_y = -9999;

	if ( TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy") ){
		beam_x = TFvtxGlobalParCntrl::get_float_par("beam_x_seed");
		beam_y = TFvtxGlobalParCntrl::get_float_par("beam_y_seed");
		vertex.setX(beam_x);
		vertex.setY(beam_y);
	} // */

	int eventnumber = -1;
	if ( evth )
		eventnumber = evth->get_EvtSequence();

	//Get BBC information
	int nBbcPmt = 0;
	short BbcPmt_arm[128];
	//float BbcPmt_time[128];
	//float BbcPmt_phi[128];
	float BbcPmt_val[128];
	float BbcPmt_x[128];
	float BbcPmt_y[128];
	float BbcPmt_z[128];
	float BbcPmt_eta[128];

	//Get beam angle correction
	//float tmp_z = -9999;
	//if(TMath::IsNaN(fvtx_z)) tmp_z = bbc_z;
	//else tmp_z = fvtx_z;
	//tmp_z=bbc_z;
	//double sumCha = 0;

	for (int ipmt=0; ipmt<128; ipmt++) {
		short adc = bbcraw->get_Adc(ipmt);
		short tdc = bbcraw->get_Tdc0(ipmt);
		float time0  = bbccalib->getHitTime0(ipmt, tdc, adc);
		float charge = bbccalib->getCharge(ipmt, adc);
		float bbcx = bbcgeo->getX(ipmt);
		float bbcy = bbcgeo->getY(ipmt);
		float bbcz = bbcgeo->getZ(ipmt);

		if (time0>0 && charge>0) {
			short iarm = 0;			//S
			if (bbcz > 0) iarm = 1;	//N
			//float phi = atan2(bbcy,bbcx);
			float val = charge;
			float rad = sqrt(bbcx*bbcx+bbcy*bbcy);
			float the = atan2(rad,bbcz); // bbcz-bbcv*10.0);
			float eta = -log(tan(0.5*the));
			//float fab_eta = fabs(eta);
			//Symmetry
			//if(ipmt==10||ipmt==21||ipmt==42||ipmt==53||ipmt==74||ipmt==85||ipmt==106||ipmt==117) continue;
			//BbcPmt_time[nBbcPmt] = time0;
			BbcPmt_arm [nBbcPmt] = iarm;
			//BbcPmt_phi [nBbcPmt] = phi;
			BbcPmt_val [nBbcPmt] = val;
			BbcPmt_x [nBbcPmt] = bbcx;
			BbcPmt_y [nBbcPmt] = bbcy;
			BbcPmt_z [nBbcPmt] = bbcz;
			BbcPmt_eta [nBbcPmt] = eta;

			//sumCha = sumCha + BbcPmt_val[nBbcPmt];

			nBbcPmt++;
		}//time , charge
	}//ipmt

	//if( sumCha>150 )	return ABORTEVENT;
	int doubleint_good = doubleint_util->isEventOK(topNode);
	float doubleint_frac = doubleint_util->calcFrac(topNode);
	//if( doubleint_frac<=0.9 ) return ABORTEVENT;
	if( nBbcPmt <= 0 )	return ABORTEVENT;

	//Get FVTX track information
	int nFvtxTrack = 0;
	//short FvtxTrack_nhits[1000];
	//short FvtxTrack_has_svxhit[1000];
	//float FvtxTrack_the[1000];
	float FvtxTrack_phi[1000];
	float FvtxTrack_eta[1000];
	//float FvtxTrack_dca_r[1000];
	float FvtxTrack_dca_x[1000];
	float FvtxTrack_dca_y[1000];
	//float FvtxTrack_x[1000];
	//float FvtxTrack_y[1000];
	//float FvtxTrack_z[1000];
	//float FvtxTrack_chi2[1000];
	//float FvtxTrack_pvalue[1000];

	TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();

	if ( trk_iter.count()<1000 ){
		//Loop over FVTX tracks in an event
		while (TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next()) {
			TFvtxCompactTrk* fvtx_trk = trk_ptr->get();
			float the = fvtx_trk->get_fvtx_theta();
			float eta = fvtx_trk->get_fvtx_eta();
			float phi = fvtx_trk->get_fvtx_phi();
			short hit = fvtx_trk->get_nhits();
			short arm = fvtx_trk->get_arm();
			float fvx = fvtx_trk->get_fvtx_vtx().getX();
			float fvy = fvtx_trk->get_fvtx_vtx().getY();
			float fvz = fvtx_trk->get_fvtx_vtx().getZ();
			float chi2 = fvtx_trk->get_chi2_ndf();

			if(the==0 && phi==0 && fvx==0 && fvy==0 && fvz==0) continue;
			if((arm==0 && eta>0) || (arm==1 && eta<0) || arm<0 || arm>1) continue;

			//FVTX track condition
			if( fabs(eta)>3 ) continue;
			if( fabs(eta)<1 ) continue;	//1.2? 1.5?
			if( hit<3 ) continue;
			if( chi2>4 || isnan(chi2) || chi2<=0 ) continue;

			float pvalue = TMath::Prob(chi2*(2*hit-5),2*hit-5);
			if ( pvalue<0.05 ) continue;

			//float dca_x = fvx+(tan(the)*cos(phi)*(fvtx_z-fvz));
			//float dca_y = fvy+(tan(the)*sin(phi)*(fvtx_z-fvz));
			float dca_x = fvx+(tan(the)*cos(phi))*(fvz-vertex.getZ())-vertex.getX();
			float dca_y = fvy+(tan(the)*sin(phi))*(fvz-vertex.getZ())-vertex.getY();
			//float dca_r = dca_x*cos(phi) + dca_y*sin(phi);

			//if ( fabs(dca_x)>2 || fabs(dca_y)>2 || fabs(dca_r)>2.0 ) continue;
			//if ( fabs(dca_x)>5 || fabs(dca_y)>5 ) continue;
			//if( dca_r<-2.85 || dca_r>2.85 ) continue;
			if(arm==0) arm=-1;
			//FvtxTrack_the[nFvtxTrack] = the;
			FvtxTrack_phi[nFvtxTrack] = phi;
			FvtxTrack_eta[nFvtxTrack] = eta;
			//FvtxTrack_chi2[nFvtxTrack] = chi2;
			//FvtxTrack_nhits[nFvtxTrack] = hit;
			//FvtxTrack_has_svxhit[nFvtxTrack] = (fvtx_trk->get_svxhit_pattern()>0) ? 1 : 0;
			//FvtxTrack_x[nFvtxTrack] = fvx;
			//FvtxTrack_y[nFvtxTrack] = fvy;
			//FvtxTrack_z[nFvtxTrack] = fvz;
			//FvtxTrack_x[nFvtxTrack] = arm*(20)*tan(the)*cos(phi);
			//FvtxTrack_y[nFvtxTrack] = arm*(20)*tan(the)*sin(phi);
			//FvtxTrack_z[nFvtxTrack] = arm*20;
 
			//FvtxTrack_pvalue[nFvtxTrack] = pvalue;
			//FvtxTrack_dca_r[nFvtxTrack] = dca_r;
			FvtxTrack_dca_x[nFvtxTrack] = dca_x;
			FvtxTrack_dca_y[nFvtxTrack] = dca_y;
			nFvtxTrack++;
		} // while
	}//if


	//Get central arm track information
	int nCntTrack = 0;
	float CntTrack_phi[1000];
	float CntTrack_eta[1000];
	float CntTrack_pt[1000];
	//float CntTrack_phi1[1000];
	//float CntTrack_phi3[1000];
	//float CntTrack_zpc1[1000];
	//float CntTrack_zpc3[1000];
	//float CntTrack_mom[1000];
	//short CntTrack_ch[1000];
	short CntTrack_arm[1000];
	//float CntTrack_pc3dz[1000];
	//float CntTrack_pc3dphi[1000];
	float CntTrack_pc3sdz[1000];
	float CntTrack_pc3sdphi[1000];

	if ( central->get_npart()<1000 ){

		for(unsigned int itrk=0;itrk<central->get_npart();itrk++)
		{
		PHSnglCentralTrack *d_cnt = central->get_track(itrk);

			float zed    = d_cnt->get_zed();
			float mom    = d_cnt->get_mom();
			float the0   = d_cnt->get_the0();
			//short ch     = d_cnt->get_charge();
			float eta0   = -log(tan(the0/2.0)); 
			if(fabs(eta0)>0.35) continue;
			float phi0   = d_cnt->get_phi0();      
			phi0		 = atan2(sin(phi0), cos(phi0)); 	  
			//float xpc1 = d_cnt->get_ppc1x();
			//float ypc1 = d_cnt->get_ppc1y();
			//float xpc3 = d_cnt->get_ppc3x();
			//float ypc3 = d_cnt->get_ppc3y();

			//float phi1 = atan2(ypc1,xpc1);
			//float phi3 = atan2(ypc3,xpc3);

			float n0 = d_cnt->get_n0();

			//float zpc1 = d_cnt->get_ppc1z();
			//float zpc3 = d_cnt->get_ppc3z();
			float pt   = mom*sin(the0);
			float pc3dphi  = d_cnt->get_pc3dphi();
			float pc3dz    = d_cnt->get_pc3dz();
			//float ecorr = d_cnt->get_ecore();
			//float sdphi = -9999;
			//float sdz   = -9999;
			short   arm   = -100;
			arm	    = d_cnt->get_dcarm();

			//Central arm track quality cut
			//bool good_trk = mom<10.0 && mom>0.02 && fabs(pc3dphi)<5.0 && fabs(pc3dz)<10.0 && (d_cnt->get_quality()==31 || d_cnt->get_quality()==63) && fabs(zed)<75 && n0<0;
			bool good_trk = mom<10.0 && mom>0.02 && fabs(pc3dphi)<5.0 && fabs(pc3dz)<20.0 && (d_cnt->get_quality()==31 || d_cnt->get_quality()==63) && fabs(zed)<75;
			//sdphi = calcsdphi(pc3dphi, arm, int(ch), mom);
			//sdz   = calcsdz  (pc3dz,   arm, int(ch), mom);

			//if ( fabs(sdphi)>2.0 ) continue;
			//if ( fabs(sdz)>2.0 ) continue;
			//if ( pt<0.2 || pt>5.0 ) continue;
			if ( pc3dphi<-1000 ) continue;
			if ( pc3dz<-1000 ) continue;
			if ( n0>0 ) continue; 
			//if ( zpc1<-900 ) continue;
			//if ( zpc3<-900 ) continue;
			if ( fabs(zed)>75. ) continue;
			if ( fabs(pc3dphi)>5 ) continue;
			if ( fabs(pc3dz)>20. ) continue;
			if ( !good_trk ) continue;
			if ( fabs(d_cnt->get_pc3sdphi())>3 ) continue;
			if ( fabs(d_cnt->get_pc3sdz())>3 ) continue;

			CntTrack_phi[nCntTrack] = phi0;
			CntTrack_eta[nCntTrack] = eta0;
			CntTrack_pt[nCntTrack]	= pt;
			//CntTrack_phi1[nCntTrack] = phi1;
			//CntTrack_phi3[nCntTrack] = phi3;
			//CntTrack_zpc1[nCntTrack] = zpc1;
			//CntTrack_zpc3[nCntTrack] = zpc3;
			//CntTrack_mom[nCntTrack] = mom;
			//CntTrack_ch[nCntTrack]	= ch;
			CntTrack_arm[nCntTrack]	= arm;
			//CntTrack_pc3dphi[nCntTrack] = pc3dphi;
			//CntTrack_pc3dz[nCntTrack] = pc3dz;
			CntTrack_pc3sdphi[nCntTrack] = d_cnt->get_pc3sdphi();
			CntTrack_pc3sdz[nCntTrack] = d_cnt->get_pc3sdz();

			nCntTrack++;

		} //CUTS on momentum, pc3sdz && pc3sdphi < 10, quality, zed, n0
	}//if

	int index = 0;

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&runnumber);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&eventnumber);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&scaled);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&centrality);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&bbc_z);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&beam_x);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&beam_y);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&fvtx_x);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&fvtx_y);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&fvtx_z);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&nFvtxTrack);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_nhits);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_has_svxhit);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_phi);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_eta);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_x);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_y);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_z);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_dca_r);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_dca_x);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_dca_y);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_chi2);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(FvtxTrack_pvalue);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&nCntTrack);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_phi);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_eta);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_pt);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_phi1);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_phi3);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_zpc1);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_zpc3);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_mom);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_ch);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_arm);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_pc3dphi);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_pc3dz);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_pc3sdphi);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(CntTrack_pc3sdz);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&nBbcPmt);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_arm);
	//((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_phi);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_val);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_x);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_y);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_z);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(BbcPmt_eta);

	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&doubleint_frac);
	((TBranch*) T->GetListOfBranches()->At(index++))->SetAddress(&doubleint_good);

	T->Fill();
	return EVENT_OK;
   
}

//==============================================================
int pp_filter::End(PHCompositeNode *topNode)
{
  cout << "pp_filter::End" << endl;

	PHTFileServer::get().cd( _filename );
	PHTFileServer::get().write( _filename );

	delete T;
	delete bbccalib;
	delete bbcgeo;
	delete doubleint_util;

  return 0; 
}


//==============================================================p+p
/*double pp_filter::calcsdphi(double dphi, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0] + [1]*x + [2]/x + [3]/sqrt(x) + [4]/(x*x) + [5]/(x*x*x) + [6]/(x*x*x*x) + [7]*(x*x*x*x*x)",0.25,6);
	TF1 *fs = new TF1("fs","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]/(x*x*x) + [7]*(x*x*x*x*x)", 0.25, 6);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;
	if(iarm==0 && ich==0)
		fm->SetParameters(0.053953, -0.00293396, 0.118234, -0.138271, -0.0396848, 0.0109908, -0.00125326, 2.68837e-07);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(-0.0263552, 0.00148585, -0.0645708, 0.0734685, 0.0229488, -0.00680684, 0.000807035, 1.70141e-07);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(0.0901418, -0.00509145, 0.199729, -0.235885, -0.0650362, 0.0171501, -0.00187679, 5.62225e-07);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(-0.0345231, 0.00166382, -0.0860921, 0.0945831, 0.0317197, -0.00936153, 0.00110374, 1.61881e-07);

	double dphimean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(0.0200671, -0.116078, 0.188911, -0.100673, 0.0135155, -0.00121969, -0.00201708, 1.24563e-06);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.00635643, -0.038373, 0.0688587, -0.039379, 0.00643928, -0.000669075, -0.0005727, 1.11442e-06);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(0.0166008, -0.0873469, 0.133341, -0.064538, 0.00697485, -0.000451859, -0.0017451, 2.85698e-07);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.00897022, -0.049434, 0.0821664, -0.043506, 0.00618132, -0.000566941, -0.000903271, 8.04322e-07);
	double dphisigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dphi-dphimean)/dphisigma;
}

double pp_filter::calcsdz(double dz, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x", 0.25, 6);
	TF1 *fs = new TF1("fs","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]/(x*x*x*x*x*x)", 0.25, 6);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;
	if(iarm==0 && ich==0)
		fm->SetParameters(0.10684, -1.42693, 2.84825, -1.77799, 0.289802, -0.0235297);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(0.121702, -1.89415, 3.93834, -2.61587, 0.473624, -0.0464364);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(0.00703618, -0.221735, 0.488234, 1.28112, 0.0581354, -0.00960117);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(0.105871, -1.38684, 2.57782, 0.0863292, 0.196061, -0.0171235);

	double dzmean = fm->Eval(mom);
	if(iarm==0 && ich==0)
		fs->SetParameters(-6.83742, 68.1483, -129.733, 83.2065, -14.9444, 2.23058, -0.140185, 0.0123563);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.417558, -4.1306, 10.1396, -5.86591, 1.68618, -0.278388, 0.0208052, -0.000402488);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(1.04376, -11.3705, 25.5366, -16.6522, 4.14993, -0.746033, 0.0589887, -0.00083927);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.239338, -2.29652, 6.73106, -3.67691, 1.31774, -0.252852, 0.021737, -0.000317913);
	double dzsigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dz-dzmean)/dzsigma;
} // */


//=============================================================p+Au
double pp_filter::calcsdphi(double dphi, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]/x + [5]/(x*x) + [6]/(x*x*x)",0.2,10.0);
	TF1 *fs = new TF1("fs","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]/x + [5]/(x*x) + [6]/(x*x*x)",0.2,10.0);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;
	if(iarm==0 && ich==0)
		fm->SetParameters(0.00426717, -0.00202627, 0.000586614, -6.1908e-05, -0.00271755, 0.000944383, -2.39895e-06);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(0.00919223, -0.00409367, 0.00101229, -7.89514e-05, -0.00709587, 0.00271149, -0.000476753);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(-0.00780958, 0.00325533, -0.000802682, 7.0065e-05, 0.00690397, -0.00302744, 0.00056355);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(0.00876433, -0.00443719, 0.00097195, -6.31313e-05, -0.00882012, 0.00364605, -0.000649617);

	double dphimean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(0.00119616, 0.000953903, -0.00035938, 5.5425e-05, 0.000497536, 0.000314669, 3.37492e-05);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.00166549, 0.000246022, -2.9702e-05, 7.93948e-06, 0.000486903, 0.000463383, -1.75067e-05);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(-0.00593436, 0.0041137, -0.000989007, 0.000100075, 0.00843891, -0.00309795, 0.000517634);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.00599454, -0.00217918, 0.000615328, -4.70661e-05, -0.00264616, 0.00162936, -0.000166532);
	double dphisigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dphi-dphimean)/dphisigma;
}

double pp_filter::calcsdz(double dz, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]/x + [5]/(x*x) + [6]/(x*x*x)",0.2,10.0);
	TF1 *fs = new TF1("fs","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]/x + [5]/(x*x) + [6]/(x*x*x)",0.2,10.0);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;

	if(iarm==0 && ich==0)
		fm->SetParameters(0.0843502, 0.155309, -0.045217, 0.00461161, 0.201241, -0.0795091, 0.0125749);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(0.017154, 0.216081, -0.0656767, 0.00687339, 0.15123, -0.00881832, -0.00800743);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(3.65253, -0.911113, 0.228005, -0.0208106, -1.12228, 0.385303, -0.0509514);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(1.79367, 0.108671, -0.0260577, 0.00217082, 0.36894, -0.182597, 0.025292);

	double dzmean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(1.25887, 0.118658, -0.0412419, 0.00777933, 0.247035, 0.165851, -0.0233439);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(1.24663, 0.127391, -0.0446334, 0.00896978, 0.2553, 0.173663, -0.0257924);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(-0.0896254, 0.832907, -0.232166, 0.0272444, 1.57114, -0.360131, 0.0514998);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(1.17657, 0.0745693, -0.0252284, 0.00572068, 0.607418, -0.0119101, 0.00541764);

	double dzsigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dz-dzmean)/dzsigma;
}
// */

/*
//==============================================================p+Al
double pp_filter::calcsdphi(double dphi, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0] + [1]*x + [2]/x + [3]/sqrt(x) + [4]/(x*x) + [5]/(x*x*x) + [6]/(x*x*x*x) + [7]/(x*x*x*x*x)",0.25,5.0);
	TF1 *fs = new TF1("fs","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]/(x*x*x*x*x)",0.25,5.0);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;
	if(iarm==0 && ich==0)
		fm->SetParameters(0.00552569, -0.000427181, 0.00392857, -0.00856794, 0.00152166, -0.00164362, 0.000634552, -7.88276e-05);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(-0.0107069, 0.00115555, -0.00939008, 0.0209184, -0.00348523, 0.00390584, -0.00149888, 0.000196745);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(-0.0665314, 0.00226561, -0.21827, 0.208051, 0.113662, -0.0522445, 0.0134901, -0.00141195);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(-0.058099, 0.00369886, -0.101342, 0.135845, 0.0197662, 0.001219, -0.00221716, 0.000372653);

	double dphimean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(0.00771713, -0.0757092, 0.150449, -0.0965422, 0.019364, -0.00322469, 0.000247463, -2.41301e-05);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.00124666, -0.0119739, 0.0300937, -0.0222497, 0.00655846, -0.00148368, 0.000156174, -2.35195e-06);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(0.0112793, -0.112046, 0.220994, -0.140718, 0.0274988, -0.00448121, 0.000341311, -3.75452e-05);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.000283973, 0.000848869, 0.000738514, 0.000538291, 0.000529033, -0.000223204, 4.97976e-05, 7.66201e-07);
	double dphisigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dphi-dphimean)/dphisigma;
}

double pp_filter::calcsdz(double dz, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]*x*x*x",0.25,5.0);
	TF1 *fs = new TF1("fs","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]*x*x*x*x",0.25,5.0);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;

	if(iarm==0 && ich==0)
		fm->SetParameters(-0.197932, 3.25742, -8.30674, 6.99787, -2.34848, 0.717242, -0.125499, 0.00901983);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(-0.337141, 4.61928, -10.6958, 8.17767, -2.29435, 0.589494, -0.0871297, 0.00534563);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(-0.181183, 3.65565, -10.4313, 11.6002, -3.83973, 1.29894, -0.248149, 0.0193011);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(0.196526, -4.76172, 13.1941, -10.2478, 4.76665, -1.67144, 0.322439, -0.0246976);

	double dzmean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(-0.310759, 6.02966, -13.5684, 11.9593, -3.29368, 0.955956, -0.160155, 0.0118375);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.162592, -2.76408, 9.29557, -7.50567, 3.30779, -1.03811, 0.183357, -0.0125307);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(0.0988538, -0.0925221, 0.806753, 1.01449, -0.286874, 0.190513, -0.0481125, 0.00529243);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.57826, -10.8832, 31.0691, -26.5573, 10.281, -3.32447, 0.594462, -0.041857);

	double dzsigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dz-dzmean)/dzsigma;
}*/
//=====================================================================
double pp_filter::getBBCX(float zvtx){
	double x_corrected = fx->Eval(zvtx);
	return x_corrected;
}
double pp_filter::getBBCY(float zvtx){
	double y_corrected = fy->Eval(zvtx);
	return y_corrected;
}

