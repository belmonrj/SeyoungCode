#include "pp_filter_ep.h"

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
#include <TH1.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <PHPoint.h>

//#include "Run15pppc3dphidzcalibsmoothpass1.h"
#include <TFvtxGlobalParCntrl.h>
#include <TFvtxCompactTrkMap.h>
#include <TFvtxCompactCoordMap.h>
#include <PHCentralTrack.h>
#include <PHSnglCentralTrack.h>
#include <DoubleInteractionUtil.h>

using namespace std;

//==============================================================
pp_filter_ep::pp_filter_ep(const char* output):
	SubsysReco("RPCALIBRATOR"),
	_filename(output),
	HistoManager(new Fun4AllHistoManager("RPCalibHisto")),
	bbccalib(NULL),
	bbcgeo(NULL),
	doubleint_util(NULL),
	runnumber(0)
{

	fx = NULL;
	fy = NULL;

	for (int ii=0; ii<10; ii++){
		h_Fvtx_phi[ii] = NULL;
	}//ii

	for (int icent=0; icent<6; icent++){
		for (int ii=0; ii<10; ii++){
			h_Bbc_qx[icent][ii] = NULL;
			h_Bbc_qy[icent][ii] = NULL;
			h_Bbc_psi[icent][ii] = NULL;
			h_Bbc_psi_cos[icent][ii] = NULL;
			h_Bbc_psi_sin[icent][ii] = NULL;
			h_Cnt_psi[icent][ii] = NULL;
			h_Cnt_qx[icent][ii] = NULL;
			h_Cnt_qy[icent][ii] = NULL;
		}//ii
	}//icent



}

//==============================================================
pp_filter_ep::~pp_filter_ep()
{

}

//==============================================================
int pp_filter_ep::Init(PHCompositeNode *topNode)
{
	cout << "pp_filter_ep::Init started..." << endl;
 
	bbccalib = new BbcCalib();
	bbcgeo   = new BbcGeo();
	doubleint_util = new DoubleInteractionUtil(); //multi vertex event

	//PHTFileServer::get().open(_filename,"RECREATE");
	//cout << "Opening TFile : " << _filename.c_str() << endl;
	
	//Beam center correction p+p
	//fx = new TF1("fx", "-0.00311056*x + 1.68313", -10, 10);
	//fy = new TF1("fy", "0.00175398*x + 0.663009", -10, 10);
	//Beam center correction p+Au
	//fx = new TF1("fx", "0.0330196*x + 2.01022", -10, 10);
	//fy = new TF1("fy", "9.68418e-04*x + 0.80678", -10, 10);
	//Beam center correction p+Al
	fx = new TF1("fx", "0*x + 0", -10, 10);	//0.0260934*x + 2.073
	fy = new TF1("fy", "0*x + 0", -10, 10);	//0.00208715*x + 0.726095

	cout << "pp_filter_ep::Init ended." << endl;
	return 0;
}

//==============================================================
int pp_filter_ep::InitRun(PHCompositeNode *topNode)
{
	cout << "pp_filter_ep::InitRun started..." << endl;

	recoConsts *rc = recoConsts::instance();
	runnumber = rc->get_IntFlag("RUNNUMBER");

	///Calibration BBC PMT
	int icalibversion = 4002;		// Calibration ID(calibrated PMTs are used)
	RunToTime* runTime = RunToTime::instance();
	PHTimeStamp* ts( runTime->getBeginTime(runnumber) );
	PHTimeStamp tstart = *ts;
	delete ts;
	int bbccalib_version = 4002;
	cout << "pp_filter_ep::InitRun, run number= " << runnumber
		 << " Version of BbcCalib = " << bbccalib_version << endl;
	bbccalib->restore(tstart, bbccalib_version);
	rc->set_IntFlag("BBCCALIBVERSION", icalibversion);  

	TFvtxGlobalParCntrl::set_pdb_run_number(runnumber);
	TFvtxGlobalParCntrl::init_run();
	TFvtxGlobalParCntrl::print();

	doubleint_util->setBbcCalib(topNode);

	//Histogram initialization
	initHisto();
	return 0; 
}

//==============================================================
int pp_filter_ep::process_event(PHCompositeNode *topNode)
{

	//Get nodes and Check
	VtxOut *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut"); //Vertex info
	TrigLvl1 *triglvl1 = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1"); //Trigger info
	PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal"); //Global info
	BbcRaw *bbcraw  = findNode::getClass<BbcRaw>(topNode,"BbcRaw"); //BBC info
	//TFvtxCompactTrkMap *trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode, "TFvtxCompactTrkMap"); //FVTX track
	TFvtxCompactCoordMap *hitfvtx_map = findNode::getClass<TFvtxCompactCoordMap>(topNode,"TFvtxCompactCoordMap");	//FVTX cluster
	PHCentralTrack *central = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack"); //Central arm track
	EventHeader *evth = findNode::getClass<EventHeader>(topNode,"EventHeader"); //Event header

	if(!global		) { cout << PHWHERE << "No PHGloval object !"		<< endl; return ABORTEVENT; }
	if(!vtxout		) { cout << PHWHERE << "Could not find VtxOut !"	<< endl; return ABORTEVENT; }
	if(!triglvl1	) { cout << PHWHERE << "Could not find TrigLvl1 !"	<< endl; return ABORTEVENT; }
	if(!hitfvtx_map	) { cout << PHWHERE << "No TFvtxCompactTrkMap object !" << endl; return ABORTEVENT; }
	if(!central		) { cout << PHWHERE << "CNT Event skip "				<< endl; return ABORTEVENT; }
	if(!bbcraw		) { cerr << PHWHERE << " Could not find BbcRaw !"		<< endl; return ABORTEVENT; }
	if(!evth		) { cerr << PHWHERE << " Could not find EventHeader !"	<< endl; return ABORTEVENT; }


	//Check vertex
	float bbc_z = global->getBbcZVertex();
	if ( fabs(bbc_z)>10.0 ) return ABORTEVENT;

	//Check trigger
	int scaled = triglvl1->get_lvl1_trigscaled();
	int itrig = -9999; 
	
	//MB BBCLL1 +-30, novertex, narrowvtx 
	if ( (scaled&0x00000001) || (scaled&0x00000002) || (scaled&0x00000010) || (scaled&0x00000008) ) itrig = 1;
	if ( itrig<0 ) return ABORTEVENT;
	short centrality = short(global->getCentrality()+0.01);
	//if ( centrality<0 || centrality>99) return ABORTEVENT;
	if ( centrality<0 || centrality>30) return ABORTEVENT;
	if ( (scaled&0x00000008) && centrality>5.5 ) return ABORTEVENT;

	PHPoint fvtx_vertex = vtxout->get_Vertex("FVTX");
	//float fvtxx = fvtx_vertex.getX();
	//float fvtxy = fvtx_vertex.getY();
	float fvtxz = fvtx_vertex.getZ();
	//if(TMath::Abs(fvtxx) > 10 || TMath::Abs(fvtxy) > 10 || TMath::Abs(fvtxz) > 30) return ABORTEVENT;

	if ( TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy") ){
		//fvtxx = TFvtxGlobalParCntrl::get_float_par("beam_x_seed");
		//fvtxy = TFvtxGlobalParCntrl::get_float_par("beam_y_seed");
	}

	//int eventnumber = -1;
	//if ( evth )
		//eventnumber = evth->get_EvtSequence();

	//z-bin
	float tmp_z = fvtxz;
	if ( isnan(tmp_z) || fabs(fvtxz)>10.0 ) tmp_z = bbc_z;
	int i_zvtx = (int)(tmp_z+10.0)/2.;
	if ( i_zvtx<0 || i_zvtx>=10 ) return ABORTEVENT;

	//cent-bin
	int i_cent = (int)(centrality-1)/5;
	if ( i_cent<0 || i_cent>=6 ) return ABORTEVENT;

	//Get BBC information
	int nBbcPmt = 0;

	//Get beam angle correction
	float sumCha = 0;
	for (int ipmt=0; ipmt<128; ipmt++){
		short adc = bbcraw->get_Adc(ipmt);
		short tdc = bbcraw->get_Tdc0(ipmt);
		float time0  = bbccalib->getHitTime0(ipmt, tdc, adc);
		float charge = bbccalib->getCharge(ipmt, adc);
		//float bbcx = bbcgeo->getX(ipmt);
		//float bbcy = bbcgeo->getY(ipmt);
		float bbcz = bbcgeo->getZ(ipmt);

		if (time0>0 && charge>0) {
			short iarm = 0;			//S
			if (bbcz > 0) iarm = 1; //N
			//float phi = atan2(bbcy,bbcx);
			//float val = charge;
			//float rad = sqrt(bbcx*bbcx+bbcy*bbcy);
			//float the = atan2(rad,bbcz); // bbcz-bbcv*10.0);
			//float eta = -log(tan(0.5*the));
			//Symmetry
			//if(ipmt==10||ipmt==21||ipmt==42||ipmt==53||ipmt==74||ipmt==85||ipmt==106||ipmt==117) continue;
			if (iarm==0) sumCha = sumCha + charge;
			nBbcPmt++;
		}//time, charge
	}//ipmt

	if( sumCha<=0 )	return ABORTEVENT;
	float doubleint_frac = doubleint_util->calcFrac(topNode);
	if( doubleint_frac<=0.9 ) return ABORTEVENT;
	if( nBbcPmt <= 0 )	return ABORTEVENT;

	//Get BBC RP
	float Qx_bbc = 0, Qy_bbc = 0, Qw_bbc = 0;
	for (int ipmt=0; ipmt<128; ipmt++){
		short adc = bbcraw->get_Adc(ipmt);
		short tdc = bbcraw->get_Tdc0(ipmt);
		float time0  = bbccalib->getHitTime0(ipmt, tdc, adc);
		float charge = bbccalib->getCharge(ipmt, adc);
		float bbcx = bbcgeo->getX(ipmt);
		float bbcy = bbcgeo->getY(ipmt);
		float bbcz = bbcgeo->getZ(ipmt);

		if (time0>0 && charge>0 && bbcz<0){
			float phi = atan2(bbcy, bbcx);

			Qx_bbc += charge * cos(2*phi);
			Qy_bbc += charge * sin(2*phi);
			Qw_bbc += charge;

			//h_Bbc_phi[i_zvtx]->Fill(ipmt, charge);
		}//time, charge
	}//ipmt

	if ( Qw_bbc>0 ){

		h_Bbc_qx[i_cent][i_zvtx]->Fill(Qx_bbc/Qw_bbc);
		h_Bbc_qy[i_cent][i_zvtx]->Fill(Qy_bbc/Qw_bbc);

		float psi_bbc = atan2(Qy_bbc, Qx_bbc)/2;
		if ( psi_bbc<0 ) psi_bbc += 2*M_PI;
		h_Bbc_psi[i_cent][i_zvtx]->Fill(2*psi_bbc);
		h_Bbc_psi_sin[i_cent][i_zvtx]->Fill(sin(2*psi_bbc));
		h_Bbc_psi_cos[i_cent][i_zvtx]->Fill(cos(2*psi_bbc));
	}

	//Get FVTX track information
	int nFvtxhit = 0;
	TFvtxCompactCoordMap::iterator _iter( hitfvtx_map->range() );
	while( TFvtxCompactCoordMap::const_pointer fvtx_ptr = _iter.next() )
    {
		TFvtxCompactCoord* fvtx_coord = fvtx_ptr->get();
		PHPoint fvtx_coord_point = fvtx_coord->get_coord_midpoint();
		
		int iarm = fvtx_coord->get_arm();
		float fvtx_x = fvtx_coord_point.getX();
		float fvtx_y = fvtx_coord_point.getY();
		float fvtx_z = fvtx_coord_point.getZ();
		//float fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
		if( (fabs(fvtx_x)>999) ||(fabs(fvtx_y)>999) || (fabs(fvtx_z)>999)) continue;
		//float fvtx_the = atan2(fvtx_r,fvtx_z-fvtx_vertex.getZ());
		float fvtx_phi = atan2(fvtx_y,fvtx_x);
		//float fvtx_eta = -log(tan(0.5*fvtx_the));
		if(iarm!=0) continue;

		//fvtx_x = fvtx_x - fvtx_vertex.getX();
		//fvtx_y = fvtx_y - fvtx_vertex.getY();
		//fvtx_phi = atan2(fvtx_y, fvtx_x);

		//float beam_angle_xz = 0.0023;
		//float beam_angle_yz = 2e-4;
		//fvtx_phi = atan2((fvtx_y - fvtx_z*tan(beam_angle_yz))*cos(beam_angle_yz), (fvtx_x - fvtx_z*tan(beam_angle_xz))*cos(beam_angle_xz));
		
		h_Fvtx_phi[i_zvtx]->Fill(fvtx_phi);
		nFvtxhit++;
	} // while

	//Get central arm track information
	int nCntTrack = 0;

	float Qx_cnt = 0, Qy_cnt = 0, Qw_cnt = 0;
	if ( central->get_npart()<1000 ){
		for(unsigned int itrk=0;itrk<central->get_npart();itrk++)
		{
			PHSnglCentralTrack *d_cnt = central->get_track(itrk);

			float zed	 = d_cnt->get_zed();
			float mom	 = d_cnt->get_mom();
			float the0   = d_cnt->get_the0();
			short ch	 = d_cnt->get_charge();
			float eta0   = -log(tan(the0/2.0)); 
			if(fabs(eta0)>0.35) continue;
			float phi0   = d_cnt->get_phi0();  
			phi0		 = atan2(sin(phi0), cos(phi0)); 	  
			//float xpc1	 = d_cnt->get_ppc1x();
			//float ypc1	 = d_cnt->get_ppc1y();
			//float xpc3	 = d_cnt->get_ppc3x();
			//float ypc3	 = d_cnt->get_ppc3y();

			//float phi1 = atan2(ypc1,xpc1);
			//float phi3 = atan2(ypc3,xpc3);

			float zpc1 = d_cnt->get_ppc1z();
			float zpc3 = d_cnt->get_ppc3z();
			float pt   = mom*sin(the0);
			float pc3dphi  = d_cnt->get_pc3dphi();
			float pc3dz	   = d_cnt->get_pc3dz();
			//float ecorr = d_cnt->get_ecore();
			float sdphi = -9999;
			float sdz   = -9999;
			short arm	= -100;
			arm = d_cnt->get_dcarm();

			//Central arm track quality cut
			bool good_trk = mom<20.0 && mom>0.0 && fabs(pc3dphi)<5.0 && fabs(pc3dz)<20.0 && (d_cnt->get_quality()==31 || d_cnt->get_quality()==63) && fabs(zed)<70;

			sdphi = calcsdphi(pc3dphi, arm, int(ch), mom);
			sdz   = calcsdz  (pc3dz,   arm, int(ch), mom);

			if ( !good_trk ) continue;
			if ( fabs(sdphi)>3.0 || pc3dphi<-1000 ) continue;
			if ( fabs(sdz)>3.0 || pc3dz<-1000 ) continue;
			if ( pt<0.2 || pt>2.0 ) continue;
			//if ( fabs(pc3dphi)>0.1 ) continue;
			//if ( fabs(pc3dz)>20. ) continue;
			if ( zpc1<-900 ) continue;
			if ( zpc3<-900 ) continue;

			Qx_cnt += pt * cos(2*phi0);
			Qy_cnt += pt * sin(2*phi0);
			Qw_cnt += pt;

			nCntTrack++;

		} //CUTS on momentum, pc3sdz && pc3sdphi < 10, quality, zed, n0
	}//if

	if ( Qw_cnt>0 ){

		h_Cnt_qx[i_cent][i_zvtx]->Fill(Qx_cnt/Qw_cnt);
		h_Cnt_qy[i_cent][i_zvtx]->Fill(Qy_cnt/Qw_cnt);

		float psi_cnt = atan2(Qy_cnt, Qx_cnt); 
		if ( psi_cnt<0 ) psi_cnt += 2*M_PI;
		h_Cnt_psi[i_cent][i_zvtx]->Fill(psi_cnt);
	}

	return EVENT_OK;
   
}

//==============================================================
int pp_filter_ep::End(PHCompositeNode *topNode)
{
  cout << "pp_filter_ep::End" << endl;

	HistoManager->dumpHistos(_filename);

	delete bbccalib;
	delete bbcgeo;
	delete doubleint_util;
	delete HistoManager;

	cout << "pp_filter_ep::End" << endl;

  return 0; 
}


//==============================================================p+p
/*double pp_filter_ep::calcsdphi(double dphi, int arm, int ch, double mom){
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

double pp_filter_ep::calcsdz(double dz, int arm, int ch, double mom){
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
}
*/

//=============================================================p+Au
double pp_filter_ep::calcsdphi(double dphi, int arm, int ch, double mom){
	TF1 *fm = new TF1("fm","[0] + [1]*x + [2]/x + [3]/sqrt(x) + [4]/(x*x) + [5]/(x*x*x) + [6]/(x*x*x*x) + [7]/(x*x*x*x*x)",0.25,5.0);
	TF1 *fs = new TF1("fs","[0]/(x*x) + [1]/x + [2]/sqrt(x) + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]/(x*x*x)",0.25,5.0);
	int iarm, ich;
	if(arm==0) iarm=0;
	else if(arm==1) iarm=1;
	else return -9999;

	if(ch<0) ich = 0;
	else if(ch>0) ich = 1;
	else return -9999;
	if(iarm==0 && ich==0)
		fm->SetParameters(-0.0274505,0.00127896,-0.0746522,0.0798844,0.0313249,-0.0119987,0.00270235,-0.000247331);
	else if(iarm==0 && ich ==1)
		fm->SetParameters(-0.00261438, 0.000455341, -0.00299343, 0.00639011, -0.00078843, 0.00102687, -0.000465706, 6.20467e-05);
	else if(iarm==1 && ich ==0)
		fm->SetParameters(-0.0629548, 0.00264491, -0.164245, 0.175842, 0.0680185, -0.0254179, 0.00544956, -0.000473159);
	else if(iarm==1 && ich ==1)
		fm->SetParameters(-0.0274505,0.00127896,-0.0746522,0.0798844,0.0313249,-0.0119987,0.00270235,-0.000247331);

	double dphimean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(0.0118113, -0.0748528, 0.126266, -0.0680178, 0.00893665, -0.000645881, -1.4843e-05, -0.00103625 );
	else if(iarm==0 && ich ==1)
		fs->SetParameters(0.00598454, -0.0452776, 0.0905342, -0.0593378, 0.0133689, -0.0025327, 0.000223452, -0.000462955 );
	else if(iarm==1 && ich ==0)
		fs->SetParameters(0.014209, -0.0966811, 0.170566, -0.0969038, 0.0148504, -0.00162158, 6.16198e-05, -0.00120976 );
	else if(iarm==1 && ich ==1)
		fs->SetParameters(0.00462378, -0.0186927, 0.0224917, -0.00264353, -0.003049, 0.00113663, -0.000115423, -0.000474393 );
	double dphisigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dphi-dphimean)/dphisigma;
}

double pp_filter_ep::calcsdz(double dz, int arm, int ch, double mom){
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
		fm->SetParameters(0.0509373, -0.929285, 2.53782, -1.86313, 0.792283, -0.216353, 0.0123453, 0.0123453 );
	else if(iarm==0 && ich ==1)
		fm->SetParameters(0.306831, -6.56066, 17.6006, -14.7313, 4.83085, -1.19801, 0.0626286, 0.0626286 );
	else if(iarm==1 && ich ==0)
		fm->SetParameters(0.013713, -0.725383, 2.56501, -0.525718, 1.19284, -0.387513, 0.0243716, 0.0243719 );
	else if(iarm==1 && ich ==1)
		fm->SetParameters(0.22579, -4.56014, 11.9295, -7.8704, 3.13667, -0.787405, 0.0419168, 0.0419168 );

	double dzmean = fm->Eval(mom);

	if(iarm==0 && ich==0)
		fs->SetParameters(-0.341536, 6.58584, -14.8224, 12.836, -3.37253, 0.891513, -0.126594, 0.00705948);
	else if(iarm==0 && ich ==1)
		fs->SetParameters(-0.341028, 7.63176, -19.2791, 18.3979, -6.73693, 2.46926, -0.51596, 0.0451075);
	else if(iarm==1 && ich ==0)
		fs->SetParameters(-0.531353, 11.2687, -28.218, 25.4609, -8.57953, 2.81715, -0.543932, 0.0453912);
	else if(iarm==1 && ich ==1)
		fs->SetParameters(-0.212152, 3.40229, -5.37012, 3.91316, 0.436902, -0.547268, 0.164582, -0.0158912);

	double dzsigma = fs->Eval(mom);
	delete fm;
	delete fs;
	return (dz-dzmean)/dzsigma;
}
/*
//==============================================================p+Al
double pp_filter_ep::calcsdphi(double dphi, int arm, int ch, double mom){
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

double pp_filter_ep::calcsdz(double dz, int arm, int ch, double mom){
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
void pp_filter_ep::initHisto(){
	cout << "Start initHisto.." << endl;
	char hname[300];
	for(int i_zvtx=0; i_zvtx<10; i_zvtx++)
	{
		sprintf(hname, "EP_hFvtx_phi_Z%d",i_zvtx);
		HistoManager->registerHisto(new TH1F(hname,hname,48,-1*M_PI, 1*M_PI));
		h_Fvtx_phi[i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));
	}

	for (int i_cent=0; i_cent<6; i_cent++){
		for(int i_zvtx=0; i_zvtx<10; i_zvtx++)
		{

			sprintf(hname, "EP_hBbc_qx_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,820,-4.1,4.1));
			h_Bbc_qx[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hBbc_qy_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,820,-4.1,4.1));
			h_Bbc_qy[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hBbc_psi_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,360,0,2*M_PI));
			h_Bbc_psi[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hBbc_psi_cos_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,360,-1,1));
			h_Bbc_psi_cos[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hBbc_psi_sin_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,360,-1,1));
			h_Bbc_psi_sin[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hCnt_psi_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,360,0,2*M_PI));
			h_Cnt_psi[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hCnt_qx_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,820,-4.1,4.1));
			h_Cnt_qx[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

			sprintf(hname, "EP_hCnt_qy_C%d_Z%d",i_cent,i_zvtx);
			HistoManager->registerHisto(new TH1F(hname,hname,820,-4.1,4.1));
			h_Cnt_qy[i_cent][i_zvtx] = static_cast<TH1*>(HistoManager->getHisto(hname));

		}
	}

}


