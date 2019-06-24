#include "pAu_Cor.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <TF1.h>

using namespace std;

//Event process
void GetFvtxTrack	(struct EVNT* cor);
void GetBbcPmt		(struct EVNT* cor);
void GetCntTrack	(struct EVNT* cor);
double getBBCX		(float bbc_z);
double getBBCY		(float bbc_z);
void GetBbcIndex	();
int  SetBbcIndex	(int Bbc_arm_, float Bbc_x_, float Bbc_y_);
//void GetRunInfo		(struct EVNT* cor);
//Calculation process
void do_Real_Mix_Fvtx_Cnt(int cent, int zvtx);
void do_Real_Mix_Bbc_Cnt (int cent, int zvtx);
void do_Real_Mix_Fvtx_S_N(int cent, int zvtx);
void do_Real_Mix_Bbc_S_N (int cent, int zvtx);
void do_Real_Mix_Bbc_Fvtx_S(int cent, int zvtx);
void do_Real_Mix_Bbc_Fvtx_N(int cent, int zvtx);
void do_Real_Mix_BbcS_FvtxN(int cent, int zvtx);
void do_Real_Mix_BbcN_FvtxS(int cent, int zvtx);

//void do_Real_Mix_Cnt_Cnt (int cent, int zvtx);
//Buffer related
void do_Buff_Cnt_Fvt_Trk_Reset(int cent, int zvtx);
void Save_event_buffer(struct EVNT* cor, int evt, int cent, int zvtx); //# of event buffer, cent, zvertex
//Histograms
void book_histograms  ();
void write_histograms ();
void book_histograms_QA	();
void write_histograms_QA();

//Beam center correction p+Au
TF1* fx = new TF1("fx", "0.0330196*x + 2.01022", -10, 10);
TF1* fy = new TF1("fy", "0.000097*x + 0.0860678", -10, 10);
//Beam center correction p+Al
//TF1* fx = new TF1("fx", "0.0260934*x + 2.073", -10, 10); //0.0260934*x + 2.073
//TF1* fy = new TF1("fy", "0.00208715*x + 0.726095", -10, 10); //0.00208715*x + 0.726095
//
//TF1* fx = new TF1("fx", "0*x+0.0", -10, 10);
//TF1* fy = new TF1("fy", "0*x+0.0", -10, 10);

//===================================================================================================
void pAu_Cor(const char *fname="infile.root",const char *oname="outfile.root") //RunSingle;
{
  TFile	*infile	= new TFile(fname, "read");
  TTree	*T	= (TTree*)infile->Get("T");

  T->SetBranchAddress("runnumber", &runnumber);
  T->SetBranchAddress("eventnumber", &eventnumber);
  T->SetBranchAddress("trigbit_scaled", &trigbit_scaled);
  T->SetBranchAddress("centrality", &centrality);
  T->SetBranchAddress("beam_x", &beam_x);
  T->SetBranchAddress("beam_y", &beam_y);
  T->SetBranchAddress("bbc_z", &bbc_z);
  T->SetBranchAddress("fvtx_x", &fvtx_x);
  T->SetBranchAddress("fvtx_y", &fvtx_y);
  T->SetBranchAddress("fvtx_z", &fvtx_z);

  T->SetBranchAddress("nFvtxTrack", &nFvtxTrack);
  T->SetBranchAddress("FvtxTrack_phi", FvtxTrack_phi);
  T->SetBranchAddress("FvtxTrack_eta", FvtxTrack_eta);
  T->SetBranchAddress("FvtxTrack_dca_x", FvtxTrack_dca_x);
  T->SetBranchAddress("FvtxTrack_dca_y", FvtxTrack_dca_y);

  T->SetBranchAddress("nCntTrack", &nCntTrack);
  T->SetBranchAddress("CntTrack_phi", CntTrack_phi);
  T->SetBranchAddress("CntTrack_eta", CntTrack_eta);
  T->SetBranchAddress("CntTrack_pt", CntTrack_pt);
  T->SetBranchAddress("CntTrack_arm", CntTrack_arm);
  T->SetBranchAddress("CntTrack_pc3sdphi", CntTrack_pc3sdphi);
  T->SetBranchAddress("CntTrack_pc3sdz", CntTrack_pc3sdz);

  T->SetBranchAddress("nBbcPmt", &nBbcPmt);
  T->SetBranchAddress("BbcPmt_arm", BbcPmt_arm);
  T->SetBranchAddress("BbcPmt_val", BbcPmt_val);
  T->SetBranchAddress("BbcPmt_x",BbcPmt_x);
  T->SetBranchAddress("BbcPmt_y",BbcPmt_y);
  T->SetBranchAddress("BbcPmt_z",BbcPmt_z);
  T->SetBranchAddress("BbcPmt_eta", BbcPmt_eta);

  T->SetBranchAddress("doubleint_frac", &doubleint_frac);
  T->SetBranchAddress("doubleint_good", &doubleint_good);

  int nentries = T->GetEntries();
  cout << "Running	: " << fname << "  :  " << nentries << endl;

  //Book histograms
  book_histograms();
  book_histograms_QA();
  //GetBbcIndex();
  //buffer reset
  for (int icent=0; icent<NCENTANA; icent++){
	for (int iz=0; iz<NZANA; iz++){
		do_Buff_Cnt_Fvt_Trk_Reset(icent,iz);
	}
  }

  //Main loop : event process : Start COR
  for(int ien=0; ien<nentries; ien++)
  {
	if ( fmod(ien,nentries/10)==0 ){ cout << ien*100./nentries << " percent done" << endl; }
	T->GetEntry(ien);
	struct EVNT COR;
	
	//Get informations
	GetFvtxTrack	(&COR);
	GetBbcPmt	(&COR);
	GetCntTrack	(&COR);

	//Event cut
	if(COR.BBC.rawPmt== -9999)	{ cout << eventnumber << endl; continue;}
	if(nBbcPmt <= 0)			{ cout << eventnumber << endl; continue;}//BbcPmt cut
	if(fabs(bbc_z) > 10.0)		{ cout << eventnumber << endl; continue;}
	if(doubleint_frac<=0.9)	continue;
	//if((trigbit_scaled&0x00000002)==0)	{continue;} //Using BBCLL1 novtx
	//if((trigbit_scaled&0x00000001)==0)	{continue;} //Using BBCLL1
	//if((trigbit_scaled&0x00000010)==0)	{continue;} //Using BBCLL1 narrow vtx
	//if((trigbit_scaled&0x00000008)==0 && centrality>5)		{continue;} //Using trigger 0x00008
	//if((trigbit_scaled&0x00000008)!=0) continue;	//not using central trigger
	//if(COR.CNT.track<=0) continue;


	float sumCha=0;
	for(int gi_trk=0; gi_trk<COR.BBC.rawPmt; gi_trk++)
	{
	  if(COR.BBC.arm[gi_trk]==0)  sumCha=COR.BBC.cha[gi_trk]+sumCha;
	}

	//Set centrality
	/*if(centrality>0 && centrality<=1) COR.cent=0;
	else if(centrality>1 && centrality<=3) COR.cent=1;
	else if(centrality>3 && centrality<=5) COR.cent=2;
	else if(centrality>5 && centrality<=15) COR.cent=3;
	else if(centrality>15 && centrality<=30) COR.cent=4;
	else if(centrality>30 && centrality<=60) COR.cent=5;
	else if(centrality>60 && centrality<=84) COR.cent=6; // */

	/*if(centrality>0 && centrality<=5) COR.cent=0;
	else if(centrality>5 && centrality<=20) COR.cent=1;
	else if(centrality>20 && centrality<=40) COR.cent=2;
	else if(centrality>40 && centrality<=60) COR.cent=3;
	else if(centrality>60 && centrality<=84) COR.cent=4; // */
	//if(COR.FVTX.goodTrack[0]>60 && sumCha>200)
	//COR.cent=0;

	if(sumCha>=90)	COR.cent=0;
	else if(sumCha>=73 && sumCha<90)	COR.cent=1;
	else if(sumCha>=64 && sumCha<73)	COR.cent=2;
	else if(sumCha>=44 && sumCha<64)	COR.cent=3;
	else if(sumCha>=25 && sumCha<44)	COR.cent=4;
	else if(sumCha>=11 && sumCha<25)	COR.cent=5;
	else if(sumCha>=1  && sumCha<11)	COR.cent=6; // */
	else continue;

	//QA histos
	h_Fvt_trk[COR.cent]->Fill(COR.FVTX.goodTrack[0]); //2=All, 1=North, 0=South
	for(int gi_trk=0; gi_trk<COR.CNT.track; gi_trk++)
	{
	  h_Cnt_pt[COR.cent]->Fill(COR.CNT.pt[gi_trk]);
	}
	h_Bbc_cha[COR.cent]->Fill(sumCha);
	h_Fvt_Bbc[COR.cent]->Fill(COR.FVTX.goodTrack[0], sumCha); //sumCha=south
	for(int gi_trk=0; gi_trk<COR.FVTX.goodTrack[2]; gi_trk++)
	{
	  h_Fvt_dca[COR.cent][COR.FVTX.arm[gi_trk]]->Fill(COR.FVTX.dca_r[gi_trk]);
	  h_Fvt_eta_dca[COR.cent][COR.FVTX.arm[gi_trk]]->Fill(COR.FVTX.dca_r[gi_trk], COR.FVTX.eta[gi_trk]);
	}
	
	//Set z vertex bin
	COR.vertex = bbc_z;	// using fvtx_z; bbc_z??
	COR.zvtx = (int)(COR.vertex+10.0)/2.;
	if((COR.zvtx<0) || (COR.zvtx>NZANA-1)) { cout << eventnumber << endl; continue; }

	//Save to the buffer
	Save_event_buffer(&COR, n_Eve_buff[COR.cent][COR.zvtx], COR.cent, COR.zvtx);

	//Calculate and save Real and Mix Histograms
	n_Eve_buff[COR.cent][COR.zvtx]++;
	if(n_Eve_buff[COR.cent][COR.zvtx] == NEVE) //number of Event mixing
	{
		do_Real_Mix_Fvtx_Cnt(COR.cent, COR.zvtx);
		do_Real_Mix_Bbc_Cnt (COR.cent, COR.zvtx);

		do_Real_Mix_Fvtx_S_N(COR.cent, COR.zvtx);
		do_Real_Mix_Bbc_S_N (COR.cent, COR.zvtx);
		//do_Real_Mix_Cnt_Cnt (COR.cent, COR.zvtx);
		do_Real_Mix_Bbc_Fvtx_S(COR.cent, COR.zvtx);
		do_Real_Mix_Bbc_Fvtx_N(COR.cent, COR.zvtx);
		do_Real_Mix_BbcS_FvtxN(COR.cent, COR.zvtx);
		do_Real_Mix_BbcN_FvtxS(COR.cent, COR.zvtx);// */

		//initialize buffer
		do_Buff_Cnt_Fvt_Trk_Reset(COR.cent, COR.zvtx);	//initialize Save_event_buffer area
	}
  } //End of main loop : event process : end of COR
  //Write histograms
  TFile *outfile = new TFile(oname, "recreate");  

  //write_histograms();
  write_histograms_QA();
  outfile->Close();
}
//===================================================================================================
void GetFvtxTrack(struct EVNT* cor)
{
  int gi_trk[3]={0};
  cor->FVTX.rawTrack = nFvtxTrack;

  for(int i_trk=0; i_trk<nFvtxTrack; i_trk++)
  {
	//if( fabs(FvtxTrack_dca_x[i_trk])>2 || fabs(FvtxTrack_dca_y[i_trk])>2 ) continue;
	FvtxTrack_dca_r[i_trk]=FvtxTrack_dca_x[i_trk]*cos(FvtxTrack_phi[i_trk]) - FvtxTrack_dca_y[i_trk]*sin(FvtxTrack_phi[i_trk]);
	if( fabs(FvtxTrack_dca_r[i_trk])>5 ) continue;

	if( fabs(FvtxTrack_eta[i_trk])>3 || fabs(FvtxTrack_eta[i_trk])<=1 ) continue;

	cor->FVTX.phi[gi_trk[2]]	= FvtxTrack_phi[i_trk];
	cor->FVTX.eta[gi_trk[2]]	= FvtxTrack_eta[i_trk];

	if(FvtxTrack_eta[i_trk]<0)
	{
		cor->FVTX.arm[gi_trk[2]]=0;
		gi_trk[0]++;
	}
	else if(FvtxTrack_eta[i_trk]>0)
	{
		cor->FVTX.arm[gi_trk[2]]=1;
		gi_trk[1]++;
	}
	cor->FVTX.dca_r[gi_trk[2]]=FvtxTrack_dca_r[i_trk];
	
	gi_trk[2]++;
  }	//End of track loop

  for(int i_arm=0; i_arm<NARM+1; i_arm++)
  {
	cor->FVTX.goodTrack[i_arm]=gi_trk[i_arm];
  }
}
//===================================================================================================
void GetBbcPmt(struct EVNT* cor)
{
  int gi_pmt=0;
  for(int i_pmt=0; i_pmt<nBbcPmt; i_pmt++)
  {
	/*float tmp_x=getBBCX(bbc_z);
	float tmp_y=getBBCY(bbc_z);
	float x_B=BbcPmt_x[i_pmt]-tmp_x;
	float y_B=BbcPmt_y[i_pmt]-tmp_y; // */
	//pmt index
	/*int pmt_index=SetBbcIndex(BbcPmt_arm[i_pmt], BbcPmt_x[i_pmt], BbcPmt_y[i_pmt]);
	h_Bbc_cha_dist[pmt_index]->Fill(BbcPmt_val[i_pmt]);
	if(BbcPmt_val[i_pmt]<0.5) continue; // */
	float x_B=0;
	if(BbcPmt_arm[i_pmt]==0) x_B=-0.5;
	else if(BbcPmt_arm[i_pmt]==1) x_B=4.1;
	cor->BBC.rawPmt		= nBbcPmt;
	cor->BBC.arm[gi_pmt]	= BbcPmt_arm[i_pmt];
	cor->BBC.phi[gi_pmt]	= atan2(BbcPmt_y[i_pmt],BbcPmt_x[i_pmt]-x_B);
	//cor->BBC.phi[gi_pmt]    = atan2(y_B,x_B);
	cor->BBC.cha[gi_pmt]	= BbcPmt_val[i_pmt];
	cor->BBC.eta[gi_pmt]	= BbcPmt_eta[i_pmt];
	gi_pmt++;

	//if(BbcPmt_val[i_pmt]<0.5) cout << BbcPmt_val[i_pmt] << endl;
  }
  cor->BBC.rawPmt=gi_pmt;
}
//===================================================================================================
void GetCntTrack(struct EVNT* cor)
{
  int gi_trk=0;  //good track index
  for(int i_trk=0; i_trk<nCntTrack; i_trk++)
  {
	//if( (CntTrack_pt[i_trk]>5) || (CntTrack_pt[i_trk]<0.2) )	continue; // not cutted :? 
	if( fabs(CntTrack_pc3sdphi[i_trk])>2 || fabs(CntTrack_pc3sdz[i_trk])>2 ) continue;
	//if( CntTrack_arm[i_trk] == 0 ) continue; //0=E,1=W

	cor->CNT.phi[gi_trk]	= CntTrack_phi[i_trk];
	cor->CNT.eta[gi_trk]	= CntTrack_eta[i_trk];
	cor->CNT.pt[gi_trk]		= CntTrack_pt[i_trk];
	//cor->CNT.phi1[gi_trk]	= CntTrack_phi1[i_trk];
	//cor->CNT.phi3[gi_trk]	= CntTrack_phi3[i_trk];
	//cor->CNT.zpc1[gi_trk]	= CntTrack_zpc1[i_trk];
	//cor->CNT.zpc3[gi_trk]	= CntTrack_zpc3[i_trk];

	/*if((0.2 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 0.5))		cor->CNT.ipt_bin[gi_trk]=0;
	else if((0.5 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 1.0 ))   cor->CNT.ipt_bin[gi_trk]=1;
	else if((1.0 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 1.5 ))   cor->CNT.ipt_bin[gi_trk]=2;
	else if((1.5 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 2.0 ))   cor->CNT.ipt_bin[gi_trk]=3;
	else if((2.0 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 3.0 ))   cor->CNT.ipt_bin[gi_trk]=4;
	else if((3.0 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 5.0 ))   cor->CNT.ipt_bin[gi_trk]=5;
	else continue;// */

	/*if((0.2 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 1.0))			cor->CNT.ipt_bin[gi_trk]=0;
	else if((1.0 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 2.0))	cor->CNT.ipt_bin[gi_trk]=1;
	else if((2.0 <= cor->CNT.pt[gi_trk]) && (cor->CNT.pt[gi_trk] < 5.0))	cor->CNT.ipt_bin[gi_trk]=2;
	else continue; // */

	//combining pT
	cor->CNT.ipt_bin[gi_trk]=0;
	gi_trk++;
  }
  cor->CNT.track=gi_trk;

}
//===================================================================================================
//void GetRunInfo(int ien, struct EVNT* cor)
//{
//  cor->RUN.Run=;
//  cor->RUN.bbc_vertex_z	= ;
//  cor->RUN.fvtx_vertex_x	= ;
//  cor->RUN.fvtx_vertex_y	= ;
//  cor->RUN.fvtx_vertex_z	= ;
//  cor->RUN.beam_vertex_x	= ;
//  cor->RUN.beam_vertex_y	= ;
//  cor->RUN.eventnumber	=;
//}
//===================================================================================================
void do_Buff_Cnt_Fvt_Trk_Reset(int cent, int zvtx) //fixed val, fixed val
{
  n_Eve_buff[cent][zvtx] = 0;
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	b_fvtx_trk[cent][zvtx][i_evt] = -9999;
	b_cnt_trk[cent][zvtx][i_evt] = -9999;
	b_bbc_pmt[cent][zvtx][i_evt] = -9999;
	for(int i_trk=0; i_trk<1024; i_trk++)
	{
		b_fvtx_phi[cent][zvtx][i_evt][i_trk] = -9999;
		b_fvtx_eta[cent][zvtx][i_evt][i_trk] = -9999;
		b_fvtx_arm[cent][zvtx][i_evt][i_trk] = -9999;

		b_cnt_phi	[cent][zvtx][i_evt][i_trk] = -9999;
		b_cnt_pt	[cent][zvtx][i_evt][i_trk] = -9999;
		b_cnt_pt_bin	[cent][zvtx][i_evt][i_trk] = -9999;
		b_cnt_eta	[cent][zvtx][i_evt][i_trk] = -9999;

		//b_cnt_zpc1	[cent][zvtx][i_evt][i_trk] = -9999;
		//b_cnt_zpc3	[cent][zvtx][i_evt][i_trk] = -9999;

		b_bbc_phi[cent][zvtx][i_evt][i_trk] = -9999;
		b_bbc_eta[cent][zvtx][i_evt][i_trk] = -9999;
		b_bbc_cha[cent][zvtx][i_evt][i_trk] = -9999;
		b_bbc_arm[cent][zvtx][i_evt][i_trk] = -9999;
	}
  }//End of i_evt loop
}
//===================================================================================================
void Save_event_buffer(struct EVNT* cor, int evt, int cent, int zvtx)
{
  b_fvtx_trk	[cent][zvtx][evt] = cor->FVTX.goodTrack[2];
  for(int i_trk=0; i_trk<b_fvtx_trk[cent][zvtx][evt]; i_trk++)
  {
	b_fvtx_phi	[cent][zvtx][evt][i_trk] = cor->FVTX.phi[i_trk];
	b_fvtx_eta	[cent][zvtx][evt][i_trk] = cor->FVTX.eta[i_trk];
	b_fvtx_arm	[cent][zvtx][evt][i_trk] = cor->FVTX.arm[i_trk];
  }
  b_cnt_trk[cent][zvtx][evt] = cor->CNT.track;
  for(int i_trk=0; i_trk<cor->CNT.track; i_trk++)
  {
	b_cnt_phi	[cent][zvtx][evt][i_trk] = cor->CNT.phi[i_trk];
	b_cnt_pt	[cent][zvtx][evt][i_trk] = cor->CNT.pt [i_trk];
	b_cnt_pt_bin[cent][zvtx][evt][i_trk] = cor->CNT.ipt_bin[i_trk];
	b_cnt_eta	[cent][zvtx][evt][i_trk] = cor->CNT.eta[i_trk];
	//b_cnt_phi1	[cent][zvtx][evt][i_trk] = cor->CNT.phi1[i_trk];
	//b_cnt_phi3	[cent][zvtx][evt][i_trk] = cor->CNT.phi3[i_trk];
	//b_cnt_zpc1	[cent][zvtx][evt][i_trk] = cor->CNT.zpc1[i_trk];
	//b_cnt_zpc3	[cent][zvtx][evt][i_trk] = cor->CNT.zpc3[i_trk];
  }

  b_bbc_pmt[cent][zvtx][evt] = cor->BBC.rawPmt;
  for(int i_pmt=0; i_pmt<cor->BBC.rawPmt; i_pmt++)
  {
	b_bbc_arm[cent][zvtx][evt][i_pmt] = cor->BBC.arm[i_pmt];
	b_bbc_phi[cent][zvtx][evt][i_pmt] = cor->BBC.phi[i_pmt];
	b_bbc_eta[cent][zvtx][evt][i_pmt] = cor->BBC.eta[i_pmt];
	b_bbc_cha[cent][zvtx][evt][i_pmt] = cor->BBC.cha[i_pmt];
  }
}
//===================================================================================================
void do_Real_Mix_Fvtx_Cnt(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_trk_cnt=0; i_trk_cnt<b_cnt_trk[cent][zvtx][i_evt]; i_trk_cnt++)
	{
		if(b_cnt_trk[cent][zvtx][i_evt]==0) continue;
		for(int j_evt=0; j_evt<NEVE; j_evt++) {
  			for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
  			{
				if(b_fvtx_trk[cent][zvtx][j_evt]==0) continue;
				float dPhi = b_cnt_phi[cent][zvtx][i_evt][i_trk_cnt]-b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
				dPhi =  atan2(sin(dPhi),cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;

				float dEta = b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]-b_cnt_eta[cent][zvtx][i_evt][i_trk_cnt]; //+-0.35, 1.5-2.5
				h_Trk_eta[0][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				h_Trk_eta[0][0]->Fill(b_cnt_eta[cent][zvtx][i_evt][i_trk_cnt]);
				if(fabs(dEta)<1.5) 
				{
				h_Trk_eta[0][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				h_Trk_eta[0][1]->Fill(b_cnt_eta[cent][zvtx][i_evt][i_trk_cnt]);
				continue;
				}
				if(i_evt==j_evt)
				{
					h_Cnt_Fvt_dPhi_dEta_Real[cent][zvtx][b_cnt_pt_bin[cent][zvtx][i_evt][i_trk_cnt]]->Fill(dPhi, dEta);
					if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]==0)
						h_Cnt_Fvt_dPhi_Real[cent][zvtx][0]->Fill(b_cnt_pt[cent][zvtx][i_evt][i_trk_cnt],dPhi);
					else if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]==1)
						h_Cnt_Fvt_dPhi_Real[cent][zvtx][1]->Fill(b_cnt_pt[cent][zvtx][i_evt][i_trk_cnt], dPhi);
					else continue;
				}
				else
				{
					h_Cnt_Fvt_dPhi_dEta_Mix[cent][zvtx][b_cnt_pt_bin[cent][zvtx][i_evt][i_trk_cnt]] -> Fill(dPhi, dEta);
					if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]==0)
						h_Cnt_Fvt_dPhi_Mix[cent][zvtx][0] -> Fill(b_cnt_pt[cent][zvtx][i_evt][i_trk_cnt], dPhi);
					else if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]==1)
						h_Cnt_Fvt_dPhi_Mix[cent][zvtx][1] -> Fill(b_cnt_pt[cent][zvtx][i_evt][i_trk_cnt],dPhi);
					else continue;
				}
		}} //FVTX j
  }} //CNT i
}//Function End


void do_Real_Mix_Bbc_Cnt(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
	{
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)	continue;
		for(int j_evt=0; j_evt<NEVE; j_evt++)
		{
			for(int j_trk_cnt=0; j_trk_cnt<b_cnt_trk[cent][zvtx][j_evt]; j_trk_cnt++)
  			{
				if(b_cnt_trk[cent][zvtx][j_evt]==0)		continue;
				float dPhi = b_cnt_phi[cent][zvtx][j_evt][j_trk_cnt]-b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.0*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.0*M_PI;

				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]-b_cnt_eta[cent][zvtx][j_evt][j_trk_cnt];
				float bbcSqr = b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]*b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc];
				h_Trk_eta[1][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[1][0]->Fill(b_cnt_eta[cent][zvtx][j_evt][j_trk_cnt]);
				if(fabs(dEta)<1.5) 
				{
				h_Trk_eta[1][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[1][1]->Fill(b_cnt_eta[cent][zvtx][j_evt][j_trk_cnt]);
				continue;
				}
				if(i_evt==j_evt)
				{
					h_Cnt_Bbc_dPhi_Real[cent][b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]] -> Fill(b_cnt_pt[cent][zvtx][j_evt][j_trk_cnt], dPhi, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
					h_Cnt_Bbc_dPhi_dEta_Real[cent][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]] -> Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
					h_Cnt_Bbc_dPhi_dEta_Real_sqr[cent][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]] -> Fill(dPhi, dEta, bbcSqr);
				}
				else
				{
					h_Cnt_Bbc_dPhi_Mix[cent][b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]] -> Fill(b_cnt_pt[cent][zvtx][j_evt][j_trk_cnt], dPhi, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
					h_Cnt_Bbc_dPhi_dEta_Mix[cent][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]] -> Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
					h_Cnt_Bbc_dPhi_dEta_Mix_sqr [cent][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]] -> Fill(dPhi, dEta, bbcSqr);
				}
		}} //CNT j
	}} //BBC i
}//Function End


void do_Real_Mix_Fvtx_S_N(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_trk_fvtx=0; i_trk_fvtx<b_fvtx_trk[cent][zvtx][i_evt]; i_trk_fvtx++)
		{
			if(b_fvtx_trk[cent][zvtx][i_evt]==0)	continue;
			for(int j_evt=i_evt; j_evt<NEVE; j_evt++)
			{  //To avoid calculate twice
				for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
				{
					if(b_fvtx_trk[cent][zvtx][j_evt]==0)	continue;
					if(b_fvtx_arm[cent][zvtx][i_evt][i_trk_fvtx] - b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]==0 )	continue;	//same arm==continue;
					if((i_evt==j_evt) && (i_trk_fvtx==j_trk_fvtx))	continue;
					float dPhi = b_fvtx_phi[cent][zvtx][i_evt][i_trk_fvtx] - b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
					dPhi = atan2(sin(dPhi), cos(dPhi));
					if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
					if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
					float dEta = b_fvtx_eta[cent][zvtx][i_evt][i_trk_fvtx] - b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx];
					if(dEta > 0) dEta = -1.*dEta;
					h_Trk_eta[2][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
					h_Trk_eta[2][0]->Fill(b_fvtx_eta[cent][zvtx][i_evt][i_trk_fvtx]);
					if(fabs(dEta)<1.5) 
					{
					h_Trk_eta[2][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
					h_Trk_eta[2][1]->Fill(b_fvtx_eta[cent][zvtx][i_evt][i_trk_fvtx]);
					continue;
					}
					if(i_evt==j_evt)
					{
						h_Fvt_S_N_dPhi_dEta_Real[cent][zvtx]->Fill(dPhi, dEta);
						h_Fvt_S_N_dPhi_Real     [cent][zvtx]->Fill(dPhi);
					}
					else
					{
						h_Fvt_S_N_dPhi_dEta_Mix [cent][zvtx]->Fill(dPhi, dEta);
						h_Fvt_S_N_dPhi_Mix      [cent][zvtx]->Fill(dPhi);
					}
			}}//End of j loop
  }}//End of i loop

}//Function End


void do_Real_Mix_Bbc_Fvtx_S(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)				continue;
		if(b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=0)	continue;
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++)
		{
			if(b_fvtx_trk[cent][zvtx][j_evt]==0)		continue;
			for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
			{
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=0)	continue;
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=0 || b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=0)	continue;
				float dPhi = b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx];
				if(dEta > 0) dEta = -1.*dEta;

				h_Trk_eta[3][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				h_Trk_eta[3][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				if(fabs(dEta)<1.5)
				{
				h_Trk_eta[3][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				h_Trk_eta[3][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				continue;
				}

				if(i_evt==j_evt)
				{
					h_Bbc_Fvt_S_dPhi_dEta_Real[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
				else
				{
					h_Bbc_Fvt_S_dPhi_dEta_Mix[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
			}
		}
	}
  }
}//Function End


void do_Real_Mix_Bbc_Fvtx_N(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)				continue;
		if(b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=1)	continue;
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++)
		{
			if(b_fvtx_trk[cent][zvtx][j_evt]==0)		continue;
			for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
			{
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=1)	continue;
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=1 || b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=1)	continue;
				float dPhi = b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx];
				if(dEta < 0) dEta = -1.*dEta;
				h_Trk_eta[4][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[4][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				if(fabs(dEta)<1.5)
				{
				h_Trk_eta[4][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[4][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				continue;
				}

				if(i_evt==j_evt)
				{
					h_Bbc_Fvt_N_dPhi_dEta_Real[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
				else
				{
					h_Bbc_Fvt_N_dPhi_dEta_Mix[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
			}
		}
	}
  }
}//Function End


void do_Real_Mix_BbcS_FvtxN(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)				continue;
		if(b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=0)	continue;
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++)
		{
			if(b_fvtx_trk[cent][zvtx][j_evt]==0)		continue;
			for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
			{
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=1)	continue;
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=1 || b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=0)	continue;
				float dPhi = b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx];
				if(dEta > 0) dEta = -1.*dEta;
				h_Trk_eta[5][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[5][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				if(fabs(dEta)<1.5)
				{
				h_Trk_eta[5][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[5][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				continue;
				}
				if(i_evt==j_evt)
				{
					h_BbcS_FvtN_dPhi_dEta_Real[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
				else
				{
					h_BbcS_FvtN_dPhi_dEta_Mix[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
			}
		}
	}
  }
}//Function End


void do_Real_Mix_BbcN_FvtxS(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++)
  {
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)				continue;
		if(b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=1)	continue;
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++)
		{
			if(b_fvtx_trk[cent][zvtx][j_evt]==0)		continue;
			for(int j_trk_fvtx=0; j_trk_fvtx<b_fvtx_trk[cent][zvtx][j_evt]; j_trk_fvtx++)
			{
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=0)	continue;
				if(b_fvtx_arm[cent][zvtx][j_evt][j_trk_fvtx]!=0 || b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]!=1)	continue;
				float dPhi = b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_phi[cent][zvtx][j_evt][j_trk_fvtx];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc] - b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx];
				if(dEta < 0) dEta = -1.*dEta;
				h_Trk_eta[6][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[6][0]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				if(fabs(dEta)<1.5)
				{
				h_Trk_eta[6][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[6][1]->Fill(b_fvtx_eta[cent][zvtx][j_evt][j_trk_fvtx]);
				continue;
				}
				if(i_evt==j_evt)
				{
					h_BbcN_FvtS_dPhi_dEta_Real[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
				else
				{
					h_BbcN_FvtS_dPhi_dEta_Mix[cent][zvtx]->Fill(dPhi, dEta, b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]);
				}
			}
		}
	}
  }
}//Function End


void do_Real_Mix_Bbc_S_N(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++) {
	for(int i_pmt_bbc=0; i_pmt_bbc<b_bbc_pmt[cent][zvtx][i_evt]; i_pmt_bbc++)
	{
		if(b_bbc_pmt[cent][zvtx][i_evt]==0)		continue;
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++) { //To avoid calculate twice
			for(int j_pmt_bbc=0; j_pmt_bbc<b_bbc_pmt[cent][zvtx][j_evt]; j_pmt_bbc++)
			{
				if(b_bbc_pmt[cent][zvtx][j_evt]==0)				continue;
				if(b_bbc_arm[cent][zvtx][i_evt][i_pmt_bbc]-b_bbc_arm[cent][zvtx][j_evt][j_pmt_bbc] < 1)		continue; // same arm==continue
				if((i_evt==j_evt) && (i_pmt_bbc==j_pmt_bbc))	continue;
				float dPhi = b_bbc_phi[cent][zvtx][i_evt][i_pmt_bbc]-b_bbc_phi[cent][zvtx][j_evt][j_pmt_bbc];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.0*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.0*M_PI;
				float dEta = b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]-b_bbc_eta[cent][zvtx][j_evt][j_pmt_bbc];
				if(dEta > 0) dEta = -1.*dEta;
				h_Trk_eta[7][0]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[7][0]->Fill(b_bbc_eta[cent][zvtx][j_evt][j_pmt_bbc]);
				if(fabs(dEta)<1.5)
				{
				h_Trk_eta[7][1]->Fill(b_bbc_eta[cent][zvtx][i_evt][i_pmt_bbc]);
				h_Trk_eta[7][1]->Fill(b_bbc_eta[cent][zvtx][j_evt][j_pmt_bbc]);
				continue;
				}
				float bbcSqr[3] = {b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]*b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc],
				b_bbc_cha[cent][zvtx][j_evt][j_pmt_bbc]*b_bbc_cha[cent][zvtx][j_evt][j_pmt_bbc],
				b_bbc_cha[cent][zvtx][i_evt][i_pmt_bbc]*b_bbc_cha[cent][zvtx][j_evt][j_pmt_bbc]}; // i2, j2, i*j
				if(i_evt==j_evt)
				{
					h_Bbc_S_N_dPhi_Real [cent]->Fill(dPhi, bbcSqr[2]);
					h_Bbc_S_N_dPhi_dEta_Real    [cent]->Fill(dPhi, dEta, bbcSqr[2]);
					h_Bbc_S_N_dPhi_dEta_Real_sqr[cent]->Fill(dPhi, dEta, bbcSqr[0]+bbcSqr[1]);
				}
				else
				{
					h_Bbc_S_N_dPhi_Mix  [cent]->Fill(dPhi, bbcSqr[2]);
					h_Bbc_S_N_dPhi_dEta_Mix     [cent]->Fill(dPhi, dEta, bbcSqr[2]);
					h_Bbc_S_N_dPhi_dEta_Mix_sqr [cent]->Fill(dPhi, dEta, bbcSqr[0]+bbcSqr[1]);
				}
		}}//End of j
	}}//End of i
}//Function End


/*void do_Real_Mix_Cnt_Cnt(int cent, int zvtx)
{
  for(int i_evt=0; i_evt<NEVE; i_evt++) {
	for(int i_trk_cnt=0; i_trk_cnt<b_cnt_trk[cent][zvtx][i_evt]; i_trk_cnt++)
	{
		if(b_cnt_trk[cent][zvtx][i_evt]==0) { cout << eventnumber << endl; continue;}
		for(int j_evt=i_evt; j_evt<NEVE; j_evt++) { //To avoid calculate twice
			for(int j_trk_cnt=0; j_trk_cnt<b_cnt_trk[cent][zvtx][j_evt]; j_trk_cnt++)
  			{
				if(b_cnt_trk[cent][zvtx][j_evt]==0) { cout << eventnumber << endl; continue;}
				if((i_evt==j_evt) && (i_trk_cnt==j_trk_cnt)) { cout << eventnumber << endl; continue;} //same event & same track=continue

				float dPc1 = b_cnt_phi1[cent][zvtx][i_evt][i_trk_cnt] - b_cnt_phi1[cent][zvtx][j_evt][j_trk_cnt];
				dPc1 = atan2(sin(dPc1),cos(dPc1));
				dPc1 = fabs(dPc1);
				float dPc3 = b_cnt_phi3[cent][zvtx][i_evt][i_trk_cnt] - b_cnt_phi3[cent][zvtx][j_evt][j_trk_cnt];
				dPc3 = atan2(sin(dPc3),cos(dPc3));
				dPc3 = fabs(dPc3);
				float dZpc1 = fabs(b_cnt_zpc1[cent][zvtx][i_evt][i_trk_cnt] - b_cnt_zpc1[cent][zvtx][j_evt][j_trk_cnt]);
				float dZpc3 = fabs(b_cnt_zpc3[cent][zvtx][i_evt][i_trk_cnt] - b_cnt_zpc3[cent][zvtx][j_evt][j_trk_cnt]);
				float Cut1 = sqrt( pow(fabs(dPc1)/0.040, 2.0)  +  pow(fabs( dZpc1 )/90.0, 2.0) );
				float Cut2 = sqrt( pow(fabs(dPc1)/0.080, 2.0)  +  pow(fabs( dZpc1 )/8.00, 2.0) );
				float Cut3 = sqrt( pow(fabs(dPc3)/0.070, 2.0)  +  pow(fabs( dZpc3 )/25.0, 2.0) );
				if( (Cut1<1.) || (Cut2<1.) || (Cut3<1.)) { cout << eventnumber << endl; continue;}

				float dPhi = b_cnt_phi[cent][zvtx][i_evt][i_trk_cnt]-b_cnt_phi[cent][zvtx][j_evt][j_trk_cnt];
				dPhi = atan2(sin(dPhi), cos(dPhi));
				if(dPhi < -0.5*M_PI) dPhi = dPhi + 2.*M_PI;
				if(dPhi >= 1.5*M_PI) dPhi = dPhi - 2.*M_PI;
				float dEta = b_cnt_eta[cent][zvtx][i_evt][i_trk_cnt]-b_cnt_eta[cent][zvtx][j_evt][j_trk_cnt];
  				if(fabs(dEta)<1.5) continue;
				if(i_evt==j_evt)
				h_Cnt_Cnt_dPhi_dEta_Real[cent][b_cnt_pt_bin[cent][zvtx][i_evt][i_trk_cnt]][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]]->Fill(dPhi, dEta);
				else
				h_Cnt_Cnt_dPhi_dEta_Mix [cent][b_cnt_pt_bin[cent][zvtx][i_evt][i_trk_cnt]][b_cnt_pt_bin[cent][zvtx][j_evt][j_trk_cnt]]->Fill(dPhi, dEta);
		}}
	}}
} // */
//===================================================================================================
void book_histograms()
{
  char hname[300];
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	sprintf(hname,"Real_dPhi_BBC_S_N_C%d", i_cent);
	h_Bbc_S_N_dPhi_Real[i_cent] = new TH1D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI);

	sprintf(hname,"Mix_dPhi_BBC_S_N_C%d", i_cent);
	h_Bbc_S_N_dPhi_Mix [i_cent] = new TH1D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI);

	sprintf(hname,"Real_dPhi_dEta_BBC_S_N_C%d", i_cent);
	h_Bbc_S_N_dPhi_dEta_Real[i_cent] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBC_eta,-8,-5);

	sprintf(hname,"Real_dPhi_dEta_BBC_S_N_sqr_C%d", i_cent);
	h_Bbc_S_N_dPhi_dEta_Real_sqr[i_cent] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBC_eta,-8,-5);

	sprintf(hname,"Mix_dPhi_dEta_BBC_S_N_C%d", i_cent);
	h_Bbc_S_N_dPhi_dEta_Mix [i_cent] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBC_eta,-8,-5);

	sprintf(hname,"Mix_dPhi_dEta_BBC_S_N_sqr_C%d", i_cent);
	h_Bbc_S_N_dPhi_dEta_Mix_sqr [i_cent] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBC_eta,-8,-5);

	for(int ipt_bin_t=0; ipt_bin_t<NPTBIN; ipt_bin_t++)
	{
		sprintf(hname,"Real_dPhi_dEta_CNT_BBC_C%d_Pt%d", i_cent, ipt_bin_t);
		h_Cnt_Bbc_dPhi_dEta_Real[i_cent][ipt_bin_t] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-5,5);

		sprintf(hname,"Real_dPhi_dEta_CNT_BBC_sqr_C%d_Pt%d", i_cent, ipt_bin_t);
		h_Cnt_Bbc_dPhi_dEta_Real_sqr[i_cent][ipt_bin_t] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-5,5);

		sprintf(hname,"Mix_dPhi_dEta_CNT_BBC_C%d_Pt%d", i_cent, ipt_bin_t);
		h_Cnt_Bbc_dPhi_dEta_Mix [i_cent][ipt_bin_t] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-5,5);

		sprintf(hname,"Mix_dPhi_dEta_CNT_BBC_sqr_C%d_Pt%d", i_cent, ipt_bin_t);
		h_Cnt_Bbc_dPhi_dEta_Mix_sqr[i_cent][ipt_bin_t] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-5,5);

		/*for(int ipt_bin_a=0; ipt_bin_a<NPTBIN; ipt_bin_a++)
		{
			sprintf(hname, "Real_dPhi_dEta_CNT_CNT_C%d_TpT%d_ApT%d", i_cent, ipt_bin_t, ipt_bin_a);
			h_Cnt_Cnt_dPhi_dEta_Real[i_cent][ipt_bin_t][ipt_bin_a] = new TH2D(hname, hname, H_BBCCNT_eta,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-0.75,0.75);

			sprintf(hname, "Mix_dPhi_dEta_CNT_CNT_C%d_TpT%d_ApT%d", i_cent, ipt_bin_t, ipt_bin_a);
			h_Cnt_Cnt_dPhi_dEta_Mix[i_cent][ipt_bin_t][ipt_bin_a] = new TH2D(hname, hname, H_BBCCNT_eta,-0.5*M_PI,1.5*M_PI, H_BBCCNT_eta,-0.75,0.75);
		}// */
	}

	for(int i_arm=0; i_arm<NARM; i_arm++)
	{
		if(i_arm==0) sprintf(hname,"Real_dPhi_CNT_BBC_C%d_South", i_cent );
		if(i_arm==1) sprintf(hname,"Real_dPhi_CNT_BBC_C%d_North", i_cent );
		h_Cnt_Bbc_dPhi_Real[i_cent][i_arm] = new TH2D(hname, hname, NPTBIN,0,5, H_BBC_phi,-0.5*M_PI,1.5*M_PI);

		if(i_arm==0) sprintf(hname,"Mix_dPhi_CNT_BBC_C%d_South", i_cent);
		if(i_arm==1) sprintf(hname,"Mix_dPhi_CNT_BBC_C%d_North", i_cent);
		h_Cnt_Bbc_dPhi_Mix [i_cent][i_arm] = new TH2D(hname, hname, NPTBIN,0,5, H_BBC_phi,-0.5*M_PI,1.5*M_PI);
	}

	for(int i_zvtx=0; i_zvtx<NZANA; i_zvtx++)
	{
		sprintf(hname,"Real_dPhi_FVT_S_N_C%d_Z%d", i_cent, i_zvtx);
		h_Fvt_S_N_dPhi_Real[i_cent][i_zvtx] = new TH1D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI);

		sprintf(hname,"Mix_dPhi_FVT_S_N_C%d_Z%d", i_cent, i_zvtx);
		h_Fvt_S_N_dPhi_Mix [i_cent][i_zvtx] = new TH1D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI);

		sprintf(hname,"Real_dPhi_dEta_FVT_S_N_C%d_Z%d", i_cent, i_zvtx);
		h_Fvt_S_N_dPhi_dEta_Real[i_cent][i_zvtx] = new TH2D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-6,0);

		sprintf(hname,"Mix_dPhi_dEta_FVT_S_N_C%d_Z%d", i_cent, i_zvtx);
		h_Fvt_S_N_dPhi_dEta_Mix [i_cent][i_zvtx] = new TH2D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-6,0);

		//Same arm combination BBC, FVTX
		sprintf(hname,"Real_dPhi_dEta_BBC_FVT_N_C%d_Z%d", i_cent, i_zvtx);
		h_Bbc_Fvt_N_dPhi_dEta_Real[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,0, 5);

		sprintf(hname,"Mix_dPhi_dEta_BBC_FVT_N_C%d_Z%d", i_cent, i_zvtx);
		h_Bbc_Fvt_N_dPhi_dEta_Mix[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,0, 5);

		sprintf(hname,"Real_dPhi_dEta_BBC_FVT_S_C%d_Z%d", i_cent, i_zvtx);
		h_Bbc_Fvt_S_dPhi_dEta_Real[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-5,0);

		sprintf(hname,"Mix_dPhi_dEta_BBC_FVT_S_C%d_Z%d", i_cent, i_zvtx);
		h_Bbc_Fvt_S_dPhi_dEta_Mix[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-5,0);

		//Different arm combination BBC, FVTX
		sprintf(hname,"Real_dPhi_dEta_BBCS_FVTN_C%d_Z%d", i_cent, i_zvtx);
		h_BbcS_FvtN_dPhi_dEta_Real[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-7.5,-3.5);

		sprintf(hname,"Mix_dPhi_dEta_BBCS_FVTN_C%d_Z%d", i_cent, i_zvtx);
		h_BbcS_FvtN_dPhi_dEta_Mix[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-7.5,-3.5);

		sprintf(hname,"Real_dPhi_dEta_BBCN_FVTS_C%d_Z%d", i_cent, i_zvtx);
		h_BbcN_FvtS_dPhi_dEta_Real[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,3.5,7.5);

		sprintf(hname,"Mix_dPhi_dEta_BBCN_FVTS_C%d_Z%d", i_cent, i_zvtx);
		h_BbcN_FvtS_dPhi_dEta_Mix[i_cent][i_zvtx] = new TH2D(hname, hname, H_BBC_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,3.5,7.5);

		for(int i_arm=0; i_arm<NARM; i_arm++)
		{
			sprintf(hname,"Real_dPhi_CNT_FVT_C%d_Z%d_Arm%d", i_cent, i_zvtx, i_arm);
			h_Cnt_Fvt_dPhi_Real[i_cent][i_zvtx][i_arm] = new TH2D(hname, hname, NPTBIN,0,5, H_FVTX_phi,-0.5*M_PI,1.5*M_PI);

			sprintf(hname,"Mix_dPhi_CNT_FVT_C%d_Z%d_Arm%d", i_cent, i_zvtx, i_arm);
			h_Cnt_Fvt_dPhi_Mix [i_cent][i_zvtx][i_arm] = new TH2D(hname, hname, NPTBIN,0,5, H_FVTX_phi,-0.5*M_PI,1.5*M_PI);
		}

		for(int ipt_bin=0;ipt_bin<NPTBIN;ipt_bin++)
		{
			if(ipt_bin<(NPTBIN/2))
			{
			sprintf(hname,"Real_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin);
			h_Cnt_Fvt_dPhi_dEta_Real[i_cent][i_zvtx][ipt_bin] = new TH2D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-4,4);
  
			sprintf(hname,"Mix_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin);
			h_Cnt_Fvt_dPhi_dEta_Mix[i_cent][i_zvtx][ipt_bin] = new TH2D(hname, hname, H_FVTX_phi,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-4,4);
			}
			else if(ipt_bin>=(NPTBIN/2))
			{
			sprintf(hname,"Real_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin);
			h_Cnt_Fvt_dPhi_dEta_Real[i_cent][i_zvtx][ipt_bin] = new TH2D(hname, hname, H_FVTX_phi/2,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-4,4);

			sprintf(hname,"Mix_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin);
			h_Cnt_Fvt_dPhi_dEta_Mix[i_cent][i_zvtx][ipt_bin] = new TH2D(hname, hname, H_FVTX_phi/2,-0.5*M_PI,1.5*M_PI, H_FVTX_eta,-4,4);
			}
		}
	}//i_zvtx
  }//i_cent
}
//===================================================================================================
void write_histograms()
{
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	h_Bbc_S_N_dPhi_Real[i_cent]->Write();
	h_Bbc_S_N_dPhi_dEta_Real[i_cent]->Write();
	h_Bbc_S_N_dPhi_dEta_Real_sqr[i_cent]->Write();
	h_Bbc_S_N_dPhi_Mix[i_cent]->Write();
	h_Bbc_S_N_dPhi_dEta_Mix[i_cent]->Write();
	h_Bbc_S_N_dPhi_dEta_Mix_sqr [i_cent]->Write();
	for(int ipt_bin_t=0; ipt_bin_t<NPTBIN; ipt_bin_t++)
	{
		h_Cnt_Bbc_dPhi_dEta_Real[i_cent][ipt_bin_t]->Write();
		h_Cnt_Bbc_dPhi_dEta_Real_sqr[i_cent][ipt_bin_t]->Write();
		h_Cnt_Bbc_dPhi_dEta_Mix [i_cent][ipt_bin_t]->Write();
		h_Cnt_Bbc_dPhi_dEta_Mix_sqr [i_cent][ipt_bin_t]->Write();
		/*for(int ipt_bin_a=0; ipt_bin_a<NPTBIN; ipt_bin_a++)
		{
			h_Cnt_Cnt_dPhi_dEta_Real[i_cent][ipt_bin_t][ipt_bin_a]->Write();
			h_Cnt_Cnt_dPhi_dEta_Mix [i_cent][ipt_bin_t][ipt_bin_a]->Write();
		}*/
	}
	for(int i_arm=0; i_arm<NARM; i_arm++)
	{
		h_Cnt_Bbc_dPhi_Real[i_cent][i_arm]->Write();
		h_Cnt_Bbc_dPhi_Mix [i_cent][i_arm]->Write();
	}

	for(int i_zvtx=0; i_zvtx<NZANA; i_zvtx++)
	{
		h_Fvt_S_N_dPhi_Real [i_cent][i_zvtx]->Write();
		h_Fvt_S_N_dPhi_dEta_Real[i_cent][i_zvtx]->Write();
		h_Fvt_S_N_dPhi_Mix  [i_cent][i_zvtx]->Write();
		h_Fvt_S_N_dPhi_dEta_Mix [i_cent][i_zvtx]->Write();
		h_Bbc_Fvt_S_dPhi_dEta_Real[i_cent][i_zvtx]->Write();
		h_Bbc_Fvt_S_dPhi_dEta_Mix[i_cent][i_zvtx]->Write();
		h_Bbc_Fvt_N_dPhi_dEta_Real[i_cent][i_zvtx]->Write();
		h_Bbc_Fvt_N_dPhi_dEta_Mix[i_cent][i_zvtx]->Write();
		h_BbcN_FvtS_dPhi_dEta_Real[i_cent][i_zvtx]->Write();
		h_BbcN_FvtS_dPhi_dEta_Mix[i_cent][i_zvtx]->Write();
		h_BbcS_FvtN_dPhi_dEta_Real[i_cent][i_zvtx]->Write();
		h_BbcS_FvtN_dPhi_dEta_Mix[i_cent][i_zvtx]->Write();

		for(int i_arm=0; i_arm<NARM; i_arm++)
		{
			h_Cnt_Fvt_dPhi_Real[i_cent][i_zvtx][i_arm]->Write();
			h_Cnt_Fvt_dPhi_Mix [i_cent][i_zvtx][i_arm]->Write();
		}
		for(int ipt_bin=0;ipt_bin<NPTBIN;ipt_bin++)
		{
			h_Cnt_Fvt_dPhi_dEta_Real[i_cent][i_zvtx][ipt_bin]->Write();
			h_Cnt_Fvt_dPhi_dEta_Mix [i_cent][i_zvtx][ipt_bin]->Write();
		}
	}//i_zvtx
  }//i_cent
}


//===================================================================================================
void book_histograms_QA	()
{
  char hname[300];
  for(int ii=0; ii<8; ii++){
	sprintf(hname, "QA_hTrk_eta%d", ii);
	h_Trk_eta[ii][0]= new TH1D(hname, hname, 250, -5, 5);

	sprintf(hname, "QA_hTrk_eta%d_cutted", ii);
	h_Trk_eta[ii][1]= new TH1D(hname, hname, 250, -5, 5);
  }
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	sprintf(hname, "QA_hCnt_pt_C%d", i_cent);
	h_Cnt_pt[i_cent]= new TH1D(hname, hname, 500, 0, 5);

	sprintf(hname, "QA_hFvt_trk_C%d", i_cent);
	h_Fvt_trk[i_cent] = new TH1D(hname, hname, 200, 0, 200);

	sprintf(hname, "QA_hBbc_cha_C%d", i_cent);
	h_Bbc_cha[i_cent] = new TH1D(hname, hname, 1000, 0, 1000);

	sprintf(hname, "QA_hFvtBbc_C%d", i_cent);
	h_Fvt_Bbc[i_cent] = new TH2D(hname, hname, 200, 0, 200, 1000, 0, 1000);

	for(int i_arm=0; i_arm<2; i_arm++)
	{
	  sprintf(hname, "QA_hFvt_dca_C%d_A%d", i_cent, i_arm);
	  h_Fvt_dca[i_cent][i_arm] = new TH1D(hname, hname, 200, -20, 20);

	  sprintf(hname, "QA_hFvt_eta_dca_C%d_A%d", i_cent, i_arm);
	  h_Fvt_eta_dca[i_cent][i_arm] = new TH2D(hname, hname, 200, -20, 20, 700, -3.5, 3.5);
	}
  }
  for(int i_arm=0; i_arm<NARM; i_arm++)
  {
	sprintf(hname, "QA_hBbc_index_A%d", i_arm);
	h_Bbc_index[i_arm] = new TH2D(hname, hname, 300,-150,150,300,-150,150);
  }
  /*for(int i_pmt=0; i_pmt<128; i_pmt++)
  {
	sprintf(hname, "QA_hBBC_dist_pmt%d", i_pmt);
	h_Bbc_cha_dist[i_pmt]=new TH1D(hname, hname, 120, -1,5);
  } // */


}

void write_histograms_QA()
{
  for(int ii=0; ii<8; ii++)	
  {
	h_Trk_eta[ii][0]->Write();
	h_Trk_eta[ii][1]->Write();
  }
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	h_Cnt_pt[i_cent]->Write();
	h_Fvt_trk[i_cent]->Write();
	h_Bbc_cha[i_cent]->Write();
	h_Fvt_Bbc[i_cent]->Write();
	for(int i_arm=0; i_arm<2; i_arm++)
	{
	  h_Fvt_dca[i_cent][i_arm]->Write();
	  h_Fvt_eta_dca[i_cent][i_arm]->Write();
	}
  }
  /*for(int i_pmt=0; i_pmt<128; i_pmt++)
  {
	h_Bbc_cha_dist[i_pmt]->Write();
  }// */  
}

//==============================================
double getBBCX(float bbc_z){
	double x_corrected = fx->Eval(bbc_z);
	return x_corrected;
}
double getBBCY(float bbc_z){
	double y_corrected = fy->Eval(bbc_z);
	return y_corrected;
}

void GetBbcIndex(){
  int bbc_index;
  float bbc_x, bbc_y;

  ifstream fbbc;
  fbbc.open("bbc_map.txt");
  while ( fbbc >> bbc_index >> bbc_x >> bbc_y ){
	cout << "READ: " << bbc_index << ", " << bbc_x << ", " << bbc_y << endl;
	if ( bbc_index<64 ){
	  //SOUTH
	  h_Bbc_index[0]->Fill(bbc_x, bbc_y, bbc_index);}
	else{
	  //NORTH
	  h_Bbc_index[1]->Fill(bbc_x, bbc_y, bbc_index);}
  }
}

int SetBbcIndex(int Bbc_arm_, float Bbc_x_, float Bbc_y_){
  int bbc_index=-999;
  bbc_index = h_Bbc_index[Bbc_arm_]->GetBinContent(h_Bbc_index[0]->FindBin(Bbc_x_,Bbc_y_));

  return bbc_index;

}
