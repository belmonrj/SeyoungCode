#include <stdio.h>
#include <cmath>
#include "def_Const.h"

void v2Fit_FVTX()
{
ifstream IF[3], IFX;
ofstream OF;
char IFNAME[4][300];
char OUTPUT[300];
const float BBCN_eta=3.5, BBCS_eta=-3.5, FVTS_eta=-(1.0+1.7)/2, FVTN_eta=2.0; //eta
float eta=0;
float x_cent[NCENTANA], x_cent_err[NCENTANA], x[NPTBIN], x_err[NPTBIN];	

//Set x_cent
for(int i=0; i<NCENTANA; i++)
{
  if(i==0) {
	x_cent[i] = cent[i]/2;
	x_cent_err[i] = cent[i]/2;
  }
  else {
	x_cent[i] = (cent[i-1]+cent[i])/2;
	x_cent_err[i] = (cent[i]-cent[i-1])/2;
  }
}

//Read x_pT
int pt=0;
sprintf(IFNAME[3], "savePar/pT_x_err.txt");
IFX.open(IFNAME[3], ios::in);
while(IFX >> pt >> x[pt] >> x_err[pt]) {	cout << x[pt] << "\t" << x_err[pt] << endl;}
IFX.close();

char SUB[3][300]; 

for(int ii=0; ii<6; ii++)
{
  if(ii%6==0) 
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "FVTN");
	sprintf(SUB[2], "FVTS");
  }
  else if(ii%6==1)
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "BBCN");
	sprintf(SUB[2], "BBCS");
  }
  else if(ii%6==2) 
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "FVTS");
	sprintf(SUB[2], "BBCS");
  }
  else if(ii%6==3)
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "BBCN");
	sprintf(SUB[2], "FVTN");
  }
  else if(ii%6==4)
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "BBCN");
	sprintf(SUB[2], "FVTS");
  }

  else if(ii%6==5)
  {
	sprintf(SUB[0], "CNT");
	sprintf(SUB[1], "FVTN");
	sprintf(SUB[2], "BBCS");
  }

  for(int combi=0; combi<3; combi++)
  {
	for(int ref=0; ref<3; ref++) //reference fitting
	{
	  cout << SUB[0] << "_" << SUB[1] << "_" <<SUB[2] << "\t" << SUB[combi] << "\t" << ref << endl;

	  if(ref==0)
	  {
		sprintf(IFNAME[2], "savePar/parameter_%s_%s.txt", SUB[0], SUB[1]); //CNT-North
		sprintf(IFNAME[0], "savePar/parameter_%s_%s.txt", SUB[1], SUB[2]); //North-South
		sprintf(IFNAME[1], "savePar/parameter_%s_%s.txt", SUB[0], SUB[2]); //CNT-South

		sprintf(OUTPUT, "Results/%s%s%s_%s.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  else if(ref==1)
	  {
		sprintf(IFNAME[2], "savePar/parameter_%s_%s_ref.txt", SUB[0], SUB[1]);
		sprintf(IFNAME[0], "savePar/parameter_%s_%s_ref.txt", SUB[1], SUB[2]);
		sprintf(IFNAME[1], "savePar/parameter_%s_%s_ref.txt", SUB[0], SUB[2]);

		sprintf(OUTPUT, "Results/%s%s%s_%s_ref.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  else if(ref==2)
	  {
		sprintf(IFNAME[2], "savePar/parameter_%s_%s_temp.txt", SUB[0], SUB[1]);
		sprintf(IFNAME[0], "savePar/parameter_%s_%s_temp.txt", SUB[1], SUB[2]);
		sprintf(IFNAME[1], "savePar/parameter_%s_%s_temp.txt", SUB[0], SUB[2]);

		sprintf(OUTPUT, "Results/%s%s%s_%s_temp.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  if(!strcmp(SUB[combi], "CNT") ) eta=0;
	  else if(!strcmp(SUB[combi], "FVTN") ) eta=FVTN_eta;
	  else if(!strcmp(SUB[combi], "FVTS") ) eta=FVTS_eta;
	  else if(!strcmp(SUB[combi], "BBCN") ) eta=BBCN_eta;
	  else if(!strcmp(SUB[combi], "BBCS") ) eta=BBCS_eta;

	  int i_cent, i_pt;
	  float tmp_c1=0, tmp_c2=0, tmp_c1_err=0, tmp_c2_err=0, tmp;
	  float c1[3][NCENTANA][NPTBIN]={0}, c2[3][NCENTANA][NPTBIN]={0}, c1_err[3][NCENTANA][NPTBIN]={0}, c2_err[3][NCENTANA][NPTBIN]={0};

	  IF[1].open(IFNAME[1], ios::in);
	  while(IF[1] >> i_cent >> i_pt >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		if(tmp_c2 < 0)	c2[1][i_cent][i_pt]=sqrt(-1.0);
		c1[1][i_cent][i_pt] = tmp_c1;
		c2[1][i_cent][i_pt] = tmp_c2;
		c1_err[1][i_cent][i_pt] = tmp_c1_err;
		c2_err[1][i_cent][i_pt] = tmp_c2_err;
	  }
	  IF[1].close();

	  IF[2].open(IFNAME[2], ios::in);
	  while(IF[2] >> i_cent >> i_pt >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		if(tmp_c2 < 0)	c2[2][i_cent][i_pt]=sqrt(-1.0);
		c1[2][i_cent][i_pt] = tmp_c1;
		c2[2][i_cent][i_pt] = tmp_c2;
		c1_err[2][i_cent][i_pt] = tmp_c1_err;
		c2_err[2][i_cent][i_pt] = tmp_c2_err;
	  }
	  IF[2].close();

	  IF[0].open(IFNAME[0], ios::in);
	  while(IF[0] >> i_cent >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		for(int tmp_i=0; tmp_i<NPTBIN; tmp_i++)
		{
		  if(tmp_c2 < 0)	c2[0][i_cent][tmp_i]=sqrt(-1.0);
		  c1[0][i_cent][tmp_i] = tmp_c1;
		  c2[0][i_cent][tmp_i] = tmp_c2;
		  c1_err[0][i_cent][tmp_i] = tmp_c1_err;
		  c2_err[0][i_cent][tmp_i] = tmp_c2_err;
		}
	  }
	  IF[0].close();

	  float v1[3][NCENTANA][NPTBIN], v1_err[3][NCENTANA][NPTBIN], v2[3][NCENTANA][NPTBIN], v2_err[3][NCENTANA][NPTBIN];

	  OF.open(OUTPUT, ios::out);	  
	  for(i_cent=0; i_cent<NCENTANA; i_cent++)
	  {
		for(i_pt=0; i_pt<NPTBIN; i_pt++)
		{
		  float y20 = pow(c2_err[1][i_cent][i_pt]/c2[1][i_cent][i_pt],2.0);	//South
		  float y21 = pow(c2_err[2][i_cent][i_pt]/c2[2][i_cent][i_pt],2.0);	//North
		  float y22 = pow(c2_err[0][i_cent][0]/c2[0][i_cent][0],2.0);

		  float y10 = pow(c1_err[1][i_cent][i_pt]/c1[1][i_cent][i_pt],2.0);	//South
		  float y11 = pow(c1_err[2][i_cent][i_pt]/c1[2][i_cent][i_pt],2.0);	//North
		  float y12 = pow(c1_err[0][i_cent][0]/c1[0][i_cent][0],2.0);

		  if(combi!=1)
		  {
			v1[combi][i_cent][i_pt] = sqrt(-c1[abs(2-combi)][i_cent][i_pt]*c1[abs(1-combi)][i_cent][i_pt]/c1[combi][i_cent][i_pt] );
			v2[combi][i_cent][i_pt] = sqrt( c2[abs(2-combi)][i_cent][i_pt]*c2[abs(1-combi)][i_cent][i_pt]/c2[combi][i_cent][i_pt] );
		  }
		  else if(combi==1)
		  {
			v1[combi][i_cent][i_pt] = sqrt(-c1[abs(1+combi)][i_cent][i_pt]*c1[abs(1-combi)][i_cent][i_pt]/c1[combi][i_cent][i_pt] );
			v2[combi][i_cent][i_pt] = sqrt( c2[abs(1+combi)][i_cent][i_pt]*c2[abs(1-combi)][i_cent][i_pt]/c2[combi][i_cent][i_pt] );
		  }
		  v1_err[combi][i_cent][i_pt] = (0.5)*v1[combi][i_cent][i_pt]*sqrt(y10+y11+y12);
		  v2_err[combi][i_cent][i_pt] = (0.5)*v2[combi][i_cent][i_pt]*sqrt(y20+y21+y22);

		  cout << i_cent << "\t" << i_pt << "\t" << v1[combi][i_cent][i_pt] << "\t" << v1_err[combi][i_cent][i_pt] << "\t" << v2[combi][i_cent][i_pt] << "\t" << v2_err[combi][i_cent][i_pt] << endl;
		  OF << i_cent << "\t" << i_pt << "\t" << v1[combi][i_cent][i_pt] << "\t" << v1_err[combi][i_cent][i_pt] << "\t" << v2[combi][i_cent][i_pt] << "\t" << v2_err[combi][i_cent][i_pt] << "\t" << eta << endl;
		}
	  }
	  OF.close();

	  TGraph *gv1[NCENTANA];
	  TGraph *gv2[NCENTANA];

	  for(i_cent=0; i_cent<NCENTANA; i_cent++)
	  {
		gv1[i_cent] = new TGraphErrors(NPTBIN, &x[0], &v1[combi][i_cent][0], &x_err[0], &v1_err[combi][i_cent][0]);
		gv1[i_cent] ->SetMarkerStyle(20);
		gv1[i_cent] ->SetMarkerSize(1.1);
		gv1[i_cent] ->SetMarkerColor(i_cent+1);
		gv1[i_cent] ->SetLineColor(i_cent+1);
		gv1[i_cent] ->SetLineStyle(2);
		gv1[i_cent] ->SetLineWidth(3);

		gv2[i_cent] = new TGraphErrors(NPTBIN, &x[0], &v2[combi][i_cent][0], &x_err[0], &v2_err[combi][i_cent][0]);
		gv2[i_cent] ->SetMarkerStyle(20);
		gv2[i_cent] ->SetMarkerSize(1.1);
		gv2[i_cent] ->SetMarkerColor(i_cent+1);
		gv2[i_cent] ->SetLineColor(i_cent+1);
		gv2[i_cent] ->SetLineStyle(2);
		gv2[i_cent] ->SetLineWidth(3);

		if(i_cent>=4)
		{
			gv1[i_cent] ->SetMarkerColor(i_cent+2);
			gv1[i_cent] ->SetLineColor(i_cent+2);
			gv2[i_cent] ->SetMarkerColor(i_cent+2);
			gv2[i_cent] ->SetLineColor(i_cent+2);
		}
	  }


	  TH1F *htmp;
	  TCanvas *c_1 = new TCanvas("c_1","v_{1} fitting", 500, 370);
	  TLegend *leg = new TLegend(0.15, 0.6, 0.48, 0.9);
	  leg->SetFillStyle(0);
      leg->SetBorderSize(0);
	  leg->SetTextSize(0.05);

	  gPad->SetMargin(0.13, 0.03, 0.13, 0.05);
	  htmp = (TH1F*)gPad->DrawFrame(0, 0.02, 4.5, 0.2); //xmin, ymin, xmax, ymax

	  htmp->GetXaxis()->SetTitle("p_{T}^{CNT} (GeV/c)");
	  htmp->GetXaxis()->SetTitleSize(0.06);
	  htmp->GetXaxis()->SetTitleFont(62);
	  htmp->GetXaxis()->SetTitleOffset(0.9);
	  htmp->GetXaxis()->CenterTitle(kTRUE);
	  htmp->GetXaxis()->SetLabelSize(0.05);
	  htmp->GetXaxis()->SetLabelFont(62);

	  htmp->GetYaxis()->SetTitle("v_{1}");
	  htmp->GetYaxis()->SetTitleSize(0.07);
	  htmp->GetYaxis()->SetTitleFont(62);
	  htmp->GetYaxis()->SetTitleOffset(0.9);
	  htmp->GetYaxis()->CenterTitle(kTRUE);
	  htmp->GetYaxis()->SetLabelSize(0.05);
	  htmp->GetYaxis()->SetLabelFont(62);
	  htmp->SetMinimum(0);
	  htmp->SetMaximum(1.);


	  gv1[0]->Draw("pe");
	  if(NCENTANA>1)
	  {
		for(int i_cent=1; i_cent<NCENTANA; i_cent++) gv1[i_cent]->Draw("same, pe");
	  }

	  leg->AddEntry(gv1[0],"p+Au BBCS central 0-1%","LP");
	  if(NCENTANA>1)
	  {
		leg->AddEntry(gv1[1],"p+Au BBCS MinBias 1-3%","LP");
		if(NCENTANA>2){
		  leg->AddEntry(gv1[2],"p+Au BBCS MinBias 3-5%","LP");
		  if(NCENTANA>3){
			leg->AddEntry(gv1[3],"p+Au BBCS MinBias 5-10%","LP");
			if(NCENTANA>4){
			  leg->AddEntry(gv1[4],"p+Au BBCS MinBias 10-20%","LP");
			  if(NCENTANA>5){
				leg->AddEntry(gv1[5],"p+Au BBCS MinBias 20-60%","LP");
				if(NCENTANA>6){
				  leg->AddEntry(gv1[6],"p+Au BBCS MinBias 60-84%","LP");} }}}}
	  } // end of if NCENTANA 
	  leg->Draw();

	  TCanvas *c_2 = new TCanvas("c_2","v_{2} fitting", 500, 370);
	  gPad->SetMargin(0.13, 0.03, 0.13, 0.05);
	  htmp = (TH1F*)gPad->DrawFrame(0, 0.02, 4.5, 0.2); //xmin, ymin, xmax, ymax

	  htmp->GetXaxis()->SetTitle("p_{T}^{CNT} (GeV/c)");
	  htmp->GetXaxis()->SetTitleSize(0.06);
	  htmp->GetXaxis()->SetTitleFont(62);
	  htmp->GetXaxis()->SetTitleOffset(0.9);
	  htmp->GetXaxis()->CenterTitle(kTRUE);
	  htmp->GetXaxis()->SetLabelSize(0.05);
	  htmp->GetXaxis()->SetLabelFont(62);

	  if(combi==0)		htmp->GetYaxis()->SetTitle("v_{2}^{CNT}");
	  else if(combi==1)	htmp->GetYaxis()->SetTitle("v_{2}^{North}");
	  else				htmp->GetYaxis()->SetTitle("v_{2}^{South}");

	  htmp->GetYaxis()->SetTitleSize(0.07);
	  htmp->GetYaxis()->SetTitleFont(62);
	  htmp->GetYaxis()->SetTitleOffset(0.9);
	  htmp->GetYaxis()->CenterTitle(kTRUE);
	  htmp->GetYaxis()->SetLabelSize(0.05);
	  htmp->GetYaxis()->SetLabelFont(62);
	  htmp->SetMinimum(0);
	  if(ref==0) htmp->SetMaximum(.8);
	  else htmp->SetMaximum(0.1);
	  if(combi!=0) htmp->SetMaximum(.1);
	  
	  gv2[0]->Draw("lpe");
	  if(NCENTANA>1)
	  {
		for(int i_cent=1; i_cent<NCENTANA; i_cent++)
		  gv2[i_cent]->Draw("same, lpe");
	  }

	  leg->SetFillStyle(0);
      leg->SetBorderSize(0);

      leg->Draw();
	  
	  char outpdf1[300];
	  char outpdf2[300];

	  if(ref==0)
	  {
		sprintf(outpdf1, "Results/v1_%s%s%s_%s.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
		sprintf(outpdf2, "Results/v2_%s%s%s_%s.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }
	  else if(ref==1)
	  {
		sprintf(outpdf1, "Results/v1_%s%s%s_%s_ref.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
		sprintf(outpdf2, "Results/v2_%s%s%s_%s_ref.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }
	  else if(ref==2)
	  {
		sprintf(outpdf1, "Results/v1_%s%s%s_%s_temp.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
		sprintf(outpdf2, "Results/v2_%s%s%s_%s_temp.pdf", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  c_1->SaveAs(outpdf1);
	  c_2->SaveAs(outpdf2);
	  delete c_1;
	  delete c_2;
	}//ref
  }//combi
}//end of loop

}
