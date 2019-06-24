#include <stdio.h>
#include <cmath>
#include "def_Const.h"

int v2Fit()
{
ifstream IF[3];
ofstream OF;
char IFNAME[4][300];
char OUTPUT[300];
char OFNAME[3][300];
char HNAME [3][300];
const float BBCN_eta=3.5, BBCS_eta=-3.5, FVTS_eta=-(1.0+1.7)/2, FVTN_eta=2.0; //eta

float eta=0;
float x_cent[NCENTANA], x_cent_err[NCENTANA], x[NPTBIN], x_err[NPTBIN];	
float v1[3][NCENTANA]={0}, v1_err[3][NCENTANA]={0}, v2[3][NCENTANA]={0}, v2_err[3][NCENTANA]={0};

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

char SUB[3][300]; 

for(int ii=0; ii<4; ii++)
{
  if(ii%4==0) 
  {
	sprintf(SUB[0], "BBCN");
	sprintf(SUB[1], "FVTN");
	sprintf(SUB[2], "FVTS");
  }
  else if(ii%4==1)
  {
	sprintf(SUB[0], "BBCN");
	sprintf(SUB[1], "FVTN");
	sprintf(SUB[2], "BBCS");
  }
  else if(ii%4==2) 
  {
	sprintf(SUB[0], "BBCN");
	sprintf(SUB[1], "FVTS");
	sprintf(SUB[2], "BBCS");
  }
  else if(ii%4==3)
  {
	sprintf(SUB[0], "FVTN");
	sprintf(SUB[1], "FVTS");
	sprintf(SUB[2], "BBCS");
  }

  for(int ref=0; ref<3; ref++)
  {
	for(int combi=0; combi<3; combi++)	//combination loop
	{
	  cout << SUB[0] << "_" << SUB[1] << "_" <<SUB[2] << "\t" << SUB[combi] << "\t" << ref << endl;

	  if(ref==0)
	  {
		sprintf(IFNAME[0], "savePar/parameter_%s_%s.txt", SUB[1], SUB[2]); //
		sprintf(IFNAME[1], "savePar/parameter_%s_%s.txt", SUB[0], SUB[2]); //
		sprintf(IFNAME[2], "savePar/parameter_%s_%s.txt", SUB[0], SUB[1]); //

		sprintf(OUTPUT, "Results/%s%s%s_%s.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  else if(ref==1)
	  {
		sprintf(IFNAME[0], "savePar/parameter_%s_%s_ref.txt", SUB[1], SUB[2]); //
		sprintf(IFNAME[1], "savePar/parameter_%s_%s_ref.txt", SUB[0], SUB[2]); //
		sprintf(IFNAME[2], "savePar/parameter_%s_%s_ref.txt", SUB[0], SUB[1]); //

		sprintf(OUTPUT, "Results/%s%s%s_%s_ref.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  else if(ref==2)
	  {
		sprintf(IFNAME[0], "savePar/parameter_%s_%s_temp.txt", SUB[1], SUB[2]); //
		sprintf(IFNAME[1], "savePar/parameter_%s_%s_temp.txt", SUB[0], SUB[2]); //
		sprintf(IFNAME[2], "savePar/parameter_%s_%s_temp.txt", SUB[0], SUB[1]); //

		sprintf(OUTPUT, "Results/%s%s%s_%s_temp.txt", SUB[0], SUB[1], SUB[2], SUB[combi]);
	  }

	  if(!strcmp(SUB[combi], "CNT") ) eta=0;
	  else if(!strcmp(SUB[combi], "FVTN") ) eta=FVTN_eta;
	  else if(!strcmp(SUB[combi], "FVTS") ) eta=FVTS_eta;
	  else if(!strcmp(SUB[combi], "BBCN") ) eta=BBCN_eta;
	  else if(!strcmp(SUB[combi], "BBCS") ) eta=BBCS_eta;

	  int i_cent;
	  float tmp_c1=0, tmp_c2=0, tmp_c1_err=0, tmp_c2_err=0, tmp;
	  float c1[3][NCENTANA]={0}, c2[3][NCENTANA]={0}, c1_err[3][NCENTANA]={0}, c2_err[3][NCENTANA]={0};

	  IF[0].open(IFNAME[0], ios::in);
	  while(IF[0] >> i_cent >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		if(tmp_c2 < 0)	c2[0][i_cent]=sqrt(-1.0);
		c1[0][i_cent] = tmp_c1;
		c2[0][i_cent] = tmp_c2;
		c1_err[0][i_cent] = tmp_c1_err;
		c2_err[0][i_cent] = tmp_c2_err;
	  }
	  IF[0].close();

	  IF[1].open(IFNAME[1], ios::in);
	  while(IF[1] >> i_cent >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		if(tmp_c2 < 0)	c2[1][i_cent]=sqrt(-1.0);
		c1[1][i_cent] = tmp_c1;
		c2[1][i_cent] = tmp_c2;
		c1_err[1][i_cent] = tmp_c1_err;
		c2_err[1][i_cent] = tmp_c2_err;
	  }
	  IF[1].close();

	  IF[2].open(IFNAME[2], ios::in);
	  while(IF[2] >> i_cent >> tmp_c1 >> tmp_c1_err >> tmp_c2 >> tmp_c2_err >> tmp >> tmp >> tmp >> tmp)
	  {
		if(tmp_c2 < 0)	c2[2][i_cent]=sqrt(-1.0);
		c1[2][i_cent] = tmp_c1;
		c2[2][i_cent] = tmp_c2;
		c1_err[2][i_cent] = tmp_c1_err;
		c2_err[2][i_cent] = tmp_c2_err;
	  }
	  IF[2].close();

	  OF.open(OUTPUT, ios::out);	  
	  for(i_cent=0; i_cent<NCENTANA; i_cent++)
	  {
		float y20 = pow(c2_err[0][i_cent]/c2[0][i_cent],2.0);	//
		float y21 = pow(c2_err[1][i_cent]/c2[1][i_cent],2.0);	//
		float y22 = pow(c2_err[2][i_cent]/c2[2][i_cent],2.0);

		float y10 = pow(c1_err[0][i_cent]/c1[0][i_cent],2.0);	//
		float y11 = pow(c1_err[1][i_cent]/c1[1][i_cent],2.0);	//
		float y12 = pow(c1_err[2][i_cent]/c1[2][i_cent],2.0);

		if(combi!=1)
		{
		  v1[combi][i_cent] = sqrt(-c1[abs(2-combi)][i_cent]*c1[abs(1-combi)][i_cent]/c1[combi][i_cent] );
		  v2[combi][i_cent] = sqrt( c2[abs(2-combi)][i_cent]*c2[abs(1-combi)][i_cent]/c2[combi][i_cent] );
		}
		else if(combi==1)
		{
		  v1[combi][i_cent] = sqrt(-c1[abs(1+combi)][i_cent]*c1[abs(1-combi)][i_cent]/c1[combi][i_cent]);
		  v2[combi][i_cent] = sqrt( c2[abs(1+combi)][i_cent]*c2[abs(1-combi)][i_cent]/c2[combi][i_cent]);
		}
		v1_err[combi][i_cent] = (0.5)*v1[combi][i_cent]*sqrt(y10+y11+y12);
		v2_err[combi][i_cent] = (0.5)*v2[combi][i_cent]*sqrt(y20+y21+y22);

		cout << i_cent << "\t" << v1[combi][i_cent] << "\t" << v1_err[combi][i_cent] << "\t" << v2[combi][i_cent] << "\t" << v2_err[combi][i_cent] << endl;
		OF	 << i_cent << "\t" << "0" << "\t" << v1[combi][i_cent] << "\t" << v1_err[combi][i_cent] << "\t" << v2[combi][i_cent] << "\t" << v2_err[combi][i_cent] << "\t" << eta << endl;
	  }
	  OF.close();

	  TGraph *gv1[3];
	  TGraph *gv2[3];
	  for(int i_combi=0; i_combi<3; i_combi++)
	  {
		gv1[i_combi] = new TGraphErrors(NCENTANA, &x_cent[0], &v1[i_combi][0], &x_cent_err[0], &v1_err[i_combi][0]);
		gv1[i_combi] ->SetMarkerStyle(20);
		gv1[i_combi] ->SetMarkerSize(1.1);
		gv1[i_combi] ->SetMarkerColor(i_combi+1);
		gv1[i_combi] ->SetLineColor(i_combi+1);
		gv1[i_combi] ->SetLineStyle(2);
		gv1[i_combi] ->SetLineWidth(3);

		gv2[i_combi] = new TGraphErrors(NCENTANA, &x_cent[0], &v2[i_combi][0], &x_cent_err[0], &v2_err[i_combi][0]);
		gv2[i_combi] ->SetMarkerStyle(20);
		gv2[i_combi] ->SetMarkerSize(1.1);
		gv2[i_combi] ->SetMarkerColor(i_combi+1);
		gv2[i_combi] ->SetLineColor(i_combi+1);
		gv2[i_combi] ->SetLineStyle(2);
		gv2[i_combi] ->SetLineWidth(3);
	  }

	  TH1F *htmp;
	  TCanvas *c_1 = new TCanvas("c_1","v_{1} fitting", 500, 370);
	  TLegend *leg = new TLegend(0.15, 0.6, 0.48, 0.9);
	  leg->SetFillStyle(0);
      leg->SetBorderSize(0);
	  leg->SetTextSize(0.05);

	  gPad->SetMargin(0.13, 0.03, 0.13, 0.05);
	  htmp = (TH1F*)gPad->DrawFrame(0, 0.02, 84, 0.2); //xmin, ymin, xmax, ymax

	  htmp->GetXaxis()->SetTitle("Centrality");
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


	  gv1[combi]->Draw("pe, same");
	  leg->AddEntry(gv1[combi],SUB[combi],"LP");

	  leg->Draw();

	  TCanvas *c_2 = new TCanvas("c_2","v_{2} fitting", 500, 370);
	  gPad->SetMargin(0.13, 0.03, 0.13, 0.05);
	  htmp = (TH1F*)gPad->DrawFrame(0, 0.02, 84, 0.2); //xmin, ymin, xmax, ymax

	  htmp->GetXaxis()->SetTitle("Centrality");
	  htmp->GetXaxis()->SetTitleSize(0.06);
	  htmp->GetXaxis()->SetTitleFont(62);
	  htmp->GetXaxis()->SetTitleOffset(0.9);
	  htmp->GetXaxis()->CenterTitle(kTRUE);
	  htmp->GetXaxis()->SetLabelSize(0.05);
	  htmp->GetXaxis()->SetLabelFont(62);

	  htmp->GetYaxis()->SetTitle("v_{2}");
	  htmp->GetYaxis()->SetTitleSize(0.07);
	  htmp->GetYaxis()->SetTitleFont(62);
	  htmp->GetYaxis()->SetTitleOffset(0.9);
	  htmp->GetYaxis()->CenterTitle(kTRUE);
	  htmp->GetYaxis()->SetLabelSize(0.05);
	  htmp->GetYaxis()->SetLabelFont(62);
	  htmp->SetMinimum(0);
	  htmp->SetMaximum(.2);
	  
	  gv2[combi]->Draw("pe, same");
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
	}//combi
  }//ref




}// end of ii
return 0;
} // end of function
