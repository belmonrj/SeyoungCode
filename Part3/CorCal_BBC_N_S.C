#include <stdio.h>
#include <cmath>
#include <iostream>

#include "sys/types.h"
#include "sys/stat.h"
#include "unistd.h"
#include "def_Const.h"

TH1F *href;
TH1F *href2;
const float PI=acos(-1.);

double fit_func(double *x, double *par)
{
  double xx = x[0];
  double yy = par[0]*(href->GetBinContent(href->FindBin(xx))) + par[1];
  return yy;
}

int CorCal_BBC_N_S()
{
  char SUB[16]="BBCN_BBCS";

  gStyle->SetOptStat(0);
  char ref_name[300];
  char SaveHisto[300];
  char output_name[300];
  sprintf(ref_name, "%s/CorHisto_%d_%s.root", ref_dir, ref_run, SUB);

  double c1[NCENTANA], c2[NCENTANA], c3[NCENTANA], c4[NCENTANA];
  double c1err[NCENTANA], c2err[NCENTANA], c3err[NCENTANA], c4err[NCENTANA];
  double c1tmp[NCENTANA], c2tmp[NCENTANA];
  double c1tmp_err[NCENTANA], c2tmp_err[NCENTANA];
  double c1ref[NCENTANA], c2ref[NCENTANA];
  double c1ref_err[NCENTANA], c2ref_err[NCENTANA];

  TFile *rootopen = TFile::Open(input_name, "READ");
  TFile *ref_open = TFile::Open(ref_name, "READ");

  cout << "Open file : " << input_name << endl;
  cout << "Open file : " << ref_name << endl;

  TCanvas *C_1  = new TCanvas("C_1", "C_1", 500, 400);
  TLine	*l_1 = new TLine(-0.5*PI, 1, 1.5*PI, 1);
  TCanvas *ComCanv;

  //Fitting function for correlation histogram
  TF1		*fit_com[15];
  TF1		*f1_com[15];
  TF1		*f2_com[15];
  TF1		*f3_com[15];
  TF1		*f4_com[15];
  
  for(int i=0; i<15; i++)
  {
	fit_com[i] = new TF1("fit", "1 + (2*[0]*cos(x)) + (2*[1]*cos(2*x)) + (2*[2]*cos(3*x)) + (2*[3]*cos(4*x))", -0.5*PI, 1.5*PI);
	f1_com[i] = new TF1("f1" , "1 + (2*[0]*cos(x))"  , -0.5*PI, 1.5*PI);
	f2_com[i] = new TF1("f2" , "1 + (2*[0]*cos(2*x))", -0.5*PI, 1.5*PI);
	f3_com[i] = new TF1("f3" , "1 + (2*[0]*cos(3*x))", -0.5*PI, 1.5*PI);
	f4_com[i] = new TF1("f4" , "1 + (2*[0]*cos(4*x))", -0.5*PI, 1.5*PI);

	fit_com[i]->SetLineColor(kAzure);
	f1_com[i]->SetLineColor(kMagenta);
	f1_com[i]->SetLineStyle(2);
	f2_com[i]->SetLineColor(kPink);
	f2_com[i]->SetLineStyle(2);
	f3_com[i]->SetLineColor(kOrange);
	f3_com[i]->SetLineStyle(2);
	f4_com[i]->SetLineColor(kTeal);
	f4_com[i]->SetLineStyle(2);
  }


  //Fitting function for template fitting
  href			= new TH1F("Reference correlation function", "href", H_BBC_phi, -0.5*PI, 1.5*PI);
  TF1	*f2_tmp = new TF1("f2_tmp" , "1 + 2*[0]*cos(2*x)", -0.5*PI, 1.5*PI);
  TH1	*f_fit;

  f_fit= new TH1D("f_fit", "f_fit", H_BBC_phi, -0.5*PI, 1.5*PI);
  f_fit->SetMarkerStyle(24);
  f_fit->SetMarkerSize(1.4);

  f2_tmp->SetLineColor(kRed);
  f2_tmp->SetLineStyle(7);

  struct stat st0 = {0};
  if(stat("ALL", &st0) == -1)  mkdir("ALL", 0700);
  struct stat st1 = {0};
  if(stat("savePar", &st1) == -1)  mkdir("savePar", 0700);
  struct stat st2 = {0};
  if(stat("ALL/combined", &st2) == -1)  mkdir("ALL/combined", 0700);
  struct stat st3 = {0};
  if(stat(SUB, &st3) == -1)  mkdir(SUB, 0700);

  ofstream SavePar, SavePar_ref, SavePar_temp;
  char param_name[300];

  TH1		*histo_original	[NCENTANA];
  TH1		*histo_original_Sub[NCENTANA];

  char hname[200];
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	if(i_cent==NCENTANA-1 && saveOption==1) sprintf(hname, "COR_BBCN_BBCS_C%d", 0);
	else sprintf(hname, "SIG_BBCN_BBCS_C%d", i_cent);
	histo_original[i_cent] = new TH1D(hname, hname, H_BBC_phi, -0.5*PI, 1.5*PI);

	sprintf(hname, "SUB_BBCN_BBCS_C%d", i_cent);
	histo_original_Sub[i_cent] = new TH1D(hname, hname, H_BBC_phi, -0.5*PI, 1.5*PI);
  }

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	histo_original[i_cent]->GetXaxis()->SetTitle("delta phi");
	histo_original[i_cent]->SetMarkerColor(kBlack);
	histo_original[i_cent]->SetLineColor(kBlack);
	histo_original[i_cent]->SetLineWidth(2);
	histo_original[i_cent]->SetMarkerStyle(20);
	histo_original[i_cent]->SetAxisRange(0.99, 1.01, "Y");

	histo_original_Sub[i_cent]->GetXaxis()->SetTitle("delta phi");
	histo_original_Sub[i_cent]->SetMarkerColor(kBlack);
	histo_original_Sub[i_cent]->SetLineColor(kBlack);
	histo_original_Sub[i_cent]->SetMarkerStyle(20);
	histo_original_Sub[i_cent]->SetAxisRange(0.9995, 1.0005, "Y");
  }

  sprintf(output_name, "CorHisto_%s.root", SUB);  
  sprintf(param_name, "savePar/parameter_%s.txt", SUB);
  TFile *outfile = new TFile(output_name, "recreate");  

  href = new TH1F("Reference correlation function", "href", H_BBC_phi, -0.5*PI, 1.5*PI);	

  SavePar.open(param_name, ios::trunc);
  cout << "Parameter save file : " << param_name << endl;

  TH2D *sumR  = new TH2D("sumR",  "Real pT sum histogram", H_BBC_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -6, 0);
  TH2D *sumM  = new TH2D("sumM",  "Mix pT sum histogram", H_BBC_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -6, 0);
  TH2D *tmpsumR  = new TH2D("tmpsumR",  " ", H_BBC_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -6, 0);
  TH2D *tmpsumM  = new TH2D("tmpsumM",  " ", H_BBC_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -6, 0);

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	  TH2D *h1 = (TH2D*)rootopen->Get(Form("Real_dPhi_dEta_BBC_S_N_C%d", i_cent));
	  TH2D *h2 = (TH2D*)rootopen->Get(Form("Mix_dPhi_dEta_BBC_S_N_C%d", i_cent));

	  for (int i_phi=0; i_phi<H_BBC_phi; i_phi++)
	  {
		float tmpaddhistoR=0, tmpaddhistoM=0, tmpaddhistoR2=0, tmpaddhistoM2=0;
		for(int ieta=0; ieta<H_FVTX_eta; ieta++)
		{
		  tmpaddhistoR  = h1->GetBinContent(i_phi+1, ieta+1);
		  tmpaddhistoM  = h2->GetBinContent(i_phi+1, ieta+1);
		  tmpsumR->SetBinContent(i_phi+1, ieta+1, tmpaddhistoR);
		  tmpsumM->SetBinContent(i_phi+1, ieta+1, tmpaddhistoM);
		} //ieta
	  }//i_phi
      sumR->Add(tmpsumR);
	  sumM->Add(tmpsumM);

	TH1D *p_R  = sumR ->ProjectionX();
	TH1D *p_M  = sumM ->ProjectionX();

	char RealHisto[300];
	char MixHisto[300];
	sprintf(RealHisto, "ALL/%s_Real_C%d.pdf", SUB, i_cent);
	sprintf(MixHisto,  "ALL/%s_Mix_C%d.pdf",  SUB, i_cent);

	sumR->Draw("colz");
	p_R->Draw("ep0");
	C_1->SaveAs(RealHisto);

	sumM->Draw("colz");
	p_M ->Draw("ep0");
	C_1->SaveAs(MixHisto);

	float sum_R = sumR->GetSum();
	float sum_M = sumM->GetSum();

	for(int i_phi=0; i_phi<H_BBC_phi; i_phi++)
	{
	  float val_R = p_R->GetBinContent(i_phi+1);
	  float val_M = p_M->GetBinContent(i_phi+1);
	  float ds_R  = p_R->GetBinError(i_phi+1);
	  float ds_M  = p_M->GetBinError(i_phi+1);
	  float sumRatio = 0;
	  float valRatio = 0;
	  float dsig = 0;

	  // Check denominator of the bin value is not zero
	  if(sum_R>0 && val_M>0)
	  {
		sumRatio = sum_M/sum_R;				 // Normalization factor
		valRatio = (val_R/val_M) * sumRatio; // bin value for corH, (val_Real/sum_Real)/(val_Mix/sum_Mix)
	  }
	  // Check denominator of the ds is not zero & Calculate propagation error "ds"
	  if(val_R!=0 && val_M!=0)
	  {
		dsig = valRatio * sqrt( (ds_R/val_R)*(ds_R/val_R) + (ds_M/val_M)*(ds_M/val_M) );
	  }
	  histo_original[i_cent]->SetBinContent(i_phi+1, valRatio);
	  histo_original[i_cent]->SetBinError(i_phi+1, dsig);
	}// i_phi End

	histo_original[i_cent]->Draw("e0");
	histo_original[i_cent]->Write();
	histo_original[i_cent]->Fit(fit_com[i_cent], "RM");
	fit_com[i_cent]->Draw("same");
	l_1 ->Draw("same");
	// Save parameter
	c1[i_cent]	=fit_com[i_cent]->GetParameter(0);
	c2[i_cent]	=fit_com[i_cent]->GetParameter(1);
	c3[i_cent]	=fit_com[i_cent]->GetParameter(2);
	c4[i_cent]	=fit_com[i_cent]->GetParameter(3);
	c1err[i_cent]	=fit_com[i_cent]->GetParError(0);
	c2err[i_cent]	=fit_com[i_cent]->GetParError(1);
	c3err[i_cent]	=fit_com[i_cent]->GetParError(2);
	float chi2[NCENTANA]={0};

	// Calculate chi^2/NDF
	for(int i_phi=0; i_phi<H_BBC_phi; i_phi++)
	{
	  double data=histo_original[i_cent]->GetBinContent(i_phi+1);
	  double tmpx=histo_original[i_cent]->GetBinCenter(i_phi+1);
	  double tmperr=histo_original[i_cent]->GetBinError(i_phi+1);
	  double fitval=fit_com[i_cent]->Eval(tmpx);
	  double tmpval=(fitval-data)*(fitval-data)/(tmperr*tmperr);
	  chi2[i_cent] = chi2[i_cent]+tmpval;
	}

	SavePar << i_cent << "\t" << c1[i_cent] << "\t" << c1err[i_cent] << "\t" << c2[i_cent] << "\t" << c2err[i_cent] << "\t" << chi2[i_cent] << "\t" << H_BBC_phi-3 << "\t" << chi2[i_cent]/(H_BBC_phi-3) << "\t" << c2[i_cent]/c1[i_cent] << "\n" ;

	f1_com[i_cent]->SetParameter(0, c1[i_cent]);
	f1_com[i_cent]->Draw("same");
	f2_com[i_cent]->SetParameter(0, c2[i_cent]);
	f2_com[i_cent]->Draw("same");
	f3_com[i_cent]->SetParameter(0, c3[i_cent]);
	f3_com[i_cent]->Draw("same");
	f4_com[i_cent]->SetParameter(0, c4[i_cent]);
	f4_com[i_cent]->Draw("same");

	sprintf(SaveHisto, "ALL/%s_Cor_C%d.pdf", SUB, i_cent);
	C_1->SaveAs(SaveHisto);
	sumR->Reset();
	sumM->Reset();
	delete p_R;
	delete p_M;
  }
  SavePar.close();

  //REFERENCE fitting
if(subtraction!=0)
{
  TH1D *REF;
  if(subtraction==1) REF = (TH1D*)ref_open->Get(Form("COR_BBCN_BBCS_C%d", 0)); //Ref
  else if(subtraction==2) REF = (TH1D*)histo_original[NCENTANA-1]; //pAu peripheral

  sprintf(param_name, "savePar/parameter_%s_ref.txt", SUB);
  SavePar_ref.open(param_name, ios::trunc);
  sprintf(param_name, "savePar/parameter_%s_temp.txt", SUB);
  SavePar_temp.open(param_name, ios::trunc);

  double PAR[3][NCENTANA];
  for(int i_phi=0; i_phi<H_BBC_phi; i_phi++)
  {
	href->SetBinContent(i_phi+1, REF->GetBinContent(i_phi+1));
	href->SetBinError(i_phi+1, REF->GetBinError(i_phi+1));
  }
  

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	TF1	*ffit = new TF1("ffit", *fit_func, -0.5*PI, 1.5*PI, 2);

	histo_original[i_cent]->Fit(ffit, "N"); //Fitting wo drawing

	PAR[0][i_cent]=ffit->GetParameter(0);
	PAR[1][i_cent]=ffit->GetParameter(1);
	//cout << PAR[0][i_cent] << "  " << PAR[1][i_cent] << endl;
	for(int i_phi=0; i_phi<H_BBC_phi; i_phi++)
	{
	  double a_tmp=PAR[0][i_cent];
	  double b_tmp=PAR[1][i_cent];
	  double c_tmp=REF->GetBinContent(i_phi+1);
	  double d_tmp=a_tmp*c_tmp + b_tmp;
	  f_fit->SetBinContent(i_phi+1, d_tmp );
	  f_fit->SetBinError(i_phi+1, 0 );
	}
	f2_tmp->SetParameter(0, (fit_com[i_cent]->GetParameter(1))-PAR[0][i_cent]*(fit_com[NCENTANA-1]->GetParameter(1)) );

	histo_original[i_cent]->Draw("ep");
	f_fit->Draw("ep, same");
	l_1->Draw("same");
	f2_tmp->Draw("same");

	sprintf(SaveHisto, "%s/GetF_%d.pdf", SUB, i_cent);
	C_1->SaveAs(SaveHisto);

  }

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	for(int i_phi=0; i_phi<H_BBC_phi; i_phi++)
	{
	  double ii = histo_original[i_cent]->GetBinContent(i_phi+1);
	  double tmp_v=REF->GetBinContent(i_phi+1);
	  double jj = (PAR[0][i_cent]*tmp_v)+ PAR[1][i_cent];
	  double err_ref=sqrt( (histo_original[i_cent]->GetBinError(i_phi+1)*histo_original[i_cent]->GetBinError(i_phi+1)) + PAR[0][i_cent]*PAR[0][i_cent]*(REF->GetBinError(i_phi+1)*REF->GetBinError(i_phi+1)) );

	  histo_original_Sub[i_cent]->SetBinContent(i_phi+1, ii-jj+1);
	  histo_original_Sub[i_cent]->SetBinError(i_phi+1, err_ref);
	  histo_original_Sub[i_cent]->SetMarkerStyle(33);
	  histo_original_Sub[i_cent]->SetLineColor(kPink);
	  histo_original_Sub[i_cent]->SetMarkerColor(kPink);
	}
	histo_original_Sub[i_cent]->Draw("p");
	l_1->Draw("same");
	histo_original_Sub[i_cent]->Fit(fit_com[i_cent], "RM");

	c1ref[i_cent]		=fit_com[i_cent]->GetParameter(0);
	c2ref[i_cent]		=fit_com[i_cent]->GetParameter(1);
	float p2			=fit_com[i_cent]->GetParameter(2);
	float p3			=fit_com[i_cent]->GetParameter(3);
	c1ref_err[i_cent]	=fit_com[i_cent]->GetParError(0);
	c2ref_err[i_cent]	=fit_com[i_cent]->GetParError(1);
	float e2			=fit_com[i_cent]->GetParError(2);
	float chi2[NCENTANA]={0};

	SavePar_ref << i_cent << "\t" << c1ref[i_cent] << "\t" << c1ref_err[i_cent]  << "\t" << c2ref[i_cent] << "\t" << c2ref_err[i_cent] << "\t" << chi2[i_cent] << "\t" <<  H_BBC_phi-3 << "\t" << chi2[i_cent]/(H_FVTX_phi-3) << "\t" << c2ref[i_cent]/c1ref[i_cent] << "\n";

	cout << c2[i_cent]-(PAR[0][i_cent]*c2[NCENTANA-1]) << "  " << c2ref[i_cent] << endl;

	SavePar_temp << i_cent << "\t" << c1ref[i_cent] << "\t" << c1ref_err[i_cent]  << "\t" << c2ref[i_cent]/PAR[1][i_cent] << "\t" << c2ref_err[i_cent] << "\t" << chi2[i_cent] << "\t" <<  H_BBC_phi-3 << "\t" << chi2[i_cent]/(H_FVTX_phi-3) << "\t" << c2ref[i_cent]/c1ref[i_cent] << "\n";


	f1_com[i_cent]->SetParameter(0, c1ref[i_cent]);
	f1_com[i_cent]->Draw("same");
	f2_com[i_cent]->SetParameter(0, c2ref[i_cent]);
	f2_com[i_cent]->Draw("same");
	f3_com[i_cent]->SetParameter(0, p2);
	f3_com[i_cent]->Draw("same");
	f4_com[i_cent]->SetParameter(0, p3);
	f4_com[i_cent]->Draw("same");	

	sprintf(SaveHisto, "%s/Subtructed_C%d.pdf", SUB, i_cent);
	C_1->SaveAs(SaveHisto);
  }
  SavePar_ref.close();
  SavePar_temp.close();
}

  //Multi canvas
  ComCanv = new TCanvas("combined", "combined", 300*NCENTANA, 600);
  ComCanv->Divide(NCENTANA);

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	sprintf(SaveHisto, "ALL/combined/%s_Combined.pdf", SUB);

	ComCanv->cd(i_cent+1);
	histo_original[i_cent]->Draw("ep");
	l_1->Draw("same");
	histo_original[i_cent]->Fit(fit_com[i_cent], "RM");

	f1_com[i_cent]->SetParameter( 0, fit_com[i_cent]->GetParameter(0) );
	f1_com[i_cent]->Draw("same");
	f2_com[i_cent]->SetParameter( 0, fit_com[i_cent]->GetParameter(1) );
	f2_com[i_cent]->Draw("same");
	f3_com[i_cent]->SetParameter( 0, fit_com[i_cent]->GetParameter(2) );
	f3_com[i_cent]->Draw("same");
	f4_com[i_cent]->SetParameter( 0, fit_com[i_cent]->GetParameter(3) );
	f4_com[i_cent]->Draw("same");
	
	ComCanv->SaveAs(SaveHisto);
  }	// i_cent End, Combined canvas drawing End

  delete ComCanv;


  return 0;
}




