#include <stdio.h>
#include <cmath>
#include <iostream>

#include "sys/types.h"
#include "sys/stat.h"
#include "unistd.h"
#include "def_Const.h"

TH1F *href;
const float PI=acos(-1.);
const int tmp_phi=H_FVTX_phi/2;	//H_FVTX_phi/2
const int Ch_bin=(NPTBIN/2);	//Ch_bin=(NPTBIN/2)+1;
int NBIN_phi=0;

double fit_func(double *x, double *par)
{
  double xx = x[0];
  double yy = par[0]*(href->GetBinContent(href->FindBin(xx))) + par[1];
  return yy;
}

int CorCal_CNT_FVTN()
{
  char SUB[16]="CNT_FVTN";

  gStyle->SetOptStat(0);
  char ref_name[300];
  char SaveHisto[300];
  char output_name[300];
  sprintf(ref_name, "%s/CorHisto_%d_%s.root", ref_dir, ref_run, SUB);

  double c1[NCENTANA][NPTBIN], c2[NCENTANA][NPTBIN], c3[NCENTANA][NPTBIN], c4[NCENTANA][NPTBIN];
  double c1err[NCENTANA][NPTBIN], c2err[NCENTANA][NPTBIN], c3err[NCENTANA][NPTBIN], c4err[NCENTANA][NPTBIN];
  double c1tmp[NCENTANA][NPTBIN], c2tmp[NCENTANA][NPTBIN];
  double c1tmp_err[NCENTANA][NPTBIN], c2tmp_err[NCENTANA][NPTBIN];
  double c1ref[NCENTANA][NPTBIN], c2ref[NCENTANA][NPTBIN];
  double c1ref_err[NCENTANA][NPTBIN], c2ref_err[NCENTANA][NPTBIN];

  TFile *rootopen = TFile::Open(input_name, "READ");
  TFile *ref_open = TFile::Open(ref_name, "READ");

  cout << "Open file : " << input_name << endl;
  cout << "Open file : " << ref_name << endl;

  TCanvas *C_1  = new TCanvas("C_1", "C_1", 500, 400);
  TLine	*l_1 = new TLine(-0.5*PI, 1, 1.5*PI, 1);
  TF1	*f2_tmp = new TF1("f2_tmp" , "1 + 2*[0]*cos(2*x)", -0.5*PI, 1.5*PI);
  TH1	*f_fit	= new TH1D("f_fit", "f_fit", H_FVTX_phi, -0.5*PI, 1.5*PI);
  f_fit->SetMarkerStyle(24);
  f_fit->SetMarkerSize(1.4);
  TH1	*f_fit2	= new TH1D("f_fit2", "f_fit2", H_FVTX_phi/2, -0.5*PI, 1.5*PI);
  f_fit2->SetMarkerStyle(24);
  f_fit2->SetMarkerSize(1.4);
  f2_tmp->SetLineColor(kRed);
  f2_tmp->SetLineStyle(7);

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

  TH1		*histo_original	[NCENTANA][NPTBIN];
  TH1		*histo_original_Sub	[NCENTANA][NPTBIN];
  char hname[200];
  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	for(int ipt_bin=0; ipt_bin<Ch_bin; ipt_bin++)
	{
	  if(i_cent==NCENTANA-1 && saveOption==1) sprintf(hname, "COR_CNT_FVTN_C%d_pT%d", 0, ipt_bin);
	  else sprintf(hname, "SIG_COR_CNT_FVTN_C%d_pT%d", i_cent, ipt_bin);
	  histo_original[i_cent][ipt_bin] = new TH1D(hname, hname, H_FVTX_phi, -0.5*PI, 1.5*PI);

	  sprintf(hname, "SUB_CNT_FVTN_C%d_pT%d", i_cent, ipt_bin);
	  histo_original_Sub[i_cent][ipt_bin] = new TH1D(hname, hname, H_FVTX_phi, -0.5*PI, 1.5*PI);
	}
	for(int ipt_bin=Ch_bin; ipt_bin<NPTBIN; ipt_bin++)
	{
	  if(i_cent==NCENTANA-1 && saveOption==1) sprintf(hname, "COR_CNT_FVTN_C%d_pT%d", 0, ipt_bin);
	  else sprintf(hname, "SIG_COR_CNT_FVTN_C%d_pT%d", i_cent, ipt_bin);
	  histo_original[i_cent][ipt_bin] = new TH1D(hname, hname, tmp_phi, -0.5*PI, 1.5*PI);

	  sprintf(hname, "SUB_CNT_FVTN_C%d_pT%d", i_cent, ipt_bin);
	  histo_original_Sub[i_cent][ipt_bin] = new TH1D(hname, hname, tmp_phi, -0.5*PI, 1.5*PI);
	}

  }

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	for(int ipt_bin=0; ipt_bin<NPTBIN; ipt_bin++)
	{
	  histo_original[i_cent][ipt_bin]->GetXaxis()->SetTitle("delta phi");
	  histo_original[i_cent][ipt_bin]->SetMarkerColor(kBlack);
	  histo_original[i_cent][ipt_bin]->SetLineColor(kBlack);
	  histo_original[i_cent][ipt_bin]->SetLineWidth(2);
	  histo_original[i_cent][ipt_bin]->SetMarkerStyle(20);
	  histo_original[i_cent][ipt_bin]->SetAxisRange(0.82, 1.3, "Y");

	  histo_original_Sub[i_cent][ipt_bin]->GetXaxis()->SetTitle("delta phi");
	  histo_original_Sub[i_cent][ipt_bin]->SetMarkerColor(kBlack);
	  histo_original_Sub[i_cent][ipt_bin]->SetLineColor(kBlack);
	  histo_original_Sub[i_cent][ipt_bin]->SetMarkerStyle(20);
	  histo_original_Sub[i_cent][ipt_bin]->SetAxisRange(0.995, 1.005, "Y");
	}
  }

  sprintf(output_name, "CorHisto_%s.root", SUB);  
  sprintf(param_name, "savePar/parameter_%s.txt", SUB);
  TFile *outfile = new TFile(output_name, "recreate");  

  SavePar.open(param_name, ios::trunc);
  cout << "Parameter save file : " << param_name << endl;

  TH2D *sumR  = new TH2D("sumR for low pT",  "Real pT sum histogram for low pT", H_FVTX_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *sumM  = new TH2D("sumM for low pT",  "Mix pT sum histogram for low pT", H_FVTX_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *tmpsumR  = new TH2D("tmpsumR for low pT",  " ", H_FVTX_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *tmpsumM  = new TH2D("tmpsumM for low pT",  " ", H_FVTX_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);

  TH2D *sumR_h = new TH2D("sumR",  "Real pT sum histogram", tmp_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *sumM_h = new TH2D("sumM",  "Mix pT sum histogram", tmp_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *tmpsumR_h = new TH2D("tmpsumR",  " ", tmp_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);
  TH2D *tmpsumM_h = new TH2D("tmpsumM",  " ", tmp_phi, -0.5*PI, 1.5*PI, H_FVTX_eta, -3, 3);

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	for(int ipt_bin=0; ipt_bin<NPTBIN; ipt_bin++)
	{
	  for(int i_zvtx=0; i_zvtx<NZANA; i_zvtx++)
	  {
		if(ipt_bin<Ch_bin)
		{
			TH2D *h1 = (TH2D*)rootopen->Get(Form("Real_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin));
			TH2D *h2 = (TH2D*)rootopen->Get(Form("Mix_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin));

			NBIN_phi = h1->GetNbinsX();
			for (int i_phi=0; i_phi<H_FVTX_phi; i_phi++)
			{
		  		float tmpaddhistoR=0, tmpaddhistoM=0, tmpaddhistoR2=0, tmpaddhistoM2=0;
		  		for(int ieta=H_FVTX_eta/2; ieta<H_FVTX_eta; ieta++)
		  		{
					tmpaddhistoR  = h1->GetBinContent(i_phi+1, ieta+1);
					tmpaddhistoM  = h2->GetBinContent(i_phi+1, ieta+1);
					tmpsumR->SetBinContent(i_phi+1, ieta+1, tmpaddhistoR);
					tmpsumM->SetBinContent(i_phi+1, ieta+1, tmpaddhistoM);
		  		} //ieta
			}//i_phi
        	sumR->Add(tmpsumR);
			sumM->Add(tmpsumM);
		}

		if(ipt_bin>=Ch_bin)
		{
			TH2D *h1 = (TH2D*)rootopen->Get(Form("Real_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin));
			TH2D *h2 = (TH2D*)rootopen->Get(Form("Mix_dPhi_dEta_CNT_FVT_C%d_Z%d_Pt%d", i_cent, i_zvtx, ipt_bin));

			NBIN_phi = h1->GetNbinsX();
			for (int i_phi=0; i_phi<H_FVTX_phi; i_phi++)
			{
		  		float tmpaddhistoR=0, tmpaddhistoM=0, tmpaddhistoR2=0, tmpaddhistoM2=0;
		  		for(int ieta=H_FVTX_eta/2; ieta<H_FVTX_eta; ieta++)
		  		{
					tmpaddhistoR  = h1->GetBinContent(i_phi+1, ieta+1);
					tmpaddhistoM  = h2->GetBinContent(i_phi+1, ieta+1);
					tmpsumR_h->SetBinContent(i_phi+1, ieta+1, tmpaddhistoR);
					tmpsumM_h->SetBinContent(i_phi+1, ieta+1, tmpaddhistoM);
		  		} //ieta
			}//i_phi
        	sumR_h->Add(tmpsumR_h);
			sumM_h->Add(tmpsumM_h);
		}

	  }//iz_ana
	  
	  char RealHisto[300];
	  char MixHisto[300];
	  sprintf(RealHisto, "ALL/%s_Real_C%d_Pt%d.pdf", SUB, i_cent, ipt_bin);
	  sprintf(MixHisto,  "ALL/%s_Mix_C%d_Pt%d.pdf",  SUB, i_cent, ipt_bin);

	  float sum_R = 0;
	  float sum_M = 0;

	  TH1D *p_R = new TH1D("p_R", "projection_R", NBIN_phi, -0.5*PI, 1.5*PI);
	  TH1D *p_M = new TH1D("p_M", "projection_M", NBIN_phi, -0.5*PI, 1.5*PI);

	  TH1D *p_R_tmp = sumR->ProjectionX();
	  TH1D *p_M_tmp = sumM->ProjectionX();
	  sum_R = sumR->GetSum();
	  sum_M = sumM->GetSum();

	  if(ipt_bin<Ch_bin)
	  {	
		TH1D *p_R_tmp = sumR->ProjectionX();
		TH1D *p_M_tmp = sumM->ProjectionX();
		sum_R = sumR->GetSum();
		sum_M = sumM->GetSum();

		for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
		{
			p_R->SetBinContent(i_phi+1, p_R_tmp->GetBinContent(i_phi+1));
			p_M->SetBinContent(i_phi+1, p_M_tmp->GetBinContent(i_phi+1));
		}

	  }
	  else if(ipt_bin>=Ch_bin)
	  {
		TH1D *p_R_tmp = sumR_h->ProjectionX();
		TH1D *p_M_tmp = sumM_h->ProjectionX();
		sum_R = sumR_h->GetSum();
		sum_M = sumM_h->GetSum();
		for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
		{
		  p_R->SetBinContent(i_phi+1, p_R_tmp->GetBinContent(i_phi+1));
		  p_M->SetBinContent(i_phi+1, p_M_tmp->GetBinContent(i_phi+1));
		}
	  }

	  p_R->Draw();
	  C_1->SaveAs(RealHisto);
	  p_M->Draw();
	  C_1->SaveAs(MixHisto);

	  for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
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
		histo_original[i_cent][ipt_bin]->SetBinContent(i_phi+1, valRatio);
		histo_original[i_cent][ipt_bin]->SetBinError(i_phi+1, dsig);
	  }// i_phi End

	  histo_original[i_cent][ipt_bin]->Draw("e0");
	  histo_original[i_cent][ipt_bin]->Write();
	  histo_original[i_cent][ipt_bin]->Fit(fit_com[i_cent], "RM");
	  l_1 ->Draw("same");
	  // Save parameter
	  c1[i_cent][ipt_bin]=fit_com[i_cent]->GetParameter(0);
	  c2[i_cent][ipt_bin]=fit_com[i_cent]->GetParameter(1);
	  c3[i_cent][ipt_bin]=fit_com[i_cent]->GetParameter(2);
	  c4[i_cent][ipt_bin]=fit_com[i_cent]->GetParameter(3);
	  c1err[i_cent][ipt_bin]=fit_com[i_cent]->GetParError(0);
	  c2err[i_cent][ipt_bin]=fit_com[i_cent]->GetParError(1);
	  c3err[i_cent][ipt_bin]=fit_com[i_cent]->GetParError(2);
	  float chi2[NCENTANA][NPTBIN]={0};

	  // Calculate chi^2/NDF
	  for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
	  {
		double data=histo_original[i_cent][ipt_bin]->GetBinContent(i_phi+1);
		double tmpx=histo_original[i_cent][ipt_bin]->GetBinCenter(i_phi+1);
		double tmperr=histo_original[i_cent][ipt_bin]->GetBinError(i_phi+1);
		double fitval=fit_com[i_cent]->Eval(tmpx);
		double tmpval=(fitval-data)*(fitval-data)/(tmperr*tmperr);
		chi2[i_cent][ipt_bin] = chi2[i_cent][ipt_bin]+tmpval;
	  }

	  SavePar << i_cent << "\t" << ipt_bin << "\t" << c1[i_cent][ipt_bin] << "\t" << c1err[i_cent][ipt_bin] << "\t" << c2[i_cent][ipt_bin] << "\t" << c2err[i_cent][ipt_bin] << "\t" << chi2[i_cent][ipt_bin] << "\t" << H_FVTX_phi-3 << "\t" << chi2[i_cent][ipt_bin]/(H_FVTX_phi-3) << "\t" << c2[i_cent][ipt_bin]/c1[i_cent][ipt_bin] << "\n" ;

	  f1_com[i_cent]->SetParameter(0, c1[i_cent][ipt_bin]);
	  f1_com[i_cent]->Draw("same");
	  f2_com[i_cent]->SetParameter(0, c2[i_cent][ipt_bin]);
	  f2_com[i_cent]->Draw("same");
	  f3_com[i_cent]->SetParameter(0, c3[i_cent][ipt_bin]);
	  f3_com[i_cent]->Draw("same");
	  f4_com[i_cent]->SetParameter(0, c4[i_cent][ipt_bin]);
	  f4_com[i_cent]->Draw("same");

	  sprintf(SaveHisto, "ALL/%s_Cor_C%d_pT%d.pdf", SUB, i_cent, ipt_bin);
	  C_1->SaveAs(SaveHisto);
	  sumR->Reset();
	  sumM->Reset();
	  sumR_h->Reset();
	  sumM_h->Reset();
	  delete p_R;
	  delete p_M;
	}//ipt_bin
  }//CNT_FVTN
  delete sumR;
  delete sumM;
  delete sumR_h;
  delete sumM_h;
  SavePar.close();



  //Subtraction
  if(subtraction!=0)
  {
  double PAR[3][NCENTANA][NPTBIN];
  TH1F	*C_ref [NPTBIN];
  TF1	*f_ref[NCENTANA][NPTBIN];

  sprintf(param_name, "savePar/parameter_%s_ref.txt", SUB);
  SavePar_ref.open(param_name, ios::trunc);
  sprintf(param_name, "savePar/parameter_%s_temp.txt", SUB);
  SavePar_temp.open(param_name, ios::trunc);

  for(int ipt_bin=0; ipt_bin<NPTBIN; ipt_bin++)
  {
	NBIN_phi = histo_original[0][ipt_bin]->GetNbinsX();
	href = new TH1F("Reference correlation function", "href", NBIN_phi, -0.5*PI, 1.5*PI);
	if(subtraction==2) C_ref [ipt_bin] = (TH1F*)histo_original[NCENTANA-1][ipt_bin];
	else if(subtraction==1) C_ref [ipt_bin] = (TH1F*)ref_open->Get(Form("COR_%s_C0_pT%d", SUB, ipt_bin));

	for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
	{
		href->SetBinContent(i_phi+1, C_ref[ipt_bin]->GetBinContent(i_phi+1));
		href->SetBinError(i_phi+1, C_ref[ipt_bin]->GetBinError(i_phi+1));
	}
	for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  	{
		f_ref[i_cent][ipt_bin] = new TF1("f_ref", *fit_func, -0.5*PI, 1.5*PI, 2);

		histo_original[i_cent][ipt_bin]->Fit(f_ref[i_cent][ipt_bin], "N");

		PAR[0][i_cent][ipt_bin]=f_ref [i_cent][ipt_bin]->GetParameter(0);
		PAR[1][i_cent][ipt_bin]=f_ref [i_cent][ipt_bin]->GetParameter(1);
		if(NBIN_phi==H_FVTX_phi){
		  for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
		  {
		  double a_tmp=PAR[0][i_cent][ipt_bin];
		  double b_tmp=PAR[1][i_cent][ipt_bin];
		  double c_tmp=C_ref[ipt_bin]->GetBinContent(i_phi+1);
		  double d_tmp=a_tmp*c_tmp + b_tmp;
		  f_fit->SetBinContent(i_phi+1, d_tmp );
		  f_fit->SetBinError(i_phi+1, 0 );
		  }
		}
		else
		{
		  for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
		  {
		  double a_tmp=PAR[0][i_cent][ipt_bin];
		  double b_tmp=PAR[1][i_cent][ipt_bin];
		  double c_tmp=C_ref[ipt_bin]->GetBinContent(i_phi+1);
		  double d_tmp=a_tmp*c_tmp + b_tmp;
		  f_fit2->SetBinContent(i_phi+1, d_tmp );
		  f_fit2->SetBinError(i_phi+1, 0 );
		  }
		}

		f2_tmp->SetParameter(0, c2[i_cent][ipt_bin]-PAR[0][i_cent][ipt_bin]*c2[NCENTANA-1][ipt_bin]);

		histo_original[i_cent][ipt_bin]->Draw("ep");
		f2_tmp->Draw("same");
		l_1->Draw("same");
		if(NBIN_phi==H_FVTX_phi)		f_fit->Draw("ep, same");
		else		f_fit2->Draw("ep, same");

		sprintf(SaveHisto, "%s/GetF_C%d_pT%d.pdf", SUB, i_cent, ipt_bin);
		C_1->SaveAs(SaveHisto);


		for(int i_phi=0; i_phi<NBIN_phi; i_phi++)
		{
		  double ii = histo_original[i_cent][ipt_bin]->GetBinContent(i_phi+1);
		  double tmp_v=C_ref[ipt_bin]->GetBinContent(i_phi+1);
		  double jj = (PAR[0][i_cent][ipt_bin]*tmp_v)+ PAR[1][i_cent][ipt_bin];
		  double err_ref=sqrt( (histo_original[i_cent][ipt_bin]->GetBinError(i_phi+1)*histo_original[i_cent][ipt_bin]->GetBinError(i_phi+1)) + PAR[0][i_cent][ipt_bin]*PAR[0][i_cent][ipt_bin]*(C_ref[ipt_bin]->GetBinError(i_phi+1)*C_ref[ipt_bin]->GetBinError(i_phi+1)) );
	  
		  histo_original_Sub[i_cent][ipt_bin]->SetBinContent(i_phi+1, ii-jj+1);
		  histo_original_Sub[i_cent][ipt_bin]->SetBinError(i_phi+1, err_ref);
		  histo_original_Sub[i_cent][ipt_bin]->SetMarkerStyle(33);
		  histo_original_Sub[i_cent][ipt_bin]->SetLineColor(kPink);
		  histo_original_Sub[i_cent][ipt_bin]->SetMarkerColor(kPink);
		}
		histo_original_Sub[i_cent][ipt_bin]->Draw("p");
		l_1->Draw("same");
		histo_original_Sub[i_cent][ipt_bin]->Fit(fit_com[i_cent], "RM");

		c1ref[i_cent][ipt_bin]	=fit_com[i_cent]->GetParameter(0);
		c2ref[i_cent][ipt_bin]	=fit_com[i_cent]->GetParameter(1);
		float p2				=fit_com[i_cent]->GetParameter(2);
		float p3				=fit_com[i_cent]->GetParameter(3);
		c1ref_err[i_cent][ipt_bin]=fit_com[i_cent]->GetParError(0);
		c2ref_err[i_cent][ipt_bin]=fit_com[i_cent]->GetParError(1);
		float e2				=fit_com[i_cent]->GetParError(2);
		float chi2[NCENTANA][NPTBIN]={0};

		SavePar_ref << i_cent << "\t" << ipt_bin << "\t" << c1ref[i_cent][ipt_bin] << "\t" << c1ref_err[i_cent][ipt_bin]  << "\t" << c2ref[i_cent][ipt_bin] << "\t" << c2ref_err[i_cent][ipt_bin] << "\t" << chi2[i_cent][ipt_bin] << "\t" <<  H_BBC_phi-3 << "\t" << chi2[i_cent][ipt_bin]/(H_BBC_phi-3) << "\t" << c2ref[i_cent][ipt_bin]/c1ref[i_cent][ipt_bin] << "\n";

		cout << c2[i_cent][ipt_bin] -(PAR[0][i_cent][ipt_bin] *c2[NCENTANA-1][ipt_bin] ) << "  " << c2ref[i_cent][ipt_bin]  << endl;

		SavePar_temp << i_cent << "\t" << ipt_bin << "\t" << c1ref[i_cent][ipt_bin] << "\t" << c1ref_err[i_cent][ipt_bin]  << "\t" << c2ref[i_cent][ipt_bin]/PAR[1][i_cent][ipt_bin] << "\t" << c2ref_err[i_cent][ipt_bin] << "\t" << chi2[i_cent][ipt_bin] << "\t" <<  NBIN_phi-3 << "\t" << chi2[i_cent][ipt_bin]/(NBIN_phi-3) << "\t" << c2ref[i_cent][ipt_bin]/c1ref[i_cent][ipt_bin] << "\n";

		f1_com[i_cent]->SetParameter(0, c1ref[i_cent][ipt_bin]);
		f1_com[i_cent]->Draw("same");
		f2_com[i_cent]->SetParameter(0, c2ref[i_cent][ipt_bin]);
		f2_com[i_cent]->Draw("same");
		f3_com[i_cent]->SetParameter(0, p2);
		f3_com[i_cent]->Draw("same");
		f4_com[i_cent]->SetParameter(0, p3);

		sprintf(SaveHisto, "%s/Subtructed_C%d_pT%d.pdf", SUB, i_cent, ipt_bin);
		C_1->SaveAs(SaveHisto);
	}//ipt_bin
	C_ref[ipt_bin]->Reset();
	href->Reset();
	delete href;
  }//i_cent

  SavePar_ref.close();
  SavePar_temp.close();

  }
//MultiCanvas

  ComCanv= new TCanvas("combined", "combined", 300*((NPTBIN+1)/2), 600*2);
  ComCanv->Divide((NPTBIN+1)/2, 2);

  for(int i_cent=0; i_cent<NCENTANA; i_cent++)
  {
	for(int ipt_bin=0; ipt_bin<NPTBIN; ipt_bin++)
	{
	  sprintf(SaveHisto, "ALL/combined/%s_Combined_%d.pdf", SUB, i_cent);

	  ComCanv->cd(ipt_bin+1);
	  histo_original[i_cent][ipt_bin]->Draw("ep");
	  l_1->Draw("same");
	  histo_original[i_cent][ipt_bin]->Fit(fit_com[ipt_bin], "RM");

	  f1_com[ipt_bin]->SetParameter( 0, fit_com[ipt_bin]->GetParameter(0) );
	  f1_com[ipt_bin]->Draw("same");
	  f2_com[ipt_bin]->SetParameter( 0, fit_com[ipt_bin]->GetParameter(1) );
	  f2_com[ipt_bin]->Draw("same");	
	  ComCanv->SaveAs(SaveHisto);
	}
  }	// i_cent End, Combined canvas drawing End

  delete ComCanv;


  return 0;
}


