const int TAXI		= 10883;
const int NCENTANA	= 7;
const int NEVE		= 8;	//Event buffer
const int NZANA		= 10;
const int NPTBIN	= 1;    //Even number would be better..
const int NARM		= 2;
//histogram const
const int H_BBC_phi = 12;
const int H_BBC_eta = 60;
const int H_BBCCNT_eta = 100;
const int H_CNT_eta = 20;
const int H_FVTX_phi= 24;
const int H_FVTX_eta= 60;
//QA const
const int H_PC3_DZ  =60;
const int H_PC3_DPHI=60;
const int H_PC3 = 50;
//int centlim[NCENTANA]={17, 12, 8, 5, 0}; //5, 20, 40, 60, 100
//int centlim[NCENTANA]={17, 12, 5, 0}; //5, 20, 60, 100
//int centlim[NCENTANA]={12, 0}; //5, 20, 60, 100
//float cent  [NCENTANA] = {20, 84};//Centrality
float cent  [NCENTANA] = {1, 3, 5, 20, 40, 60, 84};//Centrality

//Reference
const int saveOption=0; //0=no need to save Ref histo
const int subtraction=2; //1=MB, 2=peripheral
char input_name[300] = "../../dAu_dEta_sys/dAu_7m1pt_eta05_sys0.root";
char ref_dir[300] = "~/Dropbox/IPA_working/2Ptcl_Cor/CorCal_pp_novtx/20170323/";
const int ref_run = 10884;
