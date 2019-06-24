void Run(){

	//gSystem->Load("pAu_Cor_C.so");
	gROOT->ProcessLine(".L pAu_Cor.C+g");
	pAu_Cor();

	return;

}
