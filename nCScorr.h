#ifndef nCSCorr_h_
#define nCSCorr_h_

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TProfile.h"
#include <iostream>

const int nCentBins = 4;

class nCScorr{
	
	public:
		nCScorr(bool ispp);
		double getCorrection(bool ispp, int nCScand, int hiBin, double jtpt, double jteta);	
	
	private:
		TF1 *f_param_a0[nCentBins];
		TF1 *f_param_a1[nCentBins];
		TF1 *flat_corr[nCentBins];
		TF1 *f_pol1[nCentBins];
		double* centBins;
		TFile *fin;
		TFile *fres;
	    double corrpt;
};

nCScorr::nCScorr(bool ispp){
	
	if(ispp) {
		fin = new TFile("new_corr_files/ncscorrfactors_pp5TeV_Apr18.root");
	    fres = new TFile("new_corr_files/fpp_reco_res_Apr18.root");
	}
	else {
		fin = new TFile("/home/dhanush/Documents/JEC/local/corr_files_May29/corr_file_May29.root");
		//fin = new TFile("new_corr_files/ncscorrfactors_PbPb5TeV_Apr17.root");
	    //fres = new TFile("new_corr_files/f_reco_res_Apr17.root");
	}
	if(!fin) std::cout << "Input file for nCScorrs not found!! Aborting!!" << endl;
	else std::cout << "using input file: "<< fin->GetName()/* << " and " << fres->GetName() */<<std::endl;
	for(int i=0; i<nCentBins; i++){
      f_param_a0[i] = (TF1*)fin->Get(Form("f2_par0_cent%d",i))->Clone(Form("f2_par0_cent%d",i));
      f_param_a1[i] = (TF1*)fin->Get(Form("f2_par1_cent%d",i))->Clone(Form("f2_par1_cent%d",i));
      flat_corr[i] = (TF1*)fin->Get(Form("f_flatcorr_cent%d",i))->Clone(Form("f_flatcorr_cent%d",i));
      //f_pol1[i] = (TF1*)fres->Get(Form("f_reco_ratio_cent%d",i))->Clone(Form("f_reco_ratio_cent%d",i));
      f_pol1[i] = (TF1*)fin->Get(Form("f_reco_ratio_cent%d",i))->Clone(Form("f_reco_ratio_cent%d",i));
	}
	
	centBins = new double[nCentBins+1];
	
	double tempCentBins[nCentBins+1] = {0,20,60,100,200};
	
	for(int i=0; i<=nCentBins; i++){
		centBins[i] = tempCentBins[i];
	}
	
}

double nCScorr::getCorrection(bool ispp, int nCScand, int hiBin, double jtpt, double jteta){
	
	if(!fin){ std::cout << "Correction file is not loaded! Returning 1" << std::endl; return 1; }
	if(abs(jteta)>1.6) return -1;
	if(jtpt>600 || jtpt<50) return -1;
	if(hiBin>200 || hiBin<0){ std::cout << "Warning! hiBin is not between 0 and 200!! (=" << hiBin << ")" << std::endl; return -1; }
	
	int centBin=0;
	if(!ispp){
  	  while(hiBin>centBins[centBin+1] && centBin<nCentBins-1){
		centBin++;
	  }
	}

	//Apply nCS corrections
    double p0_cs2 = f_param_a0[centBin]->Eval(jtpt);
    double p1_cs2 = f_param_a1[centBin]->Eval(jtpt);
    
    double corr_factor = 1. + p1_cs2*(nCScand - p0_cs2);

    corrpt = jtpt/corr_factor;

	//Apply residual pol1 corrections to flatten 
    corrpt /= (f_pol1[centBin]->Eval(jtpt));

	//Apply residual flat iterations to properly close refpt
	corrpt /= flat_corr[centBin]->GetParameter(0);	

	return corrpt;
	
}

#endif
