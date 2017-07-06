#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>

#define nCBins 5
#define nptBins 55
#define npfbins 21
const int npt_histoBins = 15;

char saythis[500];

using namespace std;

bool isdata = false;

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};
Double_t pt_bounds[15] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 170., 220., 280., 370., 500.};

void test_spectra(){

  //TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_data_May22.root");
  TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/pptest_histos_data_Jun14.root");
  //TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_MC_Jun12.root");
  TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/JFF_corrections/pptest_histos_MC_Jun14.root"); 
  TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_data_Jun14.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_data_Jun9_fixed.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_MC_Jun26.root");
  TFile *raa_raghav = TFile::Open("/home/dhanush/Documents/JFF_corrections/JetRAA_datapoints.root");

  TGraphAsymmErrors *g_raa[nCBins];
  
  TH1D *h_reco_full_un[nCBins];
  TH1D *h_reco_full_n[nCBins];
  TH1D *h_gen_full_un[nCBins];
  TH1D *h_gen_full_n[nCBins];
  TH1D *h_reco_corr_un[nCBins];
  TH1D *h_reco_corr_n[nCBins];
  TH1D *h_reco_full_MC_un[nCBins];
  TH1D *h_reco_full_MC_n[nCBins];
  TH1D *h_reco_corr_MC_un[nCBins];
  TH1D *h_reco_corr_MC_n[nCBins];

  TH1D *h_vz_MC[nCBins];
  TH1D *h_vz_MC_nrw[nCBins];
  TH1D *h_hibin;
  TH1D *h_hibin_nrw;

  TH1D *h_eta_full_MC[nCBins];
  TH1D *h_phi_full_MC[nCBins];
  TH1D *h_eta_full_MC_nrw[nCBins];
  TH1D *h_phi_full_MC_nrw[nCBins];  

  TH1D *h_reco_ratio[nCBins];
  TH1D *h_reco_ratio_MC[nCBins];
  TH1D *h_corr_gen_ratio[nCBins];
  TH1D *h_uncorr_gen_ratio[nCBins];
  TH1D *h_dataMC_corr_ratio[nCBins];
  TH1D *h_dataMC_uncorr_ratio[nCBins];
  TH1D *h_dataMC_corr_ratio_raascaled[nCBins];
  TH1D *h_dataMC_uncorr_ratio_raascaled[nCBins];  
  TH2F *h_ncs_pt_data[nCBins];
  TH2F *h_npf_pt_data[nCBins];
  TH2F *h_ncs_pt_MC[nCBins];
  TH1D *h_ncs_data[nCBins];
  TH1D *h_npf_data[nCBins];  
  TH1D *h_ncs_MC[nCBins];

  TH1D *h_corr_corr_ratio[nCBins];
  TH1D *h_corr_corr_ratio_MC[nCBins];
  TH1D *h_uncorr_uncorr_ratio[nCBins];
  TH1D *h_uncorr_uncorr_ratio_MC[nCBins];
  TH1D *h_gen_gen_ratio[nCBins];

  TH1D *h_corr_gen_corr_gen_ratio[nCBins];
  TH1D *h_uncorr_gen_uncorr_gen_ratio[nCBins];

  TF1 *f_raa[nCBins];

  h_hibin = (TH1D*)closure_histos_PbPb_MC->Get("h_hibin")->Clone("h_hibin");
  h_hibin_nrw = (TH1D*)closure_histos_PbPb_MC->Get("h_hibin_nrw")->Clone("h_hibin_nrw");

  for(int ibin=0;ibin<nCBins;ibin++){
    
    sprintf(saythis,"h_gen_full_new_cent%d",ibin);
    h_gen_full_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_gen_full_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_new_cent%d",ibin);
    h_reco_full_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_full_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_new_cent%d",ibin);
    h_reco_corr_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_corr_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_MC_new_cent%d",ibin);
    h_reco_full_MC_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_full_MC_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_MC_new_cent%d",ibin);
    h_reco_corr_MC_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_corr_MC_n[ibin]->Sumw2();

    sprintf(saythis,"h_dataMC_uncorr_ratio_raascaled_cent%d",ibin);
    h_dataMC_uncorr_ratio_raascaled[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_dataMC_uncorr_ratio_raascaled[ibin]->Sumw2();

    sprintf(saythis,"h_dataMC_corr_ratio_raascaled_cent%d",ibin);
    h_dataMC_corr_ratio_raascaled[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_dataMC_corr_ratio_raascaled[ibin]->Sumw2();    

    if(ibin==0){
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_pp_data->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin])); 
      h_reco_full_un[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_un_"+cent[ibin]));
      h_reco_corr_un[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_un_"+cent[ibin]));    
      h_reco_full_MC_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_MC_un_"+cent[ibin]));      
      h_reco_corr_MC_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_MC_un_"+cent[ibin]));    
      if(!isdata) h_gen_full_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_cent"+cent[ibin]))->Clone((TString)("h_gen_full_un_"+cent[ibin]));
    }

    else if(ibin > 0){     
       
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_npf_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_npf_pt_cent"+cent[ibin-1]))->Clone((TString)("h_npf_pt_data_"+cent[ibin]));       
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin]));     
      h_reco_full_un[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_un_"+cent[ibin]));
      h_vz_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_vz_cent"+cent[ibin-1]))->Clone((TString)("h_vz_"+cent[ibin]));
      h_vz_MC_nrw[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_vz_nrw_cent"+cent[ibin-1]))->Clone((TString)("h_vz_nrw_"+cent[ibin]));
      h_eta_full_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_cent"+cent[ibin-1]))->Clone((TString)("h_eta_full_"+cent[ibin]));
      h_phi_full_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_phi_full_cent"+cent[ibin-1]))->Clone((TString)("h_phi_full_"+cent[ibin])); 
      h_eta_full_MC_nrw[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_nrw_cent"+cent[ibin-1]))->Clone((TString)("h_eta_full_nrw_"+cent[ibin]));
      h_phi_full_MC_nrw[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_phi_full_nrw_cent"+cent[ibin-1]))->Clone((TString)("h_phi_full_nrw_"+cent[ibin]));
      h_reco_corr_un[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_un_"+cent[ibin]));    
      h_reco_full_MC_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_MC_un_"+cent[ibin]));      
      h_reco_corr_MC_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_MC_un_"+cent[ibin]));    
      if(!isdata) h_gen_full_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_cent"+cent[ibin-1]))->Clone((TString)("h_gen_full_un_"+cent[ibin]));    
    }     

    sprintf(saythis,"f_raa_cent%d",ibin);
    f_raa[ibin] = new TF1(saythis, "pol3", 100., 300.);
    f_raa[ibin] ->SetParameter(0,0.629);
    f_raa[ibin] ->SetParameter(1,0.00392094);
    f_raa[ibin] ->SetParameter(2,-5.876e-6);
    f_raa[ibin] ->SetParameter(3,-7.4198e-9);
  }

    g_raa[1] = (TGraphAsymmErrors*)raa_raghav->Get("RAA_R4_cent0")->Clone("RAA_R4_cent0");
    g_raa[2] = (TGraphAsymmErrors*)raa_raghav->Get("RAA_R4_cent2")->Clone("RAA_R4_cent1");
    g_raa[3] = (TGraphAsymmErrors*)raa_raghav->Get("RAA_R4_cent3")->Clone("RAA_R4_cent2");            
    g_raa[4] = (TGraphAsymmErrors*)raa_raghav->Get("RAA_R4_cent4")->Clone("RAA_R4_cent3");

    for (int i=0;i<g_raa[1]->GetN();i++) g_raa[1]->GetY()[i] *= (1./0.472574);
    for (int i=0;i<g_raa[2]->GetN();i++) g_raa[2]->GetY()[i] *= (1./0.593882); 
    for (int i=0;i<g_raa[3]->GetN();i++) g_raa[3]->GetY()[i] *= (1./0.736287); 
    for (int i=0;i<g_raa[4]->GetN();i++) g_raa[4]->GetY()[i] *= (1./0.841772);      

    g_raa[1]->Fit(f_raa[1],"Q M","",100.,260.);      
    g_raa[2]->Fit(f_raa[2],"Q M","",100.,260.);
    g_raa[3]->Fit(f_raa[3],"Q M","",100.,260.);
    g_raa[4]->Fit(f_raa[4],"Q M","",100.,240.);

  for(int i=0; i<nCBins; i++){

    for(int j=1;j<15;j++){
      //cout<<h_gen_full_un[i]->GetBinContent(j)<<endl;
      h_gen_full_n[i]->SetBinContent(j,h_gen_full_un[i]->GetBinContent(j)/h_gen_full_un[i]->GetBinWidth(j));
      h_gen_full_n[i]->SetBinError(j,h_gen_full_un[i]->GetBinError(j)/h_gen_full_un[i]->GetBinWidth(j));
      h_reco_full_n[i]->SetBinContent(j,h_reco_full_un[i]->GetBinContent(j)/h_reco_full_un[i]->GetBinWidth(j));
      h_reco_full_n[i]->SetBinError(j,h_reco_full_un[i]->GetBinError(j)/h_reco_full_un[i]->GetBinWidth(j));
      h_reco_corr_n[i]->SetBinContent(j,h_reco_corr_un[i]->GetBinContent(j)/h_reco_corr_un[i]->GetBinWidth(j));
      h_reco_corr_n[i]->SetBinError(j,h_reco_corr_un[i]->GetBinError(j)/h_reco_corr_un[i]->GetBinWidth(j));
      h_reco_full_MC_n[i]->SetBinContent(j,h_reco_full_MC_un[i]->GetBinContent(j)/h_reco_full_MC_un[i]->GetBinWidth(j));
      h_reco_full_MC_n[i]->SetBinError(j,h_reco_full_MC_un[i]->GetBinError(j)/h_reco_full_MC_un[i]->GetBinWidth(j));
      h_reco_corr_MC_n[i]->SetBinContent(j,h_reco_corr_MC_un[i]->GetBinContent(j)/h_reco_corr_MC_un[i]->GetBinWidth(j));
      h_reco_corr_MC_n[i]->SetBinError(j,h_reco_corr_MC_un[i]->GetBinError(j)/h_reco_corr_MC_un[i]->GetBinWidth(j));  
    }

    h_gen_full_n[i]->Scale(1./h_gen_full_n[i]->GetBinContent(h_gen_full_n[i]->FindBin(120)));

    h_reco_full_n[i]->Scale(1./h_reco_full_n[i]->GetBinContent(h_reco_full_n[i]->FindBin(120)));
    h_reco_corr_n[i]->Scale(1./h_reco_corr_n[i]->GetBinContent(h_reco_corr_n[i]->FindBin(120)));
    h_reco_full_MC_n[i]->Scale(1./h_reco_full_MC_n[i]->GetBinContent(h_reco_full_MC_n[i]->FindBin(120)));
    h_reco_corr_MC_n[i]->Scale(1./h_reco_corr_MC_n[i]->GetBinContent(h_reco_corr_MC_n[i]->FindBin(120)));
/*
    h_reco_full_n[i]->Scale(0.472574);
    h_reco_corr_n[i]->Scale(0.472574);
    h_reco_full_MC_n[i]->Scale(0.472574);
    h_reco_corr_MC_n[i]->Scale(0.472574);
*/  
    //h_gen_full_n[i]->Scale(1./h_gen_full_n[i]->Integral());
    //h_reco_full_MC_n[i]->Scale(1./h_reco_full_MC_n[i]->Integral(100,500));
    //h_reco_corr_MC_n[i]->Scale(1./h_reco_corr_MC_n[i]->Integral(100,500));
    //if(i==3)h_reco_corr_MC_n[3]->Draw();
  
    sprintf(saythis,"h_reco_ratio_cent%d",i);  
    h_reco_ratio[i] = (TH1D*)h_reco_corr_n[i]->Clone(saythis);
    h_reco_ratio[i]->Divide(h_reco_full_n[i]);
    sprintf(saythis,"h_reco_ratio_MC_cent%d",i);  
    h_reco_ratio_MC[i] = (TH1D*)h_reco_corr_MC_n[i]->Clone(saythis);
    h_reco_ratio_MC[i]->Divide(h_reco_full_MC_n[i]);
    sprintf(saythis,"h_corr_gen_ratio_cent%d",i);  
    h_corr_gen_ratio[i] = (TH1D*)h_reco_corr_MC_n[i]->Clone(saythis);
    h_corr_gen_ratio[i]->Divide(h_gen_full_n[i]);
    sprintf(saythis,"h_uncorr_gen_ratio_cent%d",i);  
    h_uncorr_gen_ratio[i] = (TH1D*)h_reco_full_MC_n[i]->Clone(saythis);
    h_uncorr_gen_ratio[i]->Divide(h_gen_full_n[i]); 
    sprintf(saythis,"h_dataMC_corr_ratio_cent%d",i);  
    h_dataMC_corr_ratio[i] = (TH1D*)h_reco_corr_n[i]->Clone(saythis);
    h_dataMC_corr_ratio[i]->Divide(h_reco_corr_MC_n[i]);
    sprintf(saythis,"h_dataMC_uncorr_ratio_cent%d",i);  
    h_dataMC_uncorr_ratio[i] = (TH1D*)h_reco_full_n[i]->Clone(saythis);
    h_dataMC_uncorr_ratio[i]->Divide(h_reco_full_MC_n[i]);
  
    if(i>0) {
      for(int j=1;j<15;j++){
        h_dataMC_uncorr_ratio_raascaled[i]->SetBinContent(j,0.472574*(h_dataMC_uncorr_ratio[i]->GetBinContent(j)/f_raa[i]->Eval(pt_bounds[j-1])));
        h_dataMC_uncorr_ratio_raascaled[i]->SetBinError(j,0.472574*(h_dataMC_uncorr_ratio[i]->GetBinError(j)/f_raa[i]->Eval(pt_bounds[j-1])));
        h_dataMC_corr_ratio_raascaled[i]->SetBinContent(j,0.472574*(h_dataMC_corr_ratio[i]->GetBinContent(j)/f_raa[i]->Eval(pt_bounds[j-1])));
        h_dataMC_corr_ratio_raascaled[i]->SetBinError(j,0.472574*(h_dataMC_corr_ratio[i]->GetBinError(j)/f_raa[i]->Eval(pt_bounds[j-1])));        
      }

      sprintf(saythis,"h_uncorr_gen_uncorr_gen_ratio_cent%d",i);  
      h_uncorr_gen_uncorr_gen_ratio[i] = (TH1D*)h_uncorr_gen_ratio[i]->Clone(saythis);
      h_uncorr_gen_uncorr_gen_ratio[i]->Divide(h_uncorr_gen_ratio[0]);
      sprintf(saythis,"h_corr_gen_corr_gen_ratio_cent%d",i);  
      h_corr_gen_corr_gen_ratio[i] = (TH1D*)h_corr_gen_ratio[i]->Clone(saythis);
      h_corr_gen_corr_gen_ratio[i]->Divide(h_corr_gen_ratio[0]);
      sprintf(saythis,"h_uncorr_uncorr_ratio_cent%d",i);  
      h_uncorr_uncorr_ratio[i] = (TH1D*)h_reco_full_n[i]->Clone(saythis);
      h_uncorr_uncorr_ratio[i]->Divide(h_reco_full_n[0]);
      sprintf(saythis,"h_corr_corr_ratio_cent%d",i);  
      h_corr_corr_ratio[i] = (TH1D*)h_reco_corr_n[i]->Clone(saythis);
      h_corr_corr_ratio[i]->Divide(h_reco_corr_n[0]);
      sprintf(saythis,"h_uncorr_uncorr_ratio_MC_cent%d",i);  
      h_uncorr_uncorr_ratio_MC[i] = (TH1D*)h_reco_full_MC_n[i]->Clone(saythis);
      h_uncorr_uncorr_ratio_MC[i]->Divide(h_reco_full_MC_n[0]);
      sprintf(saythis,"h_corr_corr_ratio_MC_cent%d",i);  
      h_corr_corr_ratio_MC[i] = (TH1D*)h_reco_corr_MC_n[i]->Clone(saythis);
      h_corr_corr_ratio_MC[i]->Divide(h_reco_corr_MC_n[0]);
      sprintf(saythis,"h_gen_gen_ratio_cent%d",i);  
      h_gen_gen_ratio[i] = (TH1D*)h_gen_full_n[i]->Clone(saythis);
      h_gen_gen_ratio[i]->Divide(h_gen_full_n[0]);
    }
  }

  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};
  TLine *tl1 = new TLine(100,1.,500,1.);
  tl1->SetLineStyle(2);

  TLegend *leg_reco = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco->SetLineColor(0);
  leg_reco->SetFillColor(0);
  leg_reco->AddEntry(h_reco_full_n[1], "L2L3 jets", "lepf");
  leg_reco->AddEntry(h_reco_corr_n[1], "L2L3+JFF jets", "lepf");

  TLegend *leg_reco1 = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco1->SetLineColor(0);
  leg_reco1->SetFillColor(0);
  leg_reco1->AddEntry((TObject*)0, "PbPb & P+H", "");

  TLegend *leg_reco3 = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco3->SetLineColor(0);
  leg_reco3->SetFillColor(0);
  leg_reco3->AddEntry((TObject*)0, "PbPb and pp data", "");

  TLegend *leg_reco2 = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco2->SetLineColor(0);
  leg_reco2->SetFillColor(0);
  leg_reco2->AddEntry((TObject*)0, "pp & Pythia", "");

  TLegend *leg_reco4 = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco4->SetLineColor(0);
  leg_reco4->SetFillColor(0);
  leg_reco4->AddEntry((TObject*)0, "pp data", "");

  TLegend *leg_ratio = new TLegend(0.5,0.8,0.99,0.99);
  leg_ratio->SetLineColor(0);
  leg_ratio->SetFillColor(0);
  leg_ratio->AddEntry(h_reco_ratio[1], "L2L3+JFF data / L2L3 data", "lepf");
  leg_ratio->AddEntry(h_reco_ratio_MC[1], "L2L3+JFF MC / L2L3 MC", "lepf");  

  TCanvas *c_reco_full = new TCanvas("c_reco_full","",1200,480);
  c_reco_full->Divide(5,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=0; i<nCBins; i++){
    if(i==0){
    c_reco_full->cd(i+1);
    h_reco_corr_n[i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_n[i]->GetYaxis()->SetTitle("");
    h_reco_corr_n[i]->SetLineColor(kRed);
    h_reco_corr_n[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_n[i]->Draw("e1 same");
    h_reco_full_n[i]->SetLineColor(kBlue);
    h_reco_full_n[i]->Draw("e1 same");
    leg_reco->Draw("same");
    leg_reco4->Draw("same");

    c_reco_full->cd(i+6);

    h_reco_ratio[i]->GetXaxis()->SetTitle("pT");
    h_reco_ratio[i]->GetYaxis()->SetTitle("ratio");
    h_reco_ratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_reco_ratio[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_ratio[i]->SetLineColor(kBlack);
    h_reco_ratio[i]->Draw("e1 same");
    h_reco_ratio_MC[i]->SetLineColor(6);
    h_reco_ratio_MC[i]->Draw("e1 same");
    leg_ratio->Draw("same");
    leg_reco2->Draw("same");
    tl1->Draw("same");
    }

    if(i>0){
    c_reco_full->cd(i+1);
    h_reco_corr_n[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_n[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr_n[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr_n[5-i]->SetLineColor(kRed);
    h_reco_corr_n[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_n[5-i]->Draw("e1 same");
    h_reco_full_n[5-i]->SetLineColor(kBlue);
    h_reco_full_n[5-i]->Draw("e1 same");
    leg_reco3->Draw("same");
   
    c_reco_full->cd(i+6);

    h_reco_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_reco_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_reco_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_ratio[5-i]->SetLineColor(kBlack);
    h_reco_ratio[5-i]->Draw("e1 same");
    h_reco_ratio_MC[5-i]->SetLineColor(6);
    h_reco_ratio_MC[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_reco1->Draw("same");
    //if(i==1) leg_ratio->Draw("same");
    }
  }  

//////////////////////////////////////////////////////////////////

//////////////// reco gen spectra and ratio //////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_reco_gen = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco_gen->SetLineColor(0);
  leg_reco_gen->SetFillColor(0);
  leg_reco_gen->AddEntry(h_reco_full_MC_n[1], "L2L3 jets", "lepf");
  leg_reco_gen->AddEntry(h_reco_corr_MC_n[1], "L2L3+JFF jets", "lepf");
  leg_reco_gen->AddEntry(h_gen_full_n[1], "gen jets", "lepf");
  //leg_reco->AddEntry((TObject*)0, "pp & Pythia", "");

  TLegend *leg_gen_ratio = new TLegend(0.5,0.8,0.99,0.99);
  leg_gen_ratio->SetLineColor(0);
  leg_gen_ratio->SetFillColor(0);
  leg_gen_ratio->AddEntry(h_corr_gen_ratio[1], "L2L3+JFF reco / gen", "lepf");
  leg_gen_ratio->AddEntry(h_uncorr_gen_ratio[1], "L2L3 reco / gen", "lepf");

  TCanvas *c_reco_gen = new TCanvas("c_reco_gen","",1200,480);
  c_reco_gen->Divide(5,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=0; i<nCBins; i++){
    if(i==0){
    c_reco_gen->cd(i+1);
    h_reco_corr_MC_n[i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC_n[i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC_n[i]->SetLineColor(kRed);
    h_reco_corr_MC_n[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC_n[i]->Draw("e1 same");
    h_reco_full_MC_n[i]->SetLineColor(kBlue);
    h_reco_full_MC_n[i]->Draw("e1 same");
    h_gen_full_n[i]->SetLineColor(kGreen);
    h_gen_full_n[i]->Draw("e1 same");
    leg_reco2->Draw("same");
    leg_reco_gen->Draw("same");

    c_reco_gen->cd(i+6);

    h_corr_gen_ratio[i]->GetXaxis()->SetTitle("pT");
    h_corr_gen_ratio[i]->GetYaxis()->SetTitle("ratio");
    h_corr_gen_ratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_corr_gen_ratio[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_gen_ratio[i]->SetLineColor(kBlack);
    h_corr_gen_ratio[i]->Draw("e1 same");
    h_uncorr_gen_ratio[i]->SetLineColor(6);
    h_uncorr_gen_ratio[i]->Draw("e1 same");
    leg_gen_ratio->Draw("same");
    tl1->Draw("same");
    }

    if(i>0){
    c_reco_gen->cd(i+1);
    h_reco_corr_MC_n[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC_n[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC_n[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr_MC_n[5-i]->SetLineColor(kRed);
    h_reco_corr_MC_n[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC_n[5-i]->Draw("e1 same");
    h_reco_full_MC_n[5-i]->SetLineColor(kBlue);
    h_reco_full_MC_n[5-i]->Draw("e1 same");
    h_gen_full_n[5-i]->SetLineColor(kGreen);
    h_gen_full_n[5-i]->Draw("e1 same");
    leg_reco1->Draw("same");

    c_reco_gen->cd(i+6);

    h_corr_gen_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_gen_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_gen_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_corr_gen_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_gen_ratio[5-i]->SetLineColor(kBlack);
    h_corr_gen_ratio[5-i]->Draw("e1 same");
    h_uncorr_gen_ratio[5-i]->SetLineColor(6);
    h_uncorr_gen_ratio[5-i]->Draw("e1 same");
    //if(i==1) leg_gen_ratio->Draw("same");
    tl1->Draw("same");
    }
  }

//////////////////////////////////////////////////////////////////

//////////////// PbPb - pp MC //////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_gen = new TLegend(0.5,0.8,0.99,0.99);
  leg_gen->SetLineColor(0);
  leg_gen->SetFillColor(0);
  leg_gen->AddEntry(h_gen_gen_ratio[1], "P+H gen / Pythia gen", "lepf");

  TLegend *leg_rec = new TLegend(0.5,0.8,0.99,0.99);
  leg_rec->SetLineColor(0);
  leg_rec->SetFillColor(0);
  leg_rec->AddEntry(h_corr_corr_ratio_MC[1], "P+H JFF-corr reco / Pythia JFF-corr reco", "lepf");
  leg_rec->AddEntry(h_uncorr_uncorr_ratio_MC[1], "P+H L2-L3 reco / Pythia L2-L3 reco", "lepf");

  TLegend *leg_data = new TLegend(0.5,0.8,0.99,0.99);
  leg_data->SetLineColor(0);
  leg_data->SetFillColor(0);
  leg_data->AddEntry(h_corr_corr_ratio[1], "PbPb JFF+L2L3 / pp JFF+L2L3", "lepf");
  leg_data->AddEntry(h_uncorr_uncorr_ratio[1], "PbPb L2L3 / pp L2L3", "lepf");

  TCanvas *c_gen_gen = new TCanvas("c_gen_gen","",960,720);
  c_gen_gen->Divide(4,3);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    if(i>0){
    c_gen_gen->cd(i);
    h_gen_gen_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_gen_gen_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_gen_gen_ratio[5-i]->SetTitle(centVars[4-i].c_str());
    h_gen_gen_ratio[5-i]->SetLineColor(kGreen);
    h_gen_gen_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_gen_gen_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_gen_gen_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_gen->Draw("same");

    c_gen_gen->cd(i+4);
    h_corr_corr_ratio_MC[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_corr_ratio_MC[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_corr_ratio_MC[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_corr_corr_ratio_MC[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_corr_ratio_MC[5-i]->SetLineColor(kRed);
    h_corr_corr_ratio_MC[5-i]->Draw("e1 same");
    h_uncorr_uncorr_ratio_MC[5-i]->SetLineColor(kBlue);
    h_uncorr_uncorr_ratio_MC[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_rec->Draw("same");

    c_gen_gen->cd(i+8);
    h_corr_corr_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_corr_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_corr_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_corr_corr_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_corr_ratio[5-i]->SetLineColor(kRed);
    h_corr_corr_ratio[5-i]->Draw("e1 same");
    h_uncorr_uncorr_ratio[5-i]->SetLineColor(kBlue);
    h_uncorr_uncorr_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_data->Draw("same");

    }
  }    
/*
//////////////////////////////////////////////////////////////////

//////////////////////// MC-data ratio ///////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_data_MC = new TLegend(0.5,0.8,0.99,0.99);
  leg_data_MC->SetLineColor(0);
  leg_data_MC->SetFillColor(0);
  leg_data_MC->AddEntry(h_reco_full_MC_n[1], "L2L3 jets", "lepf");
  leg_data_MC->AddEntry(h_reco_corr_MC_n[1], "L2L3+JFF jets", "lepf");
  //leg_reco->AddEntry((TObject*)0, "pp & Pythia", "");

  TLegend *leg_data_ratio = new TLegend(0.5,0.8,0.99,0.99);
  leg_data_ratio->SetLineColor(0);
  leg_data_ratio->SetFillColor(0);
  leg_data_ratio->AddEntry(h_dataMC_corr_ratio[1], "L2L3+JFF data / MC", "lepf");
  leg_data_ratio->AddEntry(h_dataMC_uncorr_ratio[1], "L2L3 data / MC", "lepf");

  TCanvas *c_data_MC = new TCanvas("c_data_MC","",1200,480);
  c_data_MC->Divide(5,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=0; i<nCBins; i++){
    if(i==0){
    c_data_MC->cd(i+1);
    h_reco_corr_MC_n[i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC_n[i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC_n[i]->SetLineColor(kRed);
    h_reco_corr_MC_n[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC_n[i]->Draw("e1 same");
    h_reco_full_MC_n[i]->SetLineColor(kBlue);
    h_reco_full_MC_n[i]->Draw("e1 same");
    h_reco_corr_n[i]->SetLineColor(kRed);
    h_reco_corr_n[i]->SetMarkerStyle(7);
    h_reco_corr_n[i]->Draw("e1 same");
    h_reco_full_n[i]->SetLineColor(kBlue);
    h_reco_full_n[i]->SetMarkerStyle(7);
    h_reco_full_n[i]->Draw("e1 same");
    leg_reco2->Draw("same");
    leg_data_MC->Draw("same");

    c_data_MC->cd(i+6);

    h_dataMC_corr_ratio[i]->GetXaxis()->SetTitle("pT");
    h_dataMC_corr_ratio[i]->GetYaxis()->SetTitle("ratio");
    h_dataMC_corr_ratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_dataMC_corr_ratio[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_dataMC_corr_ratio[i]->SetLineColor(kBlack);
    h_dataMC_corr_ratio[i]->Draw("e1 same");
    h_dataMC_uncorr_ratio[i]->SetLineColor(6);
    h_dataMC_uncorr_ratio[i]->Draw("e1 same");
    leg_data_ratio->Draw("same");
    tl1->Draw("same");
    }

    if(i>0){
    c_data_MC->cd(i+1);
    h_reco_corr_MC_n[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC_n[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC_n[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr_MC_n[5-i]->SetLineColor(kRed);
    h_reco_corr_MC_n[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC_n[5-i]->Draw("e1 same");
    h_reco_full_MC_n[5-i]->SetLineColor(kBlue);
    h_reco_full_MC_n[5-i]->Draw("e1 same");
    h_reco_corr_n[i]->SetLineColor(kRed);
    h_reco_corr_n[i]->SetMarkerStyle(7);
    h_reco_corr_n[i]->Draw("e1 same");
    h_reco_full_n[i]->SetLineColor(kBlue);
    h_reco_full_n[i]->SetMarkerStyle(7);
    h_reco_full_n[i]->Draw("e1 same");
    leg_reco1->Draw("same");

    c_data_MC->cd(i+6);

    h_dataMC_corr_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_dataMC_corr_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_dataMC_corr_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_dataMC_corr_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_dataMC_corr_ratio[5-i]->SetLineColor(kBlack);
    h_dataMC_corr_ratio[5-i]->Draw("e1 same");
    h_dataMC_uncorr_ratio[5-i]->SetLineColor(6);
    h_dataMC_uncorr_ratio[5-i]->Draw("e1 same");
    //if(i==1) leg_gen_ratio->Draw("same");
    tl1->Draw("same");
    }
  }
*/
//////////////////////////////////////////////////////////////////

//////////////////////// MC-data spectra ///////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_data_MC = new TLegend(0.5,0.8,0.99,0.99);
  leg_data_MC->SetLineColor(0);
  leg_data_MC->SetFillColor(0);
  leg_data_MC->AddEntry(h_reco_full_MC_un[1], "P+H L2L3 jets", "lepf");
  leg_data_MC->AddEntry(h_reco_full_un[1], "PbPb data L2L3 jets", "lepf");
  //leg_reco->AddEntry((TObject*)0, "pp & Pythia", "");

  TLegend *leg_data_ratio = new TLegend(0.5,0.8,0.99,0.99);
  leg_data_ratio->SetLineColor(0);
  leg_data_ratio->SetFillColor(0);
  leg_data_ratio->AddEntry(h_reco_corr_MC_un[1], "L2L3+JFF P+H", "lepf");
  leg_data_ratio->AddEntry(h_reco_corr_un[1], "L2L3+JFF PbPb data", "lepf");

  TCanvas *c_data_MC_spectra = new TCanvas("c_data_MC_spectra","",960,480);
  c_data_MC_spectra->Divide(4,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    c_data_MC_spectra->cd(i);
/*    h_reco_full_MC_un[5-i]->Rebin(10);
    h_reco_full_un[5-i]->Rebin(10);
    h_reco_full_MC_un[5-i]->Scale(1./10.);
    h_reco_full_un[5-i]->Scale(1./10.);*/    
    h_reco_full_MC_un[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_full_MC_un[5-i]->GetYaxis()->SetTitle("");
    h_reco_full_MC_un[5-i]->SetLineWidth(3);
    h_reco_full_MC_un[5-i]->SetLineColor(kCyan);
    h_reco_full_MC_un[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_full_MC_un[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_full_MC_un[5-i]->Draw("e1 same");
    h_reco_full_un[5-i]->SetLineColor(kBlack);
    h_reco_full_un[5-i]->Draw("e1 same");
    leg_reco2->Draw("same");
    leg_data_MC->Draw("same");

    c_data_MC_spectra->cd(i+4);
/*    h_reco_corr_MC_un[5-i]->Rebin(10);
    h_reco_corr_un[5-i]->Rebin(10);
    h_reco_corr_MC_un[5-i]->Scale(1./10.);
    h_reco_corr_un[5-i]->Scale(1./10.);*/
    h_reco_corr_MC_un[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC_un[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC_un[5-i]->SetLineWidth(3);
    h_reco_corr_MC_un[5-i]->SetLineColor(kCyan);
    h_reco_corr_MC_un[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC_un[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr_MC_un[5-i]->Draw("e1 same");
    h_reco_corr_un[5-i]->SetLineColor(kBlack);
    h_reco_corr_un[i]->Draw("e1 same");
    leg_reco2->Draw("same");
    leg_data_ratio->Draw("same");

  }

//////////////////////////////////////////////////////////////////

//////////////////////// (reco/gen) / (reco/gen) ///////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_recogen_recogen = new TLegend(0.5,0.8,0.99,0.99);
  leg_recogen_recogen->SetLineColor(0);
  leg_recogen_recogen->SetFillColor(0);
  leg_recogen_recogen->AddEntry(h_corr_gen_corr_gen_ratio[1], "(P+H L2L3+JFF/gen) / (Pythia L2L3+JFF/gen)", "lepf");
  leg_recogen_recogen->AddEntry(h_uncorr_gen_uncorr_gen_ratio[1], "(P+H L2L3/gen) / (Pythia L2L3/gen)", "lepf");
  //leg_reco->AddEntry((TObject*)0, "pp & Pythia", "");

  TCanvas *c_recogen_recogen = new TCanvas("c_recogen_recogen","",960,240);
  c_recogen_recogen->Divide(4,1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    c_recogen_recogen->cd(i);
    h_corr_gen_corr_gen_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_gen_corr_gen_ratio[5-i]->GetYaxis()->SetTitle("");
    h_corr_gen_corr_gen_ratio[5-i]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_corr_gen_corr_gen_ratio[5-i]->SetLineColor(kRed);
    h_corr_gen_corr_gen_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_gen_corr_gen_ratio[5-i]->SetTitle(centVars[4-i].c_str());
    h_corr_gen_corr_gen_ratio[5-i]->Draw("e1 same");
    h_uncorr_gen_uncorr_gen_ratio[5-i]->SetLineColor(kBlue);
    h_uncorr_gen_uncorr_gen_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_reco2->Draw("same");
    leg_recogen_recogen->Draw("same");
  }  

//////////////////////////////////////////////////////////////////

//////////////////////// Data/MC and RAA ///////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_data_MC_raa = new TLegend(0.5,0.8,0.99,0.99);
  leg_data_MC_raa->SetLineColor(0);
  leg_data_MC_raa->SetFillColor(0);
  leg_data_MC_raa->AddEntry(h_dataMC_uncorr_ratio[1], "L2L3 data / L2L3 MC", "lepf");
  leg_data_MC_raa->AddEntry(g_raa[1], "R_{AA}", "lp");  

  TCanvas *c_data_MC_raa = new TCanvas("c_data_MC_raa","",1200,300);
  c_data_MC_raa->Divide(4,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    c_data_MC_raa->cd(i);
    h_dataMC_corr_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_dataMC_corr_ratio[5-i]->GetYaxis()->SetTitle("");
    h_dataMC_corr_ratio[5-i]->GetYaxis()->SetRangeUser(0.,1.5);
    h_dataMC_corr_ratio[5-i]->SetTitle(centVars[4-i].c_str());
    h_dataMC_corr_ratio[5-i]->SetLineColor(kRed);
    h_dataMC_corr_ratio[5-i]->GetXaxis()->SetRangeUser(100.,280.);
    h_dataMC_corr_ratio[5-i]->Draw("e1 same");
    h_dataMC_corr_ratio[5-i]->SetLineColor(kRed);
    h_dataMC_corr_ratio[5-i]->Draw("e1 same");
    h_dataMC_uncorr_ratio[5-i]->SetLineColor(kBlue);
    h_dataMC_uncorr_ratio[5-i]->Draw("e1 same");
    g_raa[5-i]->SetLineColor(kBlue);
    g_raa[5-i]->Draw("e1 same");
    tl1->Draw("same");
     
    if(i==1){ 
      TLegend *leg_data_MC_raa = new TLegend(0.5,0.8,0.99,0.99);
      leg_data_MC_raa->SetLineColor(0);
      leg_data_MC_raa->SetFillColor(0);
      //leg_data_MC_raa->AddEntry(h_dataMC_corr_ratio_raascaled[4], "L2L3+JFF data / L2L3+JFF MC R_{AA} rescaled", "lepf");
      leg_data_MC_raa->AddEntry(h_dataMC_corr_ratio[4], "L2L3+JFF data / L2L3+JFF MC", "lepf");
      leg_data_MC_raa->AddEntry(h_dataMC_uncorr_ratio[4], "L2L3 data / L2L3 MC", "lepf");
      leg_data_MC_raa->AddEntry(g_raa[4], "R_{AA} ; 50-70%", "lp");
      leg_data_MC_raa->Draw("same");
    }
/*    
    if(i==2){ 
      TLegend *leg_data_MC_raa = new TLegend(0.5,0.8,0.99,0.99);
      leg_data_MC_raa->SetLineColor(0);
      leg_data_MC_raa->SetFillColor(0);
      //leg_data_MC_raa->AddEntry((TObject*)0, "P+H weighted by R_{AA}", "");
      leg_data_MC_raa->AddEntry(h_dataMC_uncorr_ratio_raascaled[3], "L2L3 data / L2L3 MC ; 30-50%", "lepf");
      leg_data_MC_raa->AddEntry(g_raa[3], "R_{AA} ; 30-50%", "lp");
      leg_data_MC_raa->Draw("same");
    }     

    if(i==3){ 
      TLegend *leg_data_MC_raa = new TLegend(0.5,0.8,0.99,0.99);
      leg_data_MC_raa->SetLineColor(0);
      leg_data_MC_raa->SetFillColor(0);
      leg_data_MC_raa->AddEntry(h_dataMC_uncorr_ratio_raascaled[2], "L2L3 data / L2L3 MC ; 10-30%", "lepf");
      leg_data_MC_raa->AddEntry(g_raa[2], "R_{AA} ; 10-30%", "lp");
      leg_data_MC_raa->Draw("same");
    }
    if(i==4){ 
      TLegend *leg_data_MC_raa = new TLegend(0.5,0.8,0.99,0.99);
      leg_data_MC_raa->SetLineColor(0);
      leg_data_MC_raa->SetFillColor(0);
      leg_data_MC_raa->AddEntry(h_dataMC_uncorr_ratio_raascaled[1], "L2L3 data / L2L3 MC ; 0-10%", "lepf");
      leg_data_MC_raa->AddEntry(g_raa[1], "R_{AA} ; 0-5%", "lp");
      leg_data_MC_raa->Draw("same");
    }
*/
  }


//////////////////////////////////////////////////////////////////

////////////////////// eta - phi MC /////////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_eta = new TLegend(0.5,0.8,0.99,0.99);
  leg_eta->SetLineColor(0);
  leg_eta->SetFillColor(0);
  leg_eta->AddEntry(h_eta_full_MC[1], "eta P+H reweighted", "lepf");
  leg_eta->AddEntry(h_eta_full_MC_nrw[1], "eta P+H not reweighted", "lepf");

  TLegend *leg_phi = new TLegend(0.5,0.8,0.99,0.99);
  leg_phi->SetLineColor(0);
  leg_phi->SetFillColor(0);
  leg_phi->AddEntry(h_phi_full_MC[1], "phi P+H reweighted", "lepf");
  leg_phi->AddEntry(h_phi_full_MC_nrw[1], "phi P+H not reweighted", "lepf");

  TCanvas *c_eta_phi = new TCanvas("c_eta_phi","",960,480);
  c_eta_phi->Divide(4,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    if(i>0){
    c_eta_phi->cd(i);
    h_eta_full_MC[5-i]->Scale(1./h_eta_full_MC[5-i]->Integral());
    h_eta_full_MC_nrw[5-i]->Scale(1./h_eta_full_MC_nrw[5-i]->Integral());
    h_eta_full_MC[5-i]->GetXaxis()->SetTitle("#eta");
    h_eta_full_MC[5-i]->GetYaxis()->SetTitle("");
    h_eta_full_MC[5-i]->SetTitle(centVars[4-i].c_str());
    h_eta_full_MC[5-i]->SetLineColor(kRed);
    h_eta_full_MC[5-i]->GetYaxis()->SetRangeUser(0,0.15);
    //h_eta_full_MC[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_eta_full_MC[5-i]->Draw("e1 same");
    h_eta_full_MC_nrw[5-i]->SetLineColor(kBlue);
    h_eta_full_MC_nrw[5-i]->Draw("e1 same");
    //tl1->Draw("same");
    if(i==1)leg_eta->Draw("same");

    c_eta_phi->cd(i+4);
    h_phi_full_MC[5-i]->Scale(1./h_phi_full_MC[5-i]->Integral());
    h_phi_full_MC_nrw[5-i]->Scale(1./h_phi_full_MC_nrw[5-i]->Integral());
    h_phi_full_MC[5-i]->GetXaxis()->SetTitle("#phi");
    h_phi_full_MC[5-i]->GetYaxis()->SetTitle("");
    h_phi_full_MC[5-i]->GetYaxis()->SetRangeUser(0,0.15);
    //h_phi_full_MC[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_phi_full_MC[5-i]->SetLineColor(kRed);
    h_phi_full_MC[5-i]->Draw("e1 same");
    h_phi_full_MC_nrw[5-i]->SetLineColor(kBlue);
    h_phi_full_MC_nrw[5-i]->Draw("e1 same");
    if(i==1)leg_phi->Draw("same");

    }
  }

//////////////////////////////////////////////////////////////////

////////////////////// vz MC /////////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_vz = new TLegend(0.5,0.8,0.99,0.99);
  leg_vz->SetLineColor(0);
  leg_vz->SetFillColor(0);
  leg_vz->AddEntry(h_vz_MC[1], "vz P+H reweighted", "lepf");
  leg_vz->AddEntry(h_vz_MC_nrw[1], "vz P+H not reweighted", "lepf");

  TCanvas *c_vz = new TCanvas("c_vz","",960,240);
  c_vz->Divide(4,1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    if(i>0){
    c_vz->cd(i);
    h_vz_MC[5-i]->Scale(1./h_vz_MC[5-i]->Integral());
    h_vz_MC_nrw[5-i]->Scale(1./h_vz_MC_nrw[5-i]->Integral());
    h_vz_MC[5-i]->GetXaxis()->SetTitle("vz");
    h_vz_MC[5-i]->GetYaxis()->SetTitle("");
    h_vz_MC[5-i]->SetTitle(centVars[4-i].c_str());
    h_vz_MC[5-i]->SetLineColor(kRed);
    h_vz_MC[5-i]->GetYaxis()->SetRangeUser(0,0.11);
    //h_vz_MC[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_vz_MC[5-i]->Draw("e1 same");
    h_vz_MC_nrw[5-i]->SetLineColor(kBlue);
    h_vz_MC_nrw[5-i]->Draw("e1 same");
    //tl1->Draw("same");
    if(i==1)leg_vz->Draw("same");

    }
  }

//////////////////////////////////////////////////////////////////

////////////////////// hibin MC /////////////////////////////

//////////////////////////////////////////////////////////////////

  TLegend *leg_hibin = new TLegend(0.5,0.8,0.99,0.99);
  leg_hibin->SetLineColor(0);
  leg_hibin->SetFillColor(0);
  leg_hibin->AddEntry(h_hibin, "hibin P+H reweighted", "lepf");
  leg_hibin->AddEntry(h_hibin_nrw, "hibin P+H not reweighted", "lepf");

  TCanvas *c_hibin = new TCanvas("c_hibin","",400,400);
  gStyle->SetOptStat(0);

    h_hibin->Scale(1./h_hibin->Integral());
    h_hibin_nrw->Scale(1./h_hibin_nrw->Integral());
    h_hibin->GetXaxis()->SetTitle("vz");
    h_hibin->GetYaxis()->SetTitle("");
    h_hibin->SetLineColor(kRed);
    h_hibin->GetYaxis()->SetRangeUser(0,0.04);
    h_hibin->Draw("HIST same");
    h_hibin_nrw->SetLineColor(kBlue);
    h_hibin_nrw->Draw("HIST same");
    leg_hibin->Draw("same");

}