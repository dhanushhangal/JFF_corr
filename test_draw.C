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

char saythis[500];

using namespace std;

bool isdata = true;

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};

void test_draw(){

  //TFile *closure_histos_PbPb = TFile::Open("/home/dhanush/Documents/JEC/local/closure_histos_Apr18_header_final.root");	
  TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_data_May22.root");
  TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_MC_May22.root");
  TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_data_Jun9.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_MC_Jun9.root");

  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_full_all;
  TH1D *h_reco_corr[nCBins];
  TH1D *h_reco_corr_all;
  TH1D *h_reco_full_MC[nCBins];
  TH1D *h_reco_corr_MC[nCBins];
  TH1D *h_gen_full[nCBins];
  TH1D *h_uncorr_ratio[nCBins];
  TH1D *h_corr_ratio[nCBins];
  TH1D *h_reco_ratio[nCBins];
  TH1D *h_reco_ratio_MC[nCBins];
  TH2F *h_ncs_pt_data[nCBins];
  TH2F *h_npf_pt_data[nCBins];
  TH2F *h_ncs_pt_MC[nCBins];
  TH1D *h_ncs_data[nCBins];
  TH1D *h_npf_data[nCBins];  
  TH1D *h_ncs_MC[nCBins];
  TH1D *h_ratio_PbPb_pp_data;

  //closure_histos_PbPb->cd();

  h_reco_full_all = new TH1D("h_reco_full_all","",45,50.,500.);
  h_reco_full_all->Sumw2();
  h_reco_corr_all = new TH1D("h_reco_corr_all","",45,50.,500.);
  h_reco_corr_all->Sumw2();

  for(int ibin=0;ibin<nCBins;ibin++){
/*
    sprintf(saythis,"h_uncorr_ratio_cent%d",ibin);
    h_uncorr_ratio[ibin] = new TH1D(saythis,"",50,0.,500.);
    h_uncorr_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_corr_ratio_cent%d",ibin);
    h_corr_ratio[ibin] = new TH1D(saythis,"",50,0.,500.);
    h_corr_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_reco_ratio_MC_cent%d",ibin);
    h_reco_ratio_MC[ibin] = new TH1D(saythis,"",50,0.,500.);
    h_reco_ratio_MC[ibin]->Sumw2();

    sprintf(saythis,"h_reco_ratio_cent%d",ibin);
    h_reco_ratio[ibin] = new TH1D(saythis,"",50,0.,500.);
    h_reco_ratio[ibin]->Sumw2();
*/
    if(ibin==0){
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_pp_data->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin])); 
      h_reco_full[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_"+cent[ibin]));
      h_reco_corr[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_"+cent[ibin]));    
      h_reco_full_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_"+cent[ibin]));
      h_reco_corr_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_"+cent[ibin]));    
      if(!isdata) h_gen_full[ibin] = (TH1D*)closure_histos_pp->Get((TString)("h_gen_full_cent"+cent[ibin]))->Clone((TString)("h_gen_full_"+cent[ibin]));
    }

    if(ibin > 0){ 
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_npf_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_npf_pt_cent"+cent[ibin-1]))->Clone((TString)("h_npf_pt_data_"+cent[ibin]));       
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin]));
      h_reco_full[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_"+cent[ibin]));
      h_reco_corr[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_"+cent[ibin]));    
      h_reco_full_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_"+cent[ibin]));
      h_reco_corr_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_"+cent[ibin]));    
      if(!isdata) h_gen_full[ibin] = (TH1D*)closure_histos_PbPb->Get((TString)("h_gen_full_cent"+cent[ibin-1]))->Clone((TString)("h_gen_full_"+cent[ibin]));
    }
  
  }

  TLatex *l1[nCBins][nptBins];
  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};

  for(int i=0; i<nCBins; i++){

    h_reco_full[i]->Rebin(10);
    h_reco_full[i]->Scale(1./10.);
    if(i>0) h_reco_full_all->Add(h_reco_full[i]);
    h_reco_corr[i]->Rebin(10);
    h_reco_corr[i]->Scale(1./10.);
    if(i>0) h_reco_corr_all->Add(h_reco_corr[i]);
    h_reco_full_MC[i]->Rebin(10);
    h_reco_full_MC[i]->Scale(1./10.);
    h_reco_corr_MC[i]->Rebin(10);
    h_reco_corr_MC[i]->Scale(1./10.);
    if(!isdata){
      h_gen_full[i]->Rebin(10);
      h_gen_full[i]->Scale(1./10.);
    } 
/*
  	for(int j=0; j<45; j++){
      if(!isdata){
    	  h_uncorr_ratio[i]->SetBinContent(j+6,(h_reco_full[i]->GetBinContent(j+6))/(h_gen_full[i]->GetBinContent(j+6)));
    	  h_uncorr_ratio[i]->SetBinError(j+6,sqrt((h_reco_full[i]->GetBinContent(j+6))/(h_gen_full[i]->GetBinContent(j+6)))*sqrt(pow((h_reco_full[i]->GetBinError(j+6))/(h_reco_full[i]->GetBinContent(j+6)),2)+pow((h_gen_full[i]->GetBinError(j+6))/(h_gen_full[i]->GetBinContent(j+6)),2)));
    	  h_corr_ratio[i]->SetBinContent(j+6,(h_reco_corr[i]->GetBinContent(j+6))/(h_gen_full[i]->GetBinContent(j+6)));
    	  h_corr_ratio[i]->SetBinError(j+6,sqrt((h_reco_corr[i]->GetBinContent(j+6))/(h_gen_full[i]->GetBinContent(j+6)))*sqrt(pow((h_reco_corr[i]->GetBinError(j+6))/(h_reco_corr[i]->GetBinContent(j+6)),2)+pow((h_gen_full[i]->GetBinError(j+6))/(h_gen_full[i]->GetBinContent(j+6)),2)));
      }
      h_reco_ratio[i]->SetBinContent(j+6,(h_reco_corr[i]->GetBinContent(j+6))/(h_reco_full[i]->GetBinContent(j+6)));
      h_reco_ratio[i]->SetBinError(j+6,sqrt((h_reco_corr[i]->GetBinContent(j+6))/(h_reco_full[i]->GetBinContent(j+6)))*sqrt(pow((h_reco_corr[i]->GetBinError(j+6))/(h_reco_corr[i]->GetBinContent(j+6)),2)+pow((h_reco_full[i]->GetBinError(j+6))/(h_reco_full[i]->GetBinContent(j+6)),2)));
      h_reco_ratio_MC[i]->SetBinContent(j+6,(h_reco_corr_MC[i]->GetBinContent(j+6))/(h_reco_full_MC[i]->GetBinContent(j+6)));
      h_reco_ratio_MC[i]->SetBinError(j+6,sqrt((h_reco_corr_MC[i]->GetBinContent(j+6))/(h_reco_full_MC[i]->GetBinContent(j+6)))*sqrt(pow((h_reco_corr_MC[i]->GetBinError(j+6))/(h_reco_corr_MC[i]->GetBinContent(j+6)),2)+pow((h_reco_full_MC[i]->GetBinError(j+6))/(h_reco_full_MC[i]->GetBinContent(j+6)),2)));    
    }
*/  
  sprintf(saythis,"h_reco_ratio_cent%d",i);  
  h_reco_ratio[i] = (TH1D*)h_reco_corr[i]->Clone(saythis);
  h_reco_ratio[i]->Divide(h_reco_full[i]);
  sprintf(saythis,"h_reco_ratio_MC_cent%d",i);  
  h_reco_ratio_MC[i] = (TH1D*)h_reco_corr_MC[i]->Clone(saythis);
  h_reco_ratio_MC[i]->Divide(h_reco_full_MC[i]);

  }

  //h_reco_full_all->Divide(h_reco_full[0]);
  //h_reco_corr_all->Divide(h_reco_corr[0]);
/*
  for(int i=0;i<45;i++){
    h_reco_corr_all->SetBinContent(i+1,(h_reco_corr_all->GetBinContent(i+1))/(h_reco_corr_all->GetEntries()));
    h_reco_corr[0]->SetBinContent(i+1,(h_reco_corr[0]->GetBinContent(i+1))/(h_reco_corr[0]->GetEntries()));  
  }
*/
  TCanvas *c_reco_ratio = new TCanvas("c_reco_ratio","",400,400);
  c_reco_ratio->cd(1);
  //h_reco_full_all->Draw("same");
  h_reco_corr_all->SetLineColor(kBlue); 
  h_reco_corr_all->Draw("same");
  h_reco_corr[0]->SetLineColor(kRed); 
  h_reco_corr[0]->Draw("same");

  TLegend *leg_reco = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco->SetLineColor(kWhite);
  leg_reco->AddEntry(h_reco_full[0], "uncorr reco jets", "lepf");
  leg_reco->AddEntry(h_reco_corr[0], "corr reco jets", "lepf");
  leg_reco->AddEntry((TObject*)0, "pp & Pythia", "");

  TLegend *leg_reco1 = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco1->SetLineColor(kWhite);
  leg_reco1->AddEntry((TObject*)0, "PbPb & P+H", "");

  TLegend *leg_ratio = new TLegend(0.5,0.8,0.99,0.99);
  leg_ratio->SetLineColor(kWhite);
  leg_ratio->AddEntry(h_reco_ratio[0], "corr data / uncorr data", "lepf");
  leg_ratio->AddEntry(h_reco_ratio_MC[0], "corr MC / uncorr MC", "lepf");

  TCanvas *c_reco_full = new TCanvas("c_reco_full","",1200,480);
  c_reco_full->Divide(5,2);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=0; i<nCBins; i++){
    if(i==0){
    c_reco_full->cd(i+1);
    h_reco_corr[i]->GetXaxis()->SetTitle("pT");
    h_reco_corr[i]->GetYaxis()->SetTitle("");
    //h_reco_corr[i]->SetTitle(centVars[i].c_str());
    h_reco_corr[i]->SetLineColor(kRed);
    h_reco_corr[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr[i]->Draw("e1 same");
    h_reco_full[i]->SetLineColor(kBlue);
    h_reco_full[i]->Draw("e1 same");
   /*
    h_reco_corr_MC[i]->SetLineColor(kRed);
    h_reco_corr_MC[i]->SetLineStyle(2);
    h_reco_corr_MC[i]->Draw("e1 same");
    h_reco_full_MC[i]->SetLineColor(kBlue);
    h_reco_full_MC[i]->SetLineStyle(2);
    h_reco_full_MC[i]->Draw("e1 same");
    */
    leg_reco->Draw("same");

    c_reco_full->cd(i+6);

    h_reco_ratio[i]->GetXaxis()->SetTitle("pT");
    h_reco_ratio[i]->GetYaxis()->SetTitle("ratio");
    h_reco_ratio[i]->GetYaxis()->SetRangeUser(0.2,1.8);
    h_reco_ratio[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_ratio[i]->SetLineColor(kBlack);
    h_reco_ratio[i]->Draw("e1 same");
    h_reco_ratio_MC[i]->SetLineColor(6);
    h_reco_ratio_MC[i]->Draw("e1 same");
    leg_ratio->Draw("same");
    }

    if(i>0){
    c_reco_full->cd(i+1);
    h_reco_corr[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr[5-i]->SetLineColor(kRed);
    h_reco_corr[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr[5-i]->Draw("e1 same");
    h_reco_full[5-i]->SetLineColor(kBlue);
    h_reco_full[5-i]->Draw("e1 same");
    leg_reco1->Draw("same");
   /*
    h_reco_corr_MC[5-i]->SetLineColor(kRed);
    h_reco_corr_MC[5-i]->SetLineStyle(2);
    h_reco_corr_MC[5-i]->Draw("e1 same");
    h_reco_full_MC[5-i]->SetLineColor(kBlue);
    h_reco_full_MC[5-i]->SetLineStyle(2);
    h_reco_full_MC[5-i]->Draw("e1 same");
    */
    //if(i==0)leg_reco->Draw("same");

    c_reco_full->cd(i+6);

    h_reco_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_reco_ratio[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
    h_reco_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_ratio[5-i]->SetLineColor(kBlack);
    h_reco_ratio[5-i]->Draw("e1 same");
    h_reco_ratio_MC[5-i]->SetLineColor(6);
    h_reco_ratio_MC[5-i]->Draw("e1 same");
    //if(i==0) leg_ratio->Draw("same");
    }
  } 
  
  TCanvas *c_ncs = new TCanvas("c_ncs","",1200,300);
  c_ncs->Divide(4,1);
  
  for(int i=1; i<nCBins; i++){
    if(i==0){
    c_ncs->cd(i+1);
    h_ncs_data[i] = h_ncs_pt_data[i]-> ProjectionY();
    //h_ncs_data[i]->Scale(1./(h_ncs_data[i]->Integral()));
    h_ncs_MC[i] = h_ncs_pt_MC[i]-> ProjectionY();
    //h_ncs_MC[i]->Scale(1./(h_ncs_MC[i]->Integral()));
    h_ncs_data[i]->GetXaxis()->SetTitle("nPF");
    //h_reco_corr[i]->GetYaxis()->SetTitle("");
    h_ncs_data[i]->SetLineColor(kBlack);
    h_ncs_data[i]->GetXaxis()->SetRangeUser(0.,25.);
    h_ncs_data[i]->Draw("e1 same");
    h_ncs_MC[i]->SetLineColor(kPink);
    h_ncs_MC[i]->Draw("e1 same");
    }

    if(i>0){
    c_ncs->cd(i);
    h_ncs_data[5-i] = h_ncs_pt_data[5-i]-> ProfileX();
    h_npf_data[5-i] = h_npf_pt_data[5-i]-> ProfileX();
    //h_ncs_data[5-i]->Scale(1./(h_ncs_data[5-i]->Integral()));
    h_ncs_data[5-i]->Rebin(10);
    h_npf_data[5-i]->Rebin(10);
    //h_ncs_data[5-i]->Scale(1./10.);
    h_ncs_MC[5-i] = h_ncs_pt_MC[5-i]-> ProfileX();
    //h_ncs_MC[5-i]->Scale(1./(h_ncs_MC[5-i]->Integral()));
    h_ncs_MC[5-i]->Rebin(10);
    //h_ncs_MC[5-i]->Scale(1./10.);
    h_ncs_data[5-i]->GetXaxis()->SetTitle("jtpt");
    if(i==1)h_ncs_data[5-i]->GetYaxis()->SetTitle("<ncand> ; pT>2");
    h_ncs_data[5-i]->SetTitle(centVars[4-i].c_str());
    h_ncs_data[5-i]->SetLineColor(kBlack);
    h_ncs_data[5-i]->GetXaxis()->SetRangeUser(50.,500.);
    h_ncs_data[5-i]->GetYaxis()->SetRangeUser(0.,20.);
    h_ncs_data[5-i]->Draw("e1 same");
    h_ncs_MC[5-i]->SetLineColor(kPink);
    h_ncs_MC[5-i]->Draw("e1 same");
    h_npf_data[5-i]->SetLineColor(kGreen);
    h_npf_data[5-i]->Draw("e1 same");
    if(i==1){
      TLegend *leg_ncs = new TLegend(0.5,0.8,0.99,0.99);
      leg_ncs->SetLineColor(kWhite);
      leg_ncs->AddEntry(h_npf_data[4], "data nPF id=1", "lepf");
      leg_ncs->AddEntry(h_ncs_data[4], "data nCS id=1,4,5", "lepf");
      leg_ncs->AddEntry(h_ncs_MC[4], "MC nCS id=1,4,5", "lepf");
      leg_ncs->Draw("same"); 
    }
    }
  }
}