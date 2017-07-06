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

void test_draw(){

  //TFile *closure_histos_PbPb = TFile::Open("/home/dhanush/Documents/JEC/local/closure_histos_Apr18_header_final.root");	
  //TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_data_May22.root");
  TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/pptest_histos_data_Jun14.root");
  //TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/JEC/local/test_corr/pptest_histos_MC_Jun12.root");
  TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/JFF_corrections/pptest_histos_MC_Jun14.root"); 
  TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_data_Jun14.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("/home/dhanush/Documents/JFF_corrections/test_histos_MC_Jun14.root");

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

  TH1D *h_reco_ratio[nCBins];
  TH1D *h_reco_ratio_MC[nCBins];
  TH1D *h_corr_gen_ratio[nCBins];
  TH1D *h_uncorr_gen_ratio[nCBins];
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

  for(int ibin=0;ibin<5;ibin++){
    
    sprintf(saythis,"h_gen_full_n_cent%d",ibin);
    h_gen_full_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_gen_full_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_n_cent%d",ibin);
    h_reco_full_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_full_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_n_cent%d",ibin);
    h_reco_corr_n[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_corr_n[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_n_MC_cent%d",ibin);
    h_reco_full_n_MC[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_full_n_MC[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_n_MC_cent%d",ibin);
    h_reco_corr_n_MC[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_corr_n_MC[ibin]->Sumw2();

    if(ibin==0){
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_pp_data->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin])); 
      h_reco_full_un[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_un_"+cent[ibin]));
      h_reco_corr_un[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_un_"+cent[ibin]));    
      h_reco_full_MC_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_MC_un_"+cent[ibin]));      
      h_reco_corr_MC_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_MC_un_"+cent[ibin]));    
      if(!isdata) h_gen_full_un[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_cent"+cent[ibin]))->Clone((TString)("h_gen_full_un_"+cent[ibin]));
    }

    else if (ibin > 0){
    //else{
cout<<"here"<<endl;     
      h_ncs_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_data_"+cent[ibin])); 
      h_npf_pt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_npf_pt_cent"+cent[ibin-1]))->Clone((TString)("h_npf_pt_data_"+cent[ibin]));       
      h_ncs_pt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_ncs_pt_cent"+cent[ibin-1]))->Clone((TString)("h_ncs_pt_MC_"+cent[ibin]));
      h_reco_full_un[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_un_"+cent[ibin]));
      h_reco_corr_un[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_un_"+cent[ibin]));    
      h_reco_full_MC_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_full_cent"+cent[ibin-1]))->Clone((TString)("h_reco_full_MC_un_"+cent[ibin]));      
      h_reco_corr_MC_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]))->Clone((TString)("h_reco_corr_MC_un_"+cent[ibin]));    
      if(!isdata) h_gen_full_un[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_cent"+cent[ibin-1]))->Clone((TString)("h_gen_full_un_"+cent[ibin]));
    }
  if(ibin==3)h_reco_full_MC_un[3]->Draw();     
  }

  TLatex *l1[nCBins][nptBins];
  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};

  for(int i=0; i<nCBins; i++){

    for(int j=1;j<15;j++){
      //cout<<h_gen_full_un[i]->GetBinContent(j)<<endl;
      h_gen_full_n[i]->SetBinContent(j,h_gen_full_un[i]->GetBinContent(j)/h_gen_full_un[i]->GetBinWidth(j));
      h_gen_full_n[i]->SetBinError(j,h_gen_full_un[i]->GetBinError(j)/h_gen_full_un[i]->GetBinWidth(j));
      h_reco_full_n[i]->SetBinContent(j,h_reco_full_un[i]->GetBinContent(j)/h_reco_full_un[i]->GetBinWidth(j));
      h_reco_full_n[i]->SetBinError(j,h_reco_full_un[i]->GetBinError(j)/h_reco_full_un[i]->GetBinWidth(j));
      h_reco_corr_n[i]->SetBinContent(j,h_reco_corr_un[i]->GetBinContent(j)/h_reco_corr_un[i]->GetBinWidth(j));
      h_reco_corr_n[i]->SetBinError(j,h_reco_corr_un[i]->GetBinError(j)/h_reco_corr_un[i]->GetBinWidth(j));
      h_reco_full_n_MC[i]->SetBinContent(j,h_reco_full_MC_un[i]->GetBinContent(j)/h_reco_full_MC_un[i]->GetBinWidth(j));
      h_reco_full_n_MC[i]->SetBinError(j,h_reco_full_MC_un[i]->GetBinError(j)/h_reco_full_MC_un[i]->GetBinWidth(j));
      h_reco_corr_n_MC[i]->SetBinContent(j,h_reco_corr_MC_un[i]->GetBinContent(j)/h_reco_corr_MC_un[i]->GetBinWidth(j));
      h_reco_corr_n_MC[i]->SetBinError(j,h_reco_corr_MC_un[i]->GetBinError(j)/h_reco_corr_MC_un[i]->GetBinWidth(j));  
    }
  cout<<h_gen_full_un[i]->Integral()<<endl;

/*
  h_gen_full[i]->Scale(1./h_gen_full[i]->Integral());
  //h_reco_full[i]->Scale(1./h_reco_full[i]->Integral());
  //h_reco_corr[i]->Scale(1./h_reco_corr[i]->Integral());
  h_reco_full_MC[i]->Scale(1./h_reco_full_MC[i]->Integral());
  h_reco_corr_MC[i]->Scale(1./h_reco_corr_MC[i]->Integral());
*/
  h_reco_full_n[i]->Scale(1./h_reco_full_n[i]->GetBinContent(h_reco_full_n[i]->FindBin(120)));
  h_reco_corr_n[i]->Scale(1./h_reco_corr_n[i]->GetBinContent(h_reco_corr_n[i]->FindBin(120)));

/*
    h_reco_full[i]->Rebin(10);
    h_reco_full[i]->Scale(1./10.);
    //if(i>0) h_reco_full_all->Add(h_reco_full[i]);
    h_reco_corr[i]->Rebin(10);
    h_reco_corr[i]->Scale(1./10.);
    //if(i>0) h_reco_corr_all->Add(h_reco_corr[i]);
    h_reco_full_MC[i]->Rebin(10);    
    h_reco_full_MC[i]->Scale(1./10.);
    h_reco_corr_MC[i]->Rebin(10);
    h_reco_corr_MC[i]->Scale(1./10.);
    if(!isdata){
      h_gen_full[i]->Rebin(10);
      h_gen_full[i]->Scale(1./10.);
    } 
*/ 
  sprintf(saythis,"h_reco_ratio_cent%d",i);  
  h_reco_ratio[i] = (TH1D*)h_reco_corr_n[i]->Clone(saythis);
  h_reco_ratio[i]->Divide(h_reco_full_n[i]);
  sprintf(saythis,"h_reco_ratio_MC_cent%d",i);  
  h_reco_ratio_MC[i] = (TH1D*)h_reco_corr_n_MC[i]->Clone(saythis);
  h_reco_ratio_MC[i]->Divide(h_reco_full_n_MC[i]);
  sprintf(saythis,"h_corr_gen_ratio_cent%d",i);  
  h_corr_gen_ratio[i] = (TH1D*)h_reco_corr_MC[i]->Clone(saythis);
  h_corr_gen_ratio[i]->Divide(h_gen_full_n[i]);
  sprintf(saythis,"h_uncorr_gen_ratio_cent%d",i);  
  h_uncorr_gen_ratio[i] = (TH1D*)h_reco_full_n_MC[i]->Clone(saythis);
  h_uncorr_gen_ratio[i]->Divide(h_gen_full_n[i]);

  if(i>0) {
    sprintf(saythis,"h_uncorr_uncorr_ratio_cent%d",i);  
    h_uncorr_uncorr_ratio[i] = (TH1D*)h_reco_full_n[i]->Clone(saythis);
    h_uncorr_uncorr_ratio[i]->Divide(h_reco_full_n[0]);
    sprintf(saythis,"h_corr_corr_ratio_cent%d",i);  
    h_corr_corr_ratio[i] = (TH1D*)h_reco_corr[i]->Clone(saythis);
    h_corr_corr_ratio[i]->Divide(h_reco_corr[0]);
    sprintf(saythis,"h_uncorr_uncorr_ratio_MC_cent%d",i);  
    h_uncorr_uncorr_ratio_MC[i] = (TH1D*)h_reco_full_n_MC[i]->Clone(saythis);
    h_uncorr_uncorr_ratio_MC[i]->Divide(h_reco_full_n_MC[0]);
    sprintf(saythis,"h_corr_corr_ratio_MC_cent%d",i);  
    h_corr_corr_ratio_MC[i] = (TH1D*)h_reco_corr_MC[i]->Clone(saythis);
    h_corr_corr_ratio_MC[i]->Divide(h_reco_corr_MC[0]);
    sprintf(saythis,"h_gen_gen_ratio_cent%d",i);  
    h_gen_gen_ratio[i] = (TH1D*)h_gen_full_n[i]->Clone(saythis);
    h_gen_gen_ratio[i]->Divide(h_gen_full_n[0]);

  }

  //h_gen_full[i]->Scale(1./h_gen_full[i]->GetBinContent(h_gen_full[i]->FindBin(150)));
  //h_reco_full[i]->Scale(1./h_reco_full[i]->GetBinContent(h_reco_full[i]->FindBin(120)));
  //h_reco_corr[i]->Scale(1./h_reco_corr[i]->GetBinContent(h_reco_corr[i]->FindBin(120)));
  //h_reco_full_MC[i]->Scale(1./h_reco_full_MC[i]->GetBinContent(h_reco_full_MC[i]->FindBin(150)));
  //h_reco_corr_MC[i]->Scale(1./h_reco_corr_MC[i]->GetBinContent(h_reco_corr_MC[i]->FindBin(150)));

}

  //h_reco_full_all->Divide(h_reco_full[0]);
  //h_reco_corr_all->Divide(h_reco_corr[0]);
/*
  for(int i=0;i<45;i++){
    h_reco_corr_all->SetBinContent(i+1,(h_reco_corr_all->GetBinContent(i+1))/(h_reco_corr_all->GetEntries()));
    h_reco_corr[0]->SetBinContent(i+1,(h_reco_corr[0]->GetBinContent(i+1))/(h_reco_corr[0]->GetEntries()));  
  }
*/
  TLine *tl1 = new TLine(100,1.,500,1.);
  tl1->SetLineStyle(2);

  TCanvas *c_reco_ratio = new TCanvas("c_reco_ratio","",400,400);
  c_reco_ratio->cd(1);
  //h_corr_corr_ratio[4]->Draw();
  h_reco_full_n_MC[1]->Draw("same");
/*
  h_reco_corr_all->SetLineColor(kBlue); 
  h_reco_corr_all->Draw("same");
  h_reco_corr[0]->SetLineColor(kRed); 
  h_reco_corr[0]->Draw("same");
*/

  TLegend *leg_reco = new TLegend(0.5,0.8,0.99,0.99);
  leg_reco->SetLineColor(0);
  leg_reco->SetFillColor(0);
  leg_reco->AddEntry(h_reco_full[1], "L2L3 jets", "lepf");
  leg_reco->AddEntry(h_reco_corr[1], "L2L3+JFF jets", "lepf");

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
    leg_reco4->Draw("same");

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
    leg_reco2->Draw("same");
    tl1->Draw("same");
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
    leg_reco3->Draw("same");
   /*
    h_reco_corr_MC[5-i]->SetLineColor(kRed);
    h_reco_corr_MC[5-i]->SetLineStyle(2);
    h_reco_corr_MC[5-i]->Draw("e1 same");
    h_reco_full_MC[5-i]->SetLineColor(kBlue);
    h_reco_full_MC[5-i]->SetLineStyle(2);
    h_reco_full_MC[5-i]->Draw("e1 same");
    */
    //if(i==1)leg_reco->Draw("same");

    c_reco_full->cd(i+6);

    h_reco_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_reco_ratio[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
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
  leg_reco_gen->AddEntry(h_reco_full_MC[1], "L2L3 jets", "lepf");
  leg_reco_gen->AddEntry(h_reco_corr_MC[1], "L2L3+JFF jets", "lepf");
  leg_reco_gen->AddEntry(h_gen_full[1], "gen jets", "lepf");
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
    h_reco_corr_MC[i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC[i]->GetYaxis()->SetTitle("");
    //h_reco_corr[i]->SetTitle(centVars[i].c_str());
    h_reco_corr_MC[i]->SetLineColor(kRed);
    h_reco_corr_MC[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC[i]->Draw("e1 same");
    h_reco_full_MC[i]->SetLineColor(kBlue);
    h_reco_full_MC[i]->Draw("e1 same");
    h_gen_full[i]->SetLineColor(kGreen);
    h_gen_full[i]->Draw("e1 same");
    leg_reco2->Draw("same");
   /*
    h_reco_corr_MC[i]->SetLineColor(kRed);
    h_reco_corr_MC[i]->SetLineStyle(2);
    h_reco_corr_MC[i]->Draw("e1 same");
    h_reco_full_MC[i]->SetLineColor(kBlue);
    h_reco_full_MC[i]->SetLineStyle(2);
    h_reco_full_MC[i]->Draw("e1 same");
    */
    leg_reco_gen->Draw("same");

    c_reco_gen->cd(i+6);

    h_corr_gen_ratio[i]->GetXaxis()->SetTitle("pT");
    h_corr_gen_ratio[i]->GetYaxis()->SetTitle("ratio");
    h_corr_gen_ratio[i]->GetYaxis()->SetRangeUser(0.2,1.8);
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
    h_reco_corr_MC[5-i]->GetXaxis()->SetTitle("pT");
    h_reco_corr_MC[5-i]->GetYaxis()->SetTitle("");
    h_reco_corr_MC[5-i]->SetTitle(centVars[4-i].c_str());
    h_reco_corr_MC[5-i]->SetLineColor(kRed);
    h_reco_corr_MC[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_reco_corr_MC[5-i]->Draw("e1 same");
    h_reco_full_MC[5-i]->SetLineColor(kBlue);
    h_reco_full_MC[5-i]->Draw("e1 same");
    h_gen_full[5-i]->SetLineColor(kGreen);
    h_gen_full[5-i]->Draw("e1 same");
    leg_reco1->Draw("same");
   /*
    h_reco_corr_MC[5-i]->SetLineColor(kRed);
    h_reco_corr_MC[5-i]->SetLineStyle(2);
    h_reco_corr_MC[5-i]->Draw("e1 same");
    h_reco_full_MC[5-i]->SetLineColor(kBlue);
    h_reco_full_MC[5-i]->SetLineStyle(2);
    h_reco_full_MC[5-i]->Draw("e1 same");
    */
    //if(i==1)leg_reco_gen->Draw("same");

    c_reco_gen->cd(i+6);

    h_corr_gen_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_gen_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_gen_ratio[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
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
    h_gen_gen_ratio[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
    h_gen_gen_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_gen_gen_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_gen->Draw("same");

    c_gen_gen->cd(i+4);
    h_corr_corr_ratio_MC[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_corr_ratio_MC[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_corr_ratio_MC[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
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
    //h_corr_corr_ratio[5-i]->GetYaxis()->SetRangeUser(0.2,1.8);
    h_corr_corr_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_corr_ratio[5-i]->SetLineColor(kRed);
    h_corr_corr_ratio[5-i]->Draw("e1 same");
    h_uncorr_uncorr_ratio[5-i]->SetLineColor(kBlue);
    h_uncorr_uncorr_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    leg_data->Draw("same");

    }
  }

  TCanvas *c_PbPb_pp = new TCanvas("c_PbPb_pp","",1200,300);
  c_PbPb_pp->Divide(4,1);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  for(int i=1; i<nCBins; i++){

    if(i>0){
    c_PbPb_pp->cd(i);
    h_corr_corr_ratio[5-i]->GetXaxis()->SetTitle("pT");
    h_corr_corr_ratio[5-i]->GetYaxis()->SetTitle("ratio");
    h_corr_corr_ratio[5-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_corr_corr_ratio[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_corr_corr_ratio[5-i]->SetTitle(centVars[4-i].c_str());
    h_corr_corr_ratio[5-i]->SetLineColor(kRed);
    h_corr_corr_ratio[5-i]->Draw("e1 same");
    h_uncorr_uncorr_ratio[5-i]->SetLineColor(kBlue);
    h_uncorr_uncorr_ratio[5-i]->Draw("e1 same");
    tl1->Draw("same");
    if(i==1) leg_data->Draw("same");
    leg_reco3->Draw("same");

    }
  }


  TCanvas *c_ncs = new TCanvas("c_ncs","",1200,300);
  c_ncs->Divide(4,1);
  
  for(int i=1; i<nCBins; i++){
    if(i==0){
    c_ncs->cd(i+1);
    h_ncs_data[i] = h_ncs_pt_data[i]-> ProfileX();
    //h_ncs_data[i]->Scale(1./(h_ncs_data[i]->Integral()));
    h_ncs_MC[i] = h_ncs_pt_MC[i]-> ProfileX();
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
    h_ncs_data[5-i] = h_ncs_pt_data[5-i]-> ProjectionY();
    //h_npf_data[5-i] = h_npf_pt_data[5-i]-> ProjectionY();
    h_ncs_data[5-i]->Scale(1./(h_ncs_data[5-i]->Integral()));
    //h_ncs_data[5-i]->Rebin(10);
    //h_npf_data[5-i]->Rebin(10);
    //h_ncs_data[5-i]->Scale(1./10.);
    h_ncs_MC[5-i] = h_ncs_pt_MC[5-i]-> ProjectionY();
    h_ncs_MC[5-i]->Scale(1./(h_ncs_MC[5-i]->Integral()));
    //h_ncs_MC[5-i]->Rebin(10);
    //h_ncs_MC[5-i]->Scale(1./10.);
    h_ncs_data[5-i]->GetXaxis()->SetTitle("nCS cand ; pT>2");
    //if(i==1)h_ncs_data[5-i]->GetYaxis()->SetTitle("<ncand> ; pT>2");
    h_ncs_data[5-i]->SetTitle(centVars[4-i].c_str());
    h_ncs_data[5-i]->SetLineColor(kBlack);
    h_ncs_data[5-i]->GetXaxis()->SetRangeUser(0.,30.);
    h_ncs_data[5-i]->GetYaxis()->SetRangeUser(0.,.16);
    h_ncs_data[5-i]->Draw("e1 same");
    h_ncs_MC[5-i]->SetLineColor(kPink);
    h_ncs_MC[5-i]->Draw("e1 same");
    //h_npf_data[5-i]->SetLineColor(kGreen);
    //h_npf_data[5-i]->Draw("e1 same");
    if(i==1){
      TLegend *leg_ncs = new TLegend(0.5,0.8,0.99,0.99);
      leg_ncs->SetLineColor(0);
      //leg_ncs->AddEntry(h_npf_data[4], "data nPF id=1", "lepf");
      leg_ncs->AddEntry(h_ncs_data[4], "data nCS id=1,4,5", "lepf");
      leg_ncs->AddEntry(h_ncs_MC[4], "MC nCS id=1,4,5", "lepf");
      leg_ncs->Draw("same"); 
    }
    }
  }
}