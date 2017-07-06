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

#define nCBins 4
#define nptBins 55
#define npfbins 21

char saythis[500];

using namespace std;

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};

void draw_cymbal_drum(){

  TFile *closure_histos_cymbal = TFile::Open("/home/dhanush/Documents/JFF_corrections/closure_histos_Jul6_header_id145.root");
  TFile *closure_histos_drum = TFile::Open("/home/dhanush/Documents/JFF_corrections/closure_histos_Jun14_header_id145.root");

// defining histos

  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_corr[nCBins];
  TH1D *h_gen_full[nCBins];

  TH2F *h_jt_closure_ref_nocorr_cymbal[nCBins];
  TH2F *h_jt_closure_ref_ncs2_cymbal[nCBins];
  TH2F *h_jt_closure_q_nocorr_cymbal[nCBins];
  TH2F *h_jt_closure_q_ncs2_cymbal[nCBins];
  TH2F *h_jt_closure_g_nocorr_cymbal[nCBins];
  TH2F *h_jt_closure_g_ncs2_cymbal[nCBins];

  TProfile *h_jt_closure_ref_nocorr_cymbal_px[nCBins];
  TProfile *h_jt_closure_ref_ncs2_cymbal_px[nCBins];
  TProfile *h_jt_closure_q_nocorr_cymbal_px[nCBins];
  TProfile *h_jt_closure_q_ncs2_cymbal_px[nCBins];
  TProfile *h_jt_closure_g_nocorr_cymbal_px[nCBins];
  TProfile *h_jt_closure_g_ncs2_cymbal_px[nCBins];

  TH2F *h_jt_closure_ref_nocorr_drum[nCBins];
  TH2F *h_jt_closure_ref_ncs2_drum[nCBins];
  TH2F *h_jt_closure_q_nocorr_drum[nCBins];
  TH2F *h_jt_closure_q_ncs2_drum[nCBins];
  TH2F *h_jt_closure_g_nocorr_drum[nCBins];
  TH2F *h_jt_closure_g_ncs2_drum[nCBins];

  TProfile *h_jt_closure_ref_nocorr_drum_px[nCBins];
  TProfile *h_jt_closure_ref_ncs2_drum_px[nCBins];
  TProfile *h_jt_closure_q_nocorr_drum_px[nCBins];
  TProfile *h_jt_closure_q_ncs2_drum_px[nCBins];
  TProfile *h_jt_closure_g_nocorr_drum_px[nCBins];
  TProfile *h_jt_closure_g_ncs2_drum_px[nCBins];

  for(int ibin=0;ibin<nCBins;ibin++){

    h_jt_closure_ref_nocorr_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_ref_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_nocorr_cymbal_"+cent[ibin]));
    h_jt_closure_q_nocorr_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_q_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_nocorr_cymbal_"+cent[ibin]));
    h_jt_closure_g_nocorr_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_g_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_nocorr_cymbal_"+cent[ibin]));
    h_jt_closure_ref_ncs2_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_ref_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_ncs2_cymbal_"+cent[ibin]));
    h_jt_closure_q_ncs2_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_q_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_ncs2_cymbal_"+cent[ibin]));
    h_jt_closure_g_ncs2_cymbal[ibin] = (TH2F*)closure_histos_cymbal->Get((TString)("h_jt_closure_g_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_ncs2_cymbal_"+cent[ibin]));

    h_jt_closure_ref_nocorr_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_ref_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_nocorr_drum_"+cent[ibin]));
    h_jt_closure_q_nocorr_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_q_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_nocorr_drum_"+cent[ibin]));
    h_jt_closure_g_nocorr_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_g_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_nocorr_drum_"+cent[ibin]));
    h_jt_closure_ref_ncs2_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_ref_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_ncs2_drum_"+cent[ibin]));
    h_jt_closure_q_ncs2_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_q_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_ncs2_drum_"+cent[ibin]));
    h_jt_closure_g_ncs2_drum[ibin] = (TH2F*)closure_histos_drum->Get((TString)("h_jt_closure_g_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_ncs2_drum_"+cent[ibin]));    
/*    
    h_reco_full[ibin] = (TH1D*)closure_histos->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_"+cent[ibin]));
    h_reco_corr[ibin] = (TH1D*)closure_histos->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_"+cent[ibin]));
    h_gen_full[ibin] = (TH1D*)closure_histos->Get((TString)("h_gen_full_cent"+cent[ibin]))->Clone((TString)("h_gen_full_"+cent[ibin]));
*/
  }

  TLatex *l1[nCBins][nptBins];
  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};
  const string ptVars[55] = {"50<pT<60", "60<pT<70", "70<pT<80","80<pT<90", "90<pT<100", "100<pT<110", "110<pT<120","120<pT<130", "130<pT<140", "140<pT<150", "150<pT<160","160<pT<170", "170<pT<180", "180<pT<190", "190<pT<200", "200<pT<210","210<pT<220", "220<pT<230", "230<pT<240", "240<pT<250","250<pT<260", "260<pT<270", "270<pT<280", "280<pT<290","290<pT<300", "300<pT<310", "310<pT<320", "320<pT<330", "330<pT<340", "340<pT<350", "350<pT<360", "360<pT<370", "370<pT<380", "380<pT<390", "390<pT<400", "400<pT<410", "410<pT<420", "420<pT<430", "430<pT<440", "440<pT<450","450<pT<460", "460<pT<470", "470<pT<480", "480<pT<490", "490<pT<500", "500<pT<510", "510<pT<520", "520<pT<530", "530<pT<540","540<pT<550","550<pT<560", "560<pT<570", "570<pT<580", "580<pT<590","590<pT<600"};

  TLine *tl1 = new TLine(90,1.05,500,1.05);
  TLine *tl2 = new TLine(90,0.95,500,0.95);
  TLine *tl3 = new TLine(90,1.,500,1.);
  TLine *tl4 = new TLine(120,0.9.,120,1.1);
  TLine *tl5 = new TLine(120,0.95.,120,1.3);
  tl1->SetLineStyle(2);
  tl2->SetLineStyle(2);
  tl3->SetLineStyle(2);
  tl4->SetLineStyle(2);
  tl5->SetLineStyle(2);

  TCanvas *c_closure = new TCanvas("c_closure","",1200,300);
  c_closure->Divide(4,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c_closure->cd(i+1);
    
    h_jt_closure_ref_ncs2_cymbal_px[3-i] = h_jt_closure_ref_ncs2_cymbal[3-i]->ProfileX();
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->SetTitle(centVars[3-i].c_str());
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->SetLineColor(kBlack);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_ref_ncs2_cymbal_px[3-i]->Draw("e1 same");
   
    h_jt_closure_q_ncs2_cymbal_px[3-i] = h_jt_closure_q_ncs2_cymbal[3-i]->ProfileX();
    h_jt_closure_q_ncs2_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_q_ncs2_cymbal_px[3-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_ncs2_cymbal_px[3-i]->SetLineColor(kBlue);
    h_jt_closure_q_ncs2_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_q_ncs2_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_q_ncs2_cymbal_px[3-i]->Draw("e1 same");
    h_jt_closure_g_ncs2_cymbal_px[3-i] = h_jt_closure_g_ncs2_cymbal[3-i]->ProfileX();
    h_jt_closure_g_ncs2_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_g_ncs2_cymbal_px[3-i]->SetMarkerColor(kRed);
    h_jt_closure_g_ncs2_cymbal_px[3-i]->SetLineColor(kRed);
    h_jt_closure_g_ncs2_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_g_ncs2_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_g_ncs2_cymbal_px[3-i]->Draw("e1 same");

    h_jt_closure_ref_ncs2_drum_px[3-i] = h_jt_closure_ref_ncs2_drum[3-i]->ProfileX();
    h_jt_closure_ref_ncs2_drum_px[3-i]->Rebin(2);
    h_jt_closure_ref_ncs2_drum_px[3-i]->SetLineColor(kBlack);
    h_jt_closure_ref_ncs2_drum_px[3-i]->GetXaxis()->SetRangeUser(90.,500.);
    h_jt_closure_ref_ncs2_drum_px[3-i]->Draw("e0 same");
    h_jt_closure_q_ncs2_drum_px[3-i] = h_jt_closure_q_ncs2_drum[3-i]->ProfileX();
    h_jt_closure_q_ncs2_drum_px[3-i]->Rebin(2);
    h_jt_closure_q_ncs2_drum_px[3-i]->SetLineColor(kBlue);
    h_jt_closure_q_ncs2_drum_px[3-i]->Draw("e0 same");
    h_jt_closure_g_ncs2_drum_px[3-i] = h_jt_closure_g_ncs2_drum[3-i]->ProfileX();
    h_jt_closure_g_ncs2_drum_px[3-i]->Rebin(2);
    h_jt_closure_g_ncs2_drum_px[3-i]->SetLineColor(kRed);
    h_jt_closure_g_ncs2_drum_px[3-i]->Draw("e0 same");

    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");

    if(i==0){
      TLegend *leg3 = new TLegend(0.5,0.8,0.99,0.99);
      leg3->SetLineColor(0);
      leg3->SetFillColor(0);
      leg3->AddEntry((TObject*)0, "nCS ID=1,4,5", "");
      leg3->Draw("same");
    }

    if(i==1){
      TLegend *leg1 = new TLegend(0.5,0.8,0.99,0.99);
      leg1->SetLineColor(0);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_jt_closure_ref_ncs2_cymbal_px[3], "Cymbal", "lepf");
      leg1->AddEntry(h_jt_closure_ref_ncs2_drum_px[3], "Drum", "lepf");
      leg1->Draw("same");
    }

    if(i==2){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(0);
      leg2->SetFillColor(0);
      leg2->AddEntry(h_jt_closure_ref_ncs2_cymbal_px[2], "Corrected Incl. Jets", "lepf");
      leg2->AddEntry(h_jt_closure_q_ncs2_cymbal_px[2], "Corrected Quark Jets", "lepf");
      leg2->AddEntry(h_jt_closure_g_ncs2_cymbal_px[2], "Corrected Gluon Jets", "lepf");
      leg2->Draw("same");
    } 
  
  }

  TCanvas *c_closure_uncorr = new TCanvas("c_closure_uncorr","",1200,300);
  c_closure_uncorr->Divide(4,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c_closure_uncorr->cd(i+1);
    
    h_jt_closure_ref_nocorr_cymbal_px[3-i] = h_jt_closure_ref_nocorr_cymbal[3-i]->ProfileX();
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->SetTitle(centVars[3-i].c_str());
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_ref_nocorr_cymbal_px[3-i]->Draw("e1 same");
    h_jt_closure_q_nocorr_cymbal_px[3-i] = h_jt_closure_q_nocorr_cymbal[3-i]->ProfileX();
    h_jt_closure_q_nocorr_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_q_nocorr_cymbal_px[3-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_nocorr_cymbal_px[3-i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_q_nocorr_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_q_nocorr_cymbal_px[3-i]->Draw("e1 same");
    h_jt_closure_g_nocorr_cymbal_px[3-i] = h_jt_closure_g_nocorr_cymbal[3-i]->ProfileX();
    h_jt_closure_g_nocorr_cymbal_px[3-i]->Rebin(2);
    h_jt_closure_g_nocorr_cymbal_px[3-i]->SetMarkerColor(kRed);
    h_jt_closure_g_nocorr_cymbal_px[3-i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_cymbal_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_g_nocorr_cymbal_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_g_nocorr_cymbal_px[3-i]->Draw("e1 same");

    h_jt_closure_ref_nocorr_drum_px[3-i] = h_jt_closure_ref_nocorr_drum[3-i]->ProfileX();
    h_jt_closure_ref_nocorr_drum_px[3-i]->Rebin(2);
    h_jt_closure_ref_nocorr_drum_px[3-i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_drum_px[3-i]->GetXaxis()->SetRangeUser(90.,500.);
    h_jt_closure_ref_nocorr_drum_px[3-i]->Draw("e0 same");
    h_jt_closure_q_nocorr_drum_px[3-i] = h_jt_closure_q_nocorr_drum[3-i]->ProfileX();
    h_jt_closure_q_nocorr_drum_px[3-i]->Rebin(2);
    h_jt_closure_q_nocorr_drum_px[3-i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_drum_px[3-i]->Draw("e0 same");
    h_jt_closure_g_nocorr_drum_px[3-i] = h_jt_closure_g_nocorr_drum[3-i]->ProfileX();
    h_jt_closure_g_nocorr_drum_px[3-i]->Rebin(2);
    h_jt_closure_g_nocorr_drum_px[3-i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_drum_px[3-i]->Draw("e0 same");

    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");

    if(i==0){
      TLegend *leg3 = new TLegend(0.5,0.8,0.99,0.99);
      leg3->SetLineColor(0);
      leg3->SetFillColor(0);
      leg3->AddEntry((TObject*)0, "nCS ID=1,4,5", "");
      leg3->Draw("same");
    }

    if(i==1){
      TLegend *leg1 = new TLegend(0.5,0.8,0.99,0.99);
      leg1->SetLineColor(0);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_jt_closure_ref_nocorr_cymbal_px[3], "Cymbal", "lepf");
      leg1->AddEntry(h_jt_closure_ref_nocorr_drum_px[3], "Drum", "lepf");
      leg1->Draw("same");
    }

    if(i==2){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(0);
      leg2->SetFillColor(0);
      leg2->AddEntry(h_jt_closure_ref_nocorr_cymbal_px[1], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_jt_closure_q_nocorr_cymbal_px[1], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_jt_closure_g_nocorr_cymbal_px[1], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    } 
  
  }  

}