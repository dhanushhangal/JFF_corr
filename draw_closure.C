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

#define nCBins 5
#define nptBins 40
#define npfbins 21

char saythis[500];

using namespace std;

TString cent[5] = {"0","1","2","3","4"};
TString pt[41] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40"};

int jt_nbins = 41;
Double_t jt_bin_bounds[41] = {100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.};

void draw_closure(){

  TFile *closure_histos = TFile::Open("/home/dhanush/Documents/JFF_corrections/closure_histos_Jul31_drum_header_id145_rebin.root");
  TFile *closure_histos_pp = TFile::Open("/home/dhanush/Documents/JFF_corrections/ppclosure_histos_Jul28_header.root");

// defining histos

  TH1D *h_reco[nCBins][nptBins];
  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_corr[nCBins];
  TH1D *h_gen_full[nCBins];
  TH2F *h_ncs_2[nCBins][nptBins];
  TH2F *h_ncs_2_q[nCBins][nptBins];
  TH2F *h_ncs_2_g[nCBins][nptBins];
  TH2F *h_ncs_2_corr[nCBins][nptBins];
  TH2F *h_ncs_2_corr_q[nCBins][nptBins];
  TH2F *h_ncs_2_corr_g[nCBins][nptBins];
  TH2F *h_jt_closure_ref_nocorr[nCBins];
  TH2F *h_jt_closure_ref_ncs2[nCBins];
  TH2F *h_jt_closure_reco_nocorr[nCBins];
  TH2F *h_jt_closure_reco_ncs2[nCBins];
  TH2F *h_jt_closure_q_nocorr[nCBins];
  TH2F *h_jt_closure_q_ncs2[nCBins];
  TH2F *h_jt_closure_g_nocorr[nCBins];
  TH2F *h_jt_closure_g_ncs2[nCBins];
  TH2F *h_eta_closure_nocorr[nCBins];
  TH2F *h_eta_closure_ncs2[nCBins];
  TH2F *h_eta_closure_q_nocorr[nCBins];
  TH2F *h_eta_closure_q_ncs2[nCBins];
  TH2F *h_eta_closure_g_nocorr[nCBins];
  TH2F *h_eta_closure_g_ncs2[nCBins];

  TProfile *h_ncs_2_closure[nCBins][nptBins];
  TProfile *h_ncs_2_closure_corr[nCBins][nptBins];

  TH1D *h_ncs_2_dist[nCBins][nptBins];
  TH1D *h_ncs_2_dist_q[nCBins][nptBins];
  TH1D *h_ncs_2_dist_g[nCBins][nptBins];

  TH1D *h_ncs_2_res[nCBins][nptBins];
  TH1D *h_ncs_2_res_q[nCBins][nptBins];
  TH1D *h_ncs_2_res_g[nCBins][nptBins];

  TH1D *h_ncs_2_res_corr[nCBins][nptBins];
  TH1D *h_ncs_2_res_corr_q[nCBins][nptBins];
  TH1D *h_ncs_2_res_corr_g[nCBins][nptBins];

  TH1D *h_res[nCBins];
  TH1D *h_res_corr[nCBins];

  TH1D *h_mean[nCBins];
  TH1D *h_jt_closure_ref_nocorr_ppx[nCBins];
  TH1D *h_jt_closure_px[nCBins][nptBins];

  TProfile *h_ncs_closure_ref_nocorr_px[nCBins];
  TProfile *h_ncs_closure_ref_ncs2_px[nCBins];
  TProfile *h_jt_closure_ref_nocorr_px[nCBins];
  TProfile *h_jt_closure_ref_ncs2_px[nCBins];
  TProfile *h_jt_closure_reco_nocorr_px[nCBins];
  TProfile *h_jt_closure_reco_ncs2_px[nCBins];
  TProfile *h_jt_closure_q_nocorr_px[nCBins];
  TProfile *h_jt_closure_q_ncs2_px[nCBins];
  TProfile *h_jt_closure_g_nocorr_px[nCBins];
  TProfile *h_jt_closure_g_ncs2_px[nCBins];
  TProfile *h_eta_closure_nocorr_px[nCBins];
  TProfile *h_eta_closure_ncs2_px[nCBins];
  TProfile *h_eta_closure_q_nocorr_px[nCBins];
  TProfile *h_eta_closure_q_ncs2_px[nCBins];
  TProfile *h_eta_closure_g_nocorr_px[nCBins];
  TProfile *h_eta_closure_g_ncs2_px[nCBins];

  TF1 *f_ncs_1[nCBins][nptBins];
  TF1 *f_ncs_2[nCBins][nptBins];

  TF1 *f_closure[nCBins];
  TF1 *f_closure_2[nCBins];

  TF1 *f_ncs_res[nCBins][nptBins];
  TF1 *f_ncs_res_corr[nCBins][nptBins];

  TF1 *f_ncs_res_2[nCBins][nptBins];
  TF1 *f_ncs_res_corr_2[nCBins][nptBins];

  TGraphErrors *gr_p1_param_1[nCBins];
  TGraphErrors *gr_p0_param_1[nCBins];

  TGraphErrors *gr_p1_param_2[nCBins];
  TGraphErrors *gr_p0_param_2[nCBins];
  TGraphErrors *gr_ax_2[nCBins];
  int total; 

  for(int ibin=0;ibin<nCBins;ibin++){

    sprintf(saythis,"h_res_cent%d",ibin);
    h_res[ibin] = new TH1D("","",40,100.,500.);
    h_res[ibin]->Sumw2();

    sprintf(saythis,"h_res_corr_cent%d",ibin);
    h_res_corr[ibin] = new TH1D("","",40,100.,500.);
    h_res_corr[ibin]->Sumw2();

    sprintf(saythis,"h_mean_cent%d",ibin);
    h_mean[ibin] = new TH1D("","",45,50.,500.);
    h_mean[ibin]->Sumw2();

  }

  h_jt_closure_ref_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_ref_nocorr_cent"+cent[0]))->Clone((TString)("h_jt_closure_ref_nocorr_"+cent[0]));
  h_jt_closure_q_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_q_nocorr_cent"+cent[0]))->Clone((TString)("h_jt_closure_q_nocorr_"+cent[0]));
  h_jt_closure_g_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_g_nocorr_cent"+cent[0]))->Clone((TString)("h_jt_closure_g_nocorr_"+cent[0]));
  h_jt_closure_ref_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_ref_ncs2_cent"+cent[0]))->Clone((TString)("h_jt_closure_ref_ncs2_"+cent[0]));
  h_jt_closure_q_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_q_ncs2_cent"+cent[0]))->Clone((TString)("h_jt_closure_q_ncs2_"+cent[0]));
  h_jt_closure_g_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_g_ncs2_cent"+cent[0]))->Clone((TString)("h_jt_closure_g_ncs2_"+cent[0]));
  h_jt_closure_reco_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_reco_nocorr_cent"+cent[0]))->Clone((TString)("h_jt_closure_reco_nocorr_"+cent[0]));
  h_jt_closure_reco_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_jt_closure_reco_ncs2_cent"+cent[0]))->Clone((TString)("h_jt_closure_reco_ncs2_"+cent[0]));
  h_eta_closure_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_nocorr_cent"+cent[0]))->Clone((TString)("h_eta_closure_nocorr_"+cent[0]));
  h_eta_closure_q_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_q_nocorr_cent"+cent[0]))->Clone((TString)("h_eta_closure_q_nocorr_"+cent[0]));
  h_eta_closure_g_nocorr[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_g_nocorr_cent"+cent[0]))->Clone((TString)("h_eta_closure_g_nocorr_"+cent[0]));
  h_eta_closure_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_ncs2_cent"+cent[0]))->Clone((TString)("h_eta_closure_ncs2_"+cent[0]));
  h_eta_closure_q_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_q_ncs2_cent"+cent[0]))->Clone((TString)("h_eta_closure_q_ncs2_"+cent[0]));
  h_eta_closure_g_ncs2[0] = (TH2F*)closure_histos_pp->Get((TString)("h_eta_closure_g_ncs2_cent"+cent[0]))->Clone((TString)("h_eta_closure_g_ncs2_"+cent[0]));

  for(int ibin=1;ibin<nCBins;ibin++){

    h_jt_closure_ref_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_ref_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_ref_nocorr_"+cent[ibin]));
    h_jt_closure_q_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_q_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_q_nocorr_"+cent[ibin]));
    h_jt_closure_g_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_g_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_g_nocorr_"+cent[ibin]));
    h_jt_closure_ref_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_ref_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_ref_ncs2_"+cent[ibin]));
    h_jt_closure_q_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_q_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_q_ncs2_"+cent[ibin]));
    h_jt_closure_g_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_g_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_g_ncs2_"+cent[ibin]));
    h_jt_closure_reco_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_reco_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_reco_nocorr_"+cent[ibin]));
    h_jt_closure_reco_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_reco_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_jt_closure_reco_ncs2_"+cent[ibin]));
    h_eta_closure_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_nocorr_"+cent[ibin]));
    h_eta_closure_q_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_q_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_q_nocorr_"+cent[ibin]));
    h_eta_closure_g_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_g_nocorr_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_g_nocorr_"+cent[ibin]));
    h_eta_closure_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_ncs2_"+cent[ibin]));
    h_eta_closure_q_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_q_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_q_ncs2_"+cent[ibin]));
    h_eta_closure_g_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_eta_closure_g_ncs2_cent"+cent[ibin-1]))->Clone((TString)("h_eta_closure_g_ncs2_"+cent[ibin]));

  }

  for(int ibin3=0;ibin3<nptBins;ibin3++){
      
    h_ncs_2[0][ibin3] = (TH2F*)closure_histos_pp->Get((TString)("h_ncs2_nocorr_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_"+cent[0]+"_"+pt[ibin3]));
    h_ncs_2_corr[0][ibin3] = (TH2F*)closure_histos_pp->Get((TString)("h_ncs2_corr_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_corr_"+cent[0]+"_"+pt[ibin3]));    
    h_ncs_2_q[0][ibin3] = (TH2F*)closure_histos_pp->Get((TString)("h_ncs2_nocorr_q_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_q_"+cent[0]+"_"+pt[ibin3]));
    h_ncs_2_g[0][ibin3] = (TH2F*)closure_histos_pp->Get((TString)("h_ncs2_nocorr_g_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_g_"+cent[0]+"_"+pt[ibin3]));    
    h_ncs_2_res[0][ibin3] = (TH1D*)closure_histos_pp->Get((TString)("h_closure_nocorr_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_nocorr_res_"+cent[0]+"_"+pt[ibin3]));
    h_ncs_2_res_corr[0][ibin3] = (TH1D*)closure_histos_pp->Get((TString)("h_closure_ncs2_cent"+cent[0]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_ncs2_res_"+cent[0]+"_"+pt[ibin3]));

    h_ncs_2_closure[0][ibin3] = h_ncs_2[0][ibin3]->ProfileX();
    h_ncs_2_closure_corr[0][ibin3] = h_ncs_2_corr[0][ibin3]->ProfileX();
    
    h_ncs_2_dist[0][ibin3] = h_ncs_2[0][ibin3]->ProjectionX();
    h_ncs_2_dist_q[0][ibin3] = h_ncs_2_q[0][ibin3]->ProjectionX();
    h_ncs_2_dist_q[0][ibin3] -> Scale (1./(h_ncs_2_dist_q[0][ibin3])->Integral());
    h_ncs_2_dist_q[0][ibin3] -> SetLineColor (kBlue);
    h_ncs_2_dist_g[0][ibin3] = h_ncs_2_g[0][ibin3]->ProjectionX();
    h_ncs_2_dist_g[0][ibin3] -> Scale (1./(h_ncs_2_dist_g[0][ibin3])->Integral());
    h_ncs_2_dist_g[0][ibin3] -> SetLineColor (kRed);

    for(int ibin=1;ibin<nCBins;ibin++){
      h_ncs_2[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_corr[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_corr_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_corr_"+cent[ibin]+"_"+pt[ibin3]));    
      h_ncs_2_q[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_q_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_q_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_g[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_g_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_g_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_res[ibin][ibin3] = (TH1D*)closure_histos->Get((TString)("h_closure_nocorr_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_nocorr_res_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_res_corr[ibin][ibin3] = (TH1D*)closure_histos->Get((TString)("h_closure_ncs2_cent"+cent[ibin-1]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_ncs2_res_"+cent[ibin]+"_"+pt[ibin3]));
    
      h_ncs_2_closure[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProfileX();
      h_ncs_2_closure_corr[ibin][ibin3] = h_ncs_2_corr[ibin][ibin3]->ProfileX();
 
      h_ncs_2_dist[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProjectionX();
      h_ncs_2_dist_q[ibin][ibin3] = h_ncs_2_q[ibin][ibin3]->ProjectionX();
      h_ncs_2_dist_q[ibin][ibin3] -> Scale (1./(h_ncs_2_dist_q[ibin][ibin3])->Integral());
      h_ncs_2_dist_q[ibin][ibin3] -> SetLineColor (kBlue);
      h_ncs_2_dist_g[ibin][ibin3] = h_ncs_2_g[ibin][ibin3]->ProjectionX();
      h_ncs_2_dist_g[ibin][ibin3] -> Scale (1./(h_ncs_2_dist_g[ibin][ibin3])->Integral());
      h_ncs_2_dist_g[ibin][ibin3] -> SetLineColor (kRed);  

      sprintf(saythis,"f_ncs_res_cent%d_pt%d",ibin,ibin3);
      f_ncs_res[ibin][ibin3] = new TF1(saythis, "[0]*exp(-((x-[1])^2)/(2*([2]^2)))", 0., 2.);
      f_ncs_res[ibin][ibin3]->SetParameter(0,0.1);
      f_ncs_res[ibin][ibin3]->SetParameter(1,1.);
      f_ncs_res[ibin][ibin3]->SetParameter(2,0.12);

      sprintf(saythis,"f_ncs_res_2_cent%d_pt%d",ibin,ibin3);
      f_ncs_res_2[ibin][ibin3] = new TF1(saythis, "[0]*exp(-((x-[1])^2)/(2*([2]^2)))", 0., 2.);
      f_ncs_res_2[ibin][ibin3]->SetParameter(0,0.1);
      f_ncs_res_2[ibin][ibin3]->SetParameter(1,1.);
      f_ncs_res_2[ibin][ibin3]->SetParameter(2,0.12);
      f_ncs_res_2[ibin][ibin3]->SetLineColor(kBlue);

      sprintf(saythis,"f_ncs_res_corr_cent%d_pt%d",ibin,ibin3);
      f_ncs_res_corr[ibin][ibin3] = new TF1(saythis, "[0]*exp(-((x-[1])^2)/(2*([2]^2)))", 0., 2.);
      f_ncs_res_corr[ibin][ibin3]->SetParameter(0,0.1);
      f_ncs_res_corr[ibin][ibin3]->SetParameter(1,1.);
      f_ncs_res_corr[ibin][ibin3]->SetParameter(2,0.12);

      sprintf(saythis,"f_ncs_res_corr_2_cent%d_pt%d",ibin,ibin3);
      f_ncs_res_corr_2[ibin][ibin3] = new TF1(saythis, "[0]*exp(-((x-[1])^2)/(2*([2]^2)))", 0., 2.);
      f_ncs_res_corr_2[ibin][ibin3]->SetParameter(0,0.1);
      f_ncs_res_corr_2[ibin][ibin3]->SetParameter(1,1.);
      f_ncs_res_corr_2[ibin][ibin3]->SetParameter(2,0.12);
      f_ncs_res_corr_2[ibin][ibin3]->SetLineColor(kRed);

      if(ibin==3&&ibin3==0){
        f_ncs_res[ibin][ibin3]->SetParameter(0,0.12);
        f_ncs_res[ibin][ibin3]->SetParameter(1,1.1);
        f_ncs_res[ibin][ibin3]->SetParameter(2,0.17);
        f_ncs_res_2[ibin][ibin3]->SetParameter(0,0.12);
        f_ncs_res_2[ibin][ibin3]->SetParameter(1,1.1);
        f_ncs_res_2[ibin][ibin3]->SetParameter(2,0.17);
      }  
    }
  
  }

  Double_t par_res[nCBins][nptBins];
  Double_t par_res_err[nCBins][nptBins];

  Double_t par_res_corr[nCBins][nptBins];
  Double_t par_res_corr_err[nCBins][nptBins];

  for(int ibin=0;ibin<nCBins;ibin++){

    for(int ibin3=0;ibin3<nptBins;ibin3++){ 

      //h_ncs_2_res[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProjectionY();
      h_ncs_2_res[ibin][ibin3] -> Rebin(5);
      h_ncs_2_res[ibin][ibin3] -> Scale(1./5.);
      Double_t integral_nocorr = h_ncs_2_res[ibin][ibin3]->Integral();
      h_ncs_2_res[ibin][ibin3]->Scale(1./integral_nocorr);

      //h_ncs_2_res_corr[ibin][ibin3] = h_ncs_2_corr[ibin][ibin3]->ProjectionY();
      h_ncs_2_res_corr[ibin][ibin3] -> Rebin(5);
      h_ncs_2_res_corr[ibin][ibin3] -> Scale(1./5.);
      Double_t integral_corr = h_ncs_2_res_corr[ibin][ibin3]->Integral();
      h_ncs_2_res_corr[ibin][ibin3]->Scale(1./integral_corr);
    }
  }

  TLatex *l1[nCBins][nptBins];
  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};
  //const string ptVars[55] = {"50<pT<60", "60<pT<70", "70<pT<80","80<pT<90", "90<pT<100", "100<pT<110", "110<pT<120","120<pT<130", "130<pT<140", "140<pT<150", "150<pT<160","160<pT<170", "170<pT<180", "180<pT<190", "190<pT<200", "200<pT<210","210<pT<220", "220<pT<230", "230<pT<240", "240<pT<250","250<pT<260", "260<pT<270", "270<pT<280", "280<pT<290","290<pT<300", "300<pT<310", "310<pT<320", "320<pT<330", "330<pT<340", "340<pT<350", "350<pT<360", "360<pT<370", "370<pT<380", "380<pT<390", "390<pT<400", "400<pT<410", "410<pT<420", "420<pT<430", "430<pT<440", "440<pT<450","450<pT<460", "460<pT<470", "470<pT<480", "480<pT<490", "490<pT<500", "500<pT<510", "510<pT<520", "520<pT<530", "530<pT<540","540<pT<550","550<pT<560", "560<pT<570", "570<pT<580", "580<pT<590","590<pT<600"};
  const string ptVars[40] = {"100<pT<110", "110<pT<120","120<pT<130", "130<pT<140", "140<pT<150", "150<pT<160","160<pT<170", "170<pT<180", "180<pT<190", "190<pT<200", "200<pT<210","210<pT<220", "220<pT<230", "230<pT<240", "240<pT<250","250<pT<260", "260<pT<270", "270<pT<280", "280<pT<290","290<pT<300", "300<pT<310", "310<pT<320", "320<pT<330", "330<pT<340", "340<pT<350", "350<pT<360", "360<pT<370", "370<pT<380", "380<pT<390", "390<pT<400", "400<pT<410", "410<pT<420", "420<pT<430", "430<pT<440", "440<pT<450","450<pT<460", "460<pT<470", "470<pT<480", "480<pT<490", "490<pT<500"};
 
  TLine *tl1 = new TLine(100,1.05,500,1.05);
  TLine *tl2 = new TLine(100,0.95,500,0.95);
  TLine *tl3 = new TLine(100,1.,500,1.);
  TLine *tl4 = new TLine(120,0.9,120,1.1);
  TLine *tl5 = new TLine(120,1.,120,1.3);
  TLine *tl9 = new TLine(0,1.,30,1.);
  tl1->SetLineStyle(2);
  tl2->SetLineStyle(2);
  tl3->SetLineStyle(2);
  tl4->SetLineStyle(2);
  tl5->SetLineStyle(2);
  tl9->SetLineStyle(2);

  TLine *tl6 = new TLine(-1.5,1.05,1.5,1.05);
  TLine *tl7 = new TLine(-1.5,0.95,1.5,0.95);
  TLine *tl8 = new TLine(-1.5,1.,1.5,1.);
  tl6->SetLineStyle(2);
  tl7->SetLineStyle(2);
  tl8->SetLineStyle(2);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_ncs_130_140 = new TCanvas("c_ncs_130_140","",1200,240);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_130_140->Divide(5,1); 
/*
      c_ncs_130_140 ->cd(1);
      h_ncs_2_closure[0][3]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[0][3]->GetXaxis()->SetRangeUser(0.,30.);
      h_ncs_2_closure[0][3]->GetYaxis()->SetLabelSize(0.05);
      h_ncs_2_closure[0][3]->GetXaxis()->SetLabelSize(0.05);
      h_ncs_2_closure[0][3]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[0][3]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
      h_ncs_2_closure[0][3]->GetYaxis()->SetTitleSize(0.08);
      h_ncs_2_closure[0][3]->GetYaxis()->CenterTitle();
      h_ncs_2_closure[0][3]->GetXaxis()->SetTitle("nPF cand");
      h_ncs_2_closure[0][3]->GetXaxis()->SetTitleSize(0.05);
      h_ncs_2_closure[0][3]->GetXaxis()->CenterTitle();
      h_ncs_2_closure[0][3]->GetYaxis()->SetTitleOffset(1.);
      h_ncs_2_closure[0][3]->SetLineColor(kBlue);
      h_ncs_2_closure[0][3]->Draw("same");
      h_ncs_2_closure_corr[0][3]->SetLineColor(kRed);
      h_ncs_2_closure_corr[0][3]->Draw("same");
      tl9->Draw("same");
*/
      c_ncs_130_140 ->cd(1);

      h_ncs_2_dist_q[0][3]->GetYaxis()->SetRangeUser(0.,0.21);
      h_ncs_2_dist_q[0][3]->GetXaxis()->SetRangeUser(0.,20.);
      h_ncs_2_dist_q[0][3]->GetYaxis()->SetLabelSize(0.05);
      h_ncs_2_dist_q[0][3]->GetXaxis()->SetLabelSize(0.05);
      h_ncs_2_dist_q[0][3]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_dist_q[0][3]->GetYaxis()->SetTitle("fraction of inclusive jets");
      h_ncs_2_dist_q[0][3]->GetYaxis()->SetTitleSize(0.08);
      h_ncs_2_dist_q[0][3]->GetYaxis()->SetTitleOffset(1.);
      h_ncs_2_dist_q[0][3]->GetYaxis()->CenterTitle();
      h_ncs_2_dist_q[0][3]->GetXaxis()->SetTitle("nPF cand");
      h_ncs_2_dist_q[0][3]->GetXaxis()->CenterTitle();
      h_ncs_2_dist_q[0][3]->GetXaxis()->SetTitleSize(0.08);
      h_ncs_2_dist_q[0][3]->GetXaxis()->SetTitleOffset(1.);
      h_ncs_2_dist_q[0][3]->GetXaxis()->CenterTitle();
      h_ncs_2_dist_q[0][3]->Draw("same");
      h_ncs_2_dist_g[0][3]->Draw("same");
      cout<<"here"<<endl;
      TLegend *leg0 = new TLegend(0.25,0.75,0.75,0.9);
      leg0 ->SetLineColor(kWhite);
      leg0 ->SetFillColor(0);
      leg0->AddEntry(h_ncs_2_dist_q[0][3],"quark jets","lepf");
      leg0->AddEntry(h_ncs_2_dist_g[0][3],"gluon jets","lepf");
      //legend ->AddEntry((TObject*)0, ptVars[3].c_str(), "");
      leg0 ->Draw("same");

      TLegend *leg1 = new TLegend(0.25,0.75,0.75,0.9);
      leg1 ->SetLineColor(kWhite);
      leg1 ->SetFillColor(0);
      leg1 ->AddEntry((TObject*)0,"Pythia","");
      leg1 ->Draw("same");
/*
      TLegend *leg2 = new TLegend(0.25,0.75,0.75,0.9);
      leg2 ->SetLineColor(kWhite);
      leg2 ->SetFillColor(0);
      leg2 ->AddEntry(h_ncs_2_closure[0][3],"Pre-correction jets", "lepf");
      leg2 ->AddEntry(h_ncs_2_closure_corr[0][3],"Corrected jets", "lepf");
      leg2 ->Draw("same");
*/  

      for(int i=1;i<nCBins;i++){
/*  
      c_ncs_130_140 ->cd(i+1);
      h_ncs_2_closure[5-i][3]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[5-i][3]->GetXaxis()->SetRangeUser(0.,30.);
      h_ncs_2_closure[5-i][3]->GetYaxis()->SetLabelSize(0.05);
      h_ncs_2_closure[5-i][3]->GetXaxis()->SetLabelSize(0.05);
      h_ncs_2_closure[5-i][3]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[5-i][3]->SetTitle(centVars[4-i].c_str());
      h_ncs_2_closure[5-i][3]->GetXaxis()->SetTitle("nCS cand");
      h_ncs_2_closure[5-i][3]->GetXaxis()->SetTitleSize(0.05);
      h_ncs_2_closure[5-i][3]->GetXaxis()->CenterTitle();
      h_ncs_2_closure[5-i][3]->SetLineColor(kBlue);
      h_ncs_2_closure[5-i][3]->Draw("same");
      h_ncs_2_closure_corr[5-i][3]->SetLineColor(kRed);
      h_ncs_2_closure_corr[5-i][3]->Draw("same");
      tl9->Draw("same");
*/
      c_ncs_130_140 ->cd(i+1);
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->SetRangeUser(0.,0.21);
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->SetRangeUser(0.,20.);
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->SetLabelSize(0.05);
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->SetLabelSize(0.05);
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_dist_q[5-i][3]->SetTitle(centVars[4-i].c_str());
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->SetTitleSize(0.08);
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->SetTitleOffset(1.);
      h_ncs_2_dist_q[5-i][3]->GetYaxis()->CenterTitle();
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->SetTitle("nPF cand");
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->SetTitleSize(0.08);
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->CenterTitle();
      h_ncs_2_dist_q[5-i][3]->GetXaxis()->SetTitleOffset(1.);
      h_ncs_2_dist_q[5-i][3]->Draw("same");
      h_ncs_2_dist_g[5-i][3]->Draw("same");

      if(i==1){
      TLegend *leg2 = new TLegend(0.25,0.75,0.75,0.9);
      leg2 ->SetLineColor(kWhite);
      leg2 ->SetFillColor(0);
      leg2 ->AddEntry((TObject*)0, ptVars[3].c_str(), "");
      leg2 ->Draw("same");
      }
      if(i==2){
      TLegend *legend1 = new TLegend(0.25,0.75,0.75,0.9);
      legend1 ->SetLineColor(kWhite);
      legend1 ->SetFillColor(0);
      legend1 ->AddEntry((TObject*)0, "|#eta| < 1.6", "");
      legend1 ->Draw("same");
      }
      TLegend *leg2 = new TLegend(0.25,0.75,0.75,0.9);
      leg2 ->SetLineColor(kWhite);
      leg2 ->SetFillColor(0);
      leg2->AddEntry((TObject*)0,"P+H","");
      leg2 ->Draw("same");
/*      
      if(i==1){
      TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[3].c_str(), "");
      legend ->Draw("same");
      }
      if(i==2){
      TLegend *legend1 = new TLegend(0.25,0.75,0.75,0.9);
      legend1 ->SetLineColor(kWhite);
      legend1 ->AddEntry((TObject*)0, "|#eta| < 1.6", "");
      legend1 ->Draw("same");
      }
      TLegend *leg2 = new TLegend(0.25,0.75,0.75,0.9);
      leg2 ->SetLineColor(kWhite);
      leg2->AddEntry((TObject*)0,"P+H","");
      //legend ->AddEntry((TObject*)0, ptVars[3].c_str(), "");
      leg2 ->Draw("same"); 
*/
    }

/// gen pt resolution plots

  TCanvas *c_res[10];
  for(int c=0;c<10;c++){
    sprintf(saythis,"c_res_%d",c);
    c_res[c] = new TCanvas(saythis,"",900,1200);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    c_res[c]->Divide(4,5);

    for(int j=(5*c); j<(5+(5*c)); j++){  
    for(int k=0; k<10; k++){  
      for(int i=0; i<4; i++){

      c_res[c]->cd(i+(4*k)+1);
      h_ncs_2_res[4-i][j]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[4-i][j]->GetBinContent(h_ncs_2_res_corr[4-i][j]->GetMaximumBin()));
      //h_ncs_2_res[3-i][j]->GetXaxis()->SetRangeUser(0.,2.);
      h_ncs_2_res[4-i][j]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[4-i][j]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[4-i][j]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[4-i][j]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[4-i][j]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[4-i][j]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[4-i][j]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[4-i][j]->SetLineColor(kRed);
      h_ncs_2_res_corr[4-i][j]->Draw("same");
      h_ncs_2_res[4-i][j]->SetLineColor(kBlue);
      h_ncs_2_res[4-i][j]->Draw("same");

      h_ncs_2_res[4-i][j]->Fit(f_ncs_res[4-i][j],"Q M R N","",0.,2.);
      h_ncs_2_res[4-i][j]->Fit(f_ncs_res_2[4-i][j],"Q M R","",(f_ncs_res[4-i][j]->GetParameter(1))*0.7,(f_ncs_res[4-i][j]->GetParameter(1))*1.3);
      par_res[4-i][j] = fabs(f_ncs_res_2[4-i][j]->GetParameter(2));
      par_res_err[4-i][j] = fabs(f_ncs_res_2[4-i][j]->GetParError(2));

      h_ncs_2_res_corr[4-i][j]->Fit(f_ncs_res_corr[4-i][j],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[4-i][j]->Fit(f_ncs_res_corr_2[4-i][j],"Q M R","sames",(f_ncs_res_corr[4-i][j]->GetParameter(1))*0.7,(f_ncs_res_corr[4-i][j]->GetParameter(1))*1.3);
      par_res_corr[4-i][j] = fabs(f_ncs_res_corr_2[4-i][j]->GetParameter(2));
      par_res_corr_err[4-i][j] = fabs(f_ncs_res_corr_2[4-i][j]->GetParError(2));

      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->SetFillColor(0);
      legend ->AddEntry((TObject*)0, ptVars[k+(5*c)].c_str(), "");
      legend ->Draw("same");

      }
    }
    }
  }

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin3=0; ibin3<nptBins; ibin3++){
      h_res[ibin]->SetBinContent(ibin3+1,par_res[ibin][ibin3]);
      h_res[ibin]->SetBinError(ibin3+1,par_res_err[ibin][ibin3]);

      h_res_corr[ibin]->SetBinContent(ibin3+1,par_res_corr[ibin][ibin3]);
      h_res_corr[ibin]->SetBinError(ibin3+1,par_res_corr_err[ibin][ibin3]);
/*
      h_mean[ibin]->SetBinContent(ibin3+1,par_mean[ibin][ibin3]);
      h_mean[ibin]->SetBinError(ibin3+1,par_mean_err[ibin][ibin3]);

      h_mean_corr[ibin]->SetBinContent(ibin3+1,par_mean_corr[ibin][ibin3]);
      h_mean_corr[ibin]->SetBinError(ibin3+1,par_mean_corr_err[ibin][ibin3]);  
*/
    }
  }
/*
  TCanvas *c_res_final = new TCanvas("c_res_final","",1200,240);
  gStyle->SetOptFit(1);
  c_res_final->Divide(5,1);
  for(int i=0; i<nCBins; i++){
    c_res_final->cd(i+2);
    h_res[4-i]->Rebin(4);
    h_res[4-i]->Scale(1./4.);
    h_res_corr[4-i]->Rebin(4);
    h_res_corr[4-i]->Scale(1./4.);
    h_res[4-i]->GetXaxis()->SetRangeUser(100.,350.);
    h_res[4-i]->GetYaxis()->SetRangeUser(0.05,0.25);
    h_res[4-i]->SetMarkerColor(4);
    h_res[4-i]->SetMarkerSize(0.5);
    h_res[4-i]->SetMarkerStyle(21);
    h_res_corr[4-i]->SetMarkerSize(0.5);
    h_res_corr[4-i]->SetMarkerStyle(21);
    h_res[4-i]->SetTitle(centVars[3-i].c_str());
    h_res[4-i]->GetYaxis()->SetTitle("#sigma (recopT/genpT)");
    h_res[4-i]->GetYaxis()->SetLabelSize(0.05);
    h_res[4-i]->GetXaxis()->SetTitle("gen pT");
    h_res[4-i]->SetMarkerColor(kBlack);
    h_res[4-i]->SetLineColor(kBlack);
    h_res[4-i]->Draw("same");
    h_res_corr[4-i]->SetMarkerColor(kRed);
    h_res_corr[4-i]->SetLineColor(kRed);
    h_res_corr[4-i]->Draw("same");
    
    if(i==2){
      TLegend *legend1 = new TLegend(0.25,0.75,0.75,0.9);
      legend1 ->SetLineColor(kWhite);
      legend1 ->AddEntry((TObject*)0, "|#eta| < 1.6", "");
      legend1 ->Draw("same");
    }

      TLegend *leg2 = new TLegend(0.25,0.75,0.75,0.9);
      leg2 ->SetLineColor(kWhite);
      leg2->AddEntry((TObject*)0,"P+H","");
      //legend ->AddEntry((TObject*)0, ptVars[8].c_str(), "");
      leg2 ->Draw("same");

    if(i==3){
      TLegend *leg_res = new TLegend(0.5,0.8,0.99,0.99);
      leg_res->SetLineColor(kWhite);
      leg_res->AddEntry(h_res[0][0], "Pre-correction jets", "lepf");
      leg_res->AddEntry(h_res_corr[0][0], "Corrected jets", "lepf");
      leg_res->Draw("same");
    }
  }
*/
////

//////////// pt closures

  TCanvas *c_closure = new TCanvas("c_closure","",1200,240);
  c_closure->Divide(5,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<5; i++){
  
  if(i==0){
    c_closure->cd(i+1);
    
    h_jt_closure_ref_ncs2_px[i] = h_jt_closure_ref_ncs2[i]->ProfileX();
    h_jt_closure_ref_ncs2_px[i]->Rebin(2);
    h_jt_closure_ref_ncs2_px[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_ncs2_px[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_ncs2_px[i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_ncs2_px[i]->GetXaxis()->CenterTitle();
    h_jt_closure_ref_ncs2_px[i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_ncs2_px[i]->GetYaxis()->SetTitleSize(0.08);
    h_jt_closure_ref_ncs2_px[i]->GetYaxis()->CenterTitle();
    h_jt_closure_ref_ncs2_px[i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_ncs2_px[i]->SetLineColor(kBlack);
    h_jt_closure_ref_ncs2_px[i]->SetMarkerStyle(24);
    h_jt_closure_ref_ncs2_px[i]->SetMarkerSize(0.4);
    h_jt_closure_ref_ncs2_px[i]->Draw("e1 same");
        
    h_jt_closure_q_ncs2_px[i] = h_jt_closure_q_ncs2[i]->ProfileX();
    h_jt_closure_q_ncs2_px[i]->Rebin(2);
    h_jt_closure_q_ncs2_px[i]->SetMarkerColor(kBlue);
    h_jt_closure_q_ncs2_px[i]->SetLineColor(kBlue);
    h_jt_closure_q_ncs2_px[i]->SetMarkerStyle(24);
    h_jt_closure_q_ncs2_px[i]->SetMarkerSize(0.4);
    h_jt_closure_q_ncs2_px[i]->Draw("e1 same");
    h_jt_closure_g_ncs2_px[i] = h_jt_closure_g_ncs2[i]->ProfileX();
    h_jt_closure_g_ncs2_px[i]->Rebin(2);
    h_jt_closure_g_ncs2_px[i]->SetMarkerColor(kRed);
    h_jt_closure_g_ncs2_px[i]->SetLineColor(kRed);
    h_jt_closure_g_ncs2_px[i]->SetMarkerStyle(24);
    h_jt_closure_g_ncs2_px[i]->SetMarkerSize(0.4);
    h_jt_closure_g_ncs2_px[i]->Draw("e1 same");
    h_jt_closure_ref_nocorr_px[i] = h_jt_closure_ref_nocorr[i]->ProfileX();
    
    h_jt_closure_ref_nocorr_px[i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_nocorr_px[i]->Draw("e0 same");
    h_jt_closure_q_nocorr_px[i] = h_jt_closure_q_nocorr[i]->ProfileX();
    h_jt_closure_q_nocorr_px[i]->Rebin(2);
    h_jt_closure_q_nocorr_px[i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_px[i]->Draw("e0 same");
    h_jt_closure_g_nocorr_px[i] = h_jt_closure_g_nocorr[i]->ProfileX();
    h_jt_closure_g_nocorr_px[i]->Rebin(2);
    h_jt_closure_g_nocorr_px[i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_px[i]->Draw("e0 same");

    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");  
    
    TLegend *leg0 = new TLegend(0.4,0.75,0.85,0.89);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry((TObject*)0, "Pythia", "");
      leg0->Draw("same");
      TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->AddEntry((TObject*)0, "| #eta | < 1.6", "");
      leg->Draw("same"); 

  }
  else if(i>0){ 
    c_closure->cd(i+1);

    h_jt_closure_ref_ncs2_px[5-i] = h_jt_closure_ref_ncs2[5-i]->ProfileX();
    h_jt_closure_ref_ncs2_px[5-i]->Rebin(2);
    h_jt_closure_ref_ncs2_px[5-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_ncs2_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_ncs2_px[5-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_ncs2_px[5-i]->GetXaxis()->CenterTitle();
    h_jt_closure_ref_ncs2_px[5-i]->SetTitle(centVars[4-i].c_str());
    h_jt_closure_ref_ncs2_px[5-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_ncs2_px[5-i]->SetLineColor(kBlack);
    h_jt_closure_ref_ncs2_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_ref_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_ref_ncs2_px[5-i]->Draw("e1 same");
    h_jt_closure_q_ncs2_px[5-i] = h_jt_closure_q_ncs2[5-i]->ProfileX();
    h_jt_closure_q_ncs2_px[5-i]->Rebin(2);
    h_jt_closure_q_ncs2_px[5-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_ncs2_px[5-i]->SetLineColor(kBlue);
    h_jt_closure_q_ncs2_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_q_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_q_ncs2_px[5-i]->Draw("e1 same");
    h_jt_closure_g_ncs2_px[5-i] = h_jt_closure_g_ncs2[5-i]->ProfileX();
    h_jt_closure_g_ncs2_px[5-i]->Rebin(2);
    h_jt_closure_g_ncs2_px[5-i]->SetMarkerColor(kRed);
    h_jt_closure_g_ncs2_px[5-i]->SetLineColor(kRed);
    h_jt_closure_g_ncs2_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_g_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_g_ncs2_px[5-i]->Draw("e1 same");
    h_jt_closure_ref_nocorr_px[5-i] = h_jt_closure_ref_nocorr[5-i]->ProfileX();
    h_jt_closure_ref_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[5-i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_nocorr_px[5-i]->Draw("e0 same");
    h_jt_closure_q_nocorr_px[5-i] = h_jt_closure_q_nocorr[5-i]->ProfileX();
    h_jt_closure_q_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_q_nocorr_px[5-i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_px[5-i]->Draw("e0 same");
    h_jt_closure_g_nocorr_px[5-i] = h_jt_closure_g_nocorr[5-i]->ProfileX();
    h_jt_closure_g_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_g_nocorr_px[5-i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_px[5-i]->Draw("e0 same");
    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");

    if(i==1){
      TLegend *leg1 = new TLegend(0.2,0.15,0.85,0.4);
      leg1->SetLineColor(kWhite);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_jt_closure_ref_ncs2_px[4], "Corrected Incl. Jets", "lepf");
      leg1->AddEntry(h_jt_closure_q_ncs2_px[4], "Corrected Quark Jets", "lepf");
      leg1->AddEntry(h_jt_closure_g_ncs2_px[4], "Corrected Gluon Jets", "lepf");
      leg1->Draw("same");
    }

    if(i==2){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(kWhite);
      leg2->SetFillColor(0);
      leg2->AddEntry(h_jt_closure_ref_nocorr_px[3], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_jt_closure_q_nocorr_px[3], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_jt_closure_g_nocorr_px[3], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    }

    TLegend *leg3 = new TLegend(0.4,0.75,0.85,0.89);
      leg3->SetFillColor(0);
      leg3->SetLineColor(kWhite);
      leg3->AddEntry((TObject*)0, "P+H", "");
      leg3->Draw("same");    
    }  
  }

//////////// pt closures

  TCanvas *c_closure_nocorr = new TCanvas("c_closure_nocorr","",1200,240);
  c_closure_nocorr->Divide(5,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<5; i++){
  
  if(i==0){
    c_closure_nocorr->cd(i+1);
    
    h_jt_closure_ref_nocorr_px[i] = h_jt_closure_ref_nocorr[i]->ProfileX();
    h_jt_closure_ref_nocorr_px[i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_nocorr_px[i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_nocorr_px[i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_nocorr_px[i]->GetXaxis()->CenterTitle();
    h_jt_closure_ref_nocorr_px[i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_nocorr_px[i]->GetYaxis()->SetTitleSize(0.08);
    h_jt_closure_ref_nocorr_px[i]->GetYaxis()->CenterTitle();
    h_jt_closure_ref_nocorr_px[i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_nocorr_px[i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[i]->SetMarkerStyle(24);
    h_jt_closure_ref_nocorr_px[i]->SetMarkerSize(0.4);
    h_jt_closure_ref_nocorr_px[i]->Draw("e1 same");
        
    h_jt_closure_q_nocorr_px[i] = h_jt_closure_q_nocorr[i]->ProfileX();
    h_jt_closure_q_nocorr_px[i]->Rebin(2);
    h_jt_closure_q_nocorr_px[i]->SetMarkerColor(kBlue);
    h_jt_closure_q_nocorr_px[i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_px[i]->SetMarkerStyle(24);
    h_jt_closure_q_nocorr_px[i]->SetMarkerSize(0.4);
    h_jt_closure_q_nocorr_px[i]->Draw("e1 same");
    h_jt_closure_g_nocorr_px[i] = h_jt_closure_g_nocorr[i]->ProfileX();
    h_jt_closure_g_nocorr_px[i]->Rebin(2);
    h_jt_closure_g_nocorr_px[i]->SetMarkerColor(kRed);
    h_jt_closure_g_nocorr_px[i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_px[i]->SetMarkerStyle(24);
    h_jt_closure_g_nocorr_px[i]->SetMarkerSize(0.4);
    h_jt_closure_g_nocorr_px[i]->Draw("e1 same");
    h_jt_closure_ref_nocorr_px[i] = h_jt_closure_ref_nocorr[i]->ProfileX();
    
    h_jt_closure_ref_nocorr_px[i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[i]->GetXaxis()->SetRangeUser(100.,500.);
    //h_jt_closure_ref_nocorr_px[i]->Draw("e1 same");
    h_jt_closure_q_nocorr_px[i] = h_jt_closure_q_nocorr[i]->ProfileX();
    h_jt_closure_q_nocorr_px[i]->Rebin(2);
    h_jt_closure_q_nocorr_px[i]->SetLineColor(kBlue);
    //h_jt_closure_q_nocorr_px[i]->Draw("e1 same");
    h_jt_closure_g_nocorr_px[i] = h_jt_closure_g_nocorr[i]->ProfileX();
    h_jt_closure_g_nocorr_px[i]->Rebin(2);
    h_jt_closure_g_nocorr_px[i]->SetLineColor(kRed);
    //h_jt_closure_g_nocorr_px[i]->Draw("e1 same");

    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");  
    
    TLegend *leg0 = new TLegend(0.4,0.75,0.85,0.89);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry((TObject*)0, "Pythia", "");
      leg0->Draw("same");
      TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->AddEntry((TObject*)0, "| #eta | < 1.6", "");
      leg->Draw("same"); 

  }
  else if(i>0){ 
    c_closure_nocorr->cd(i+1);

    h_jt_closure_ref_nocorr_px[5-i] = h_jt_closure_ref_nocorr[5-i]->ProfileX();
    h_jt_closure_ref_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[5-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_nocorr_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_nocorr_px[5-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_nocorr_px[5-i]->GetXaxis()->CenterTitle();
    h_jt_closure_ref_nocorr_px[5-i]->SetTitle(centVars[4-i].c_str());
    h_jt_closure_ref_nocorr_px[5-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_nocorr_px[5-i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_ref_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_ref_nocorr_px[5-i]->Draw("e1 same");
    h_jt_closure_q_nocorr_px[5-i] = h_jt_closure_q_nocorr[5-i]->ProfileX();
    h_jt_closure_q_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_q_nocorr_px[5-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_nocorr_px[5-i]->SetLineColor(kBlue);
    h_jt_closure_q_nocorr_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_q_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_q_nocorr_px[5-i]->Draw("e1 same");
    h_jt_closure_g_nocorr_px[5-i] = h_jt_closure_g_nocorr[5-i]->ProfileX();
    h_jt_closure_g_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_g_nocorr_px[5-i]->SetMarkerColor(kRed);
    h_jt_closure_g_nocorr_px[5-i]->SetLineColor(kRed);
    h_jt_closure_g_nocorr_px[5-i]->SetMarkerStyle(24);
    h_jt_closure_g_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_jt_closure_g_nocorr_px[5-i]->Draw("e1 same");
    h_jt_closure_ref_nocorr_px[5-i] = h_jt_closure_ref_nocorr[5-i]->ProfileX();
    h_jt_closure_ref_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_ref_nocorr_px[5-i]->SetLineColor(kBlack);
    h_jt_closure_ref_nocorr_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    //h_jt_closure_ref_nocorr_px[5-i]->Draw("e1 same");
    h_jt_closure_q_nocorr_px[5-i] = h_jt_closure_q_nocorr[5-i]->ProfileX();
    h_jt_closure_q_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_q_nocorr_px[5-i]->SetLineColor(kBlue);
    //h_jt_closure_q_nocorr_px[5-i]->Draw("e1 same");
    h_jt_closure_g_nocorr_px[5-i] = h_jt_closure_g_nocorr[5-i]->ProfileX();
    h_jt_closure_g_nocorr_px[5-i]->Rebin(2);
    h_jt_closure_g_nocorr_px[5-i]->SetLineColor(kRed);
    //h_jt_closure_g_nocorr_px[5-i]->Draw("e1 same");
    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");

    if(i==1){
      TLegend *leg1 = new TLegend(0.1,0.15,0.85,0.25);
      leg1->SetLineColor(kWhite);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_jt_closure_ref_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg1->AddEntry(h_jt_closure_q_nocorr_px[4], "Pre-correction Quark Jets", "lepf");
      leg1->AddEntry(h_jt_closure_g_nocorr_px[4], "Pre-correction Gluon Jets", "lepf");
      leg1->Draw("same");
    }
/*
    if(i==1){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(kWhite);
      leg2->SetFillColor(0);
      leg2->AddEntry(h_jt_closure_ref_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_jt_closure_q_nocorr_px[4], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_jt_closure_g_nocorr_px[4], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    }
*/
    TLegend *leg3 = new TLegend(0.4,0.75,0.85,0.89);
      leg3->SetFillColor(0);
      leg3->SetLineColor(kWhite);
      leg3->AddEntry((TObject*)0, "P+H", "");
      leg3->Draw("same");    
    }  
  }  

//////////////////nocorr eta closure

  TCanvas *c_closure_eta_nocorr = new TCanvas("c_closure_eta_nocorr","",1200,240);
  c_closure_eta_nocorr->Divide(5,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<5; i++){
  
  if(i==0){
    c_closure_eta_nocorr->cd(i+1);
    
    h_eta_closure_nocorr_px[i] = h_eta_closure_nocorr[i]->ProfileX();
    //h_eta_closure_nocorr_px[i]->Rebin(2);
    h_eta_closure_nocorr_px[i]->GetYaxis()->SetRangeUser(0.85,1.12);
    h_eta_closure_nocorr_px[i]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_eta_closure_nocorr_px[i]->GetXaxis()->SetTitle("#eta");
    h_eta_closure_nocorr_px[i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_eta_closure_nocorr_px[i]->SetMarkerColor(kBlack);
    h_eta_closure_nocorr_px[i]->SetLineColor(kBlack);
    h_eta_closure_nocorr_px[i]->SetMarkerStyle(24);
    h_eta_closure_nocorr_px[i]->SetMarkerSize(0.4);
    h_eta_closure_nocorr_px[i]->Draw("e1 same");
        
    h_eta_closure_q_nocorr_px[i] = h_eta_closure_q_nocorr[i]->ProfileX();
    //h_eta_closure_q_nocorr_px[i]->Rebin(2);
    h_eta_closure_q_nocorr_px[i]->SetMarkerColor(kBlue);
    h_eta_closure_q_nocorr_px[i]->SetLineColor(kBlue);
    h_eta_closure_q_nocorr_px[i]->SetMarkerStyle(24);
    h_eta_closure_q_nocorr_px[i]->SetMarkerSize(0.4);
    h_eta_closure_q_nocorr_px[i]->Draw("e1 same");
    h_eta_closure_g_nocorr_px[i] = h_eta_closure_g_nocorr[i]->ProfileX();
    //h_eta_closure_g_nocorr_px[i]->Rebin(2);
    h_eta_closure_g_nocorr_px[i]->SetMarkerColor(kRed);
    h_eta_closure_g_nocorr_px[i]->SetLineColor(kRed);
    h_eta_closure_g_nocorr_px[i]->SetMarkerStyle(24);
    h_eta_closure_g_nocorr_px[i]->SetMarkerSize(0.4);
    h_eta_closure_g_nocorr_px[i]->Draw("e1 same");
    
    tl6->Draw("same");
    tl7->Draw("same");
    tl8->Draw("same");
    //tl4->Draw("same");  

    TLegend *leg0 = new TLegend(0.4,0.75,0.85,0.89);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry((TObject*)0, "Pythia", "");
      leg0->Draw("same");
      TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->AddEntry((TObject*)0, "Gen p_{T} > 120 GeV", "");
      leg->Draw("same"); 

  }

  else if(i>0){ 
    c_closure_eta_nocorr->cd(i+1);

    h_eta_closure_nocorr_px[5-i] = h_eta_closure_nocorr[5-i]->ProfileX();
    //h_eta_closure_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_nocorr_px[5-i]->GetYaxis()->SetRangeUser(0.85,1.12);
    h_eta_closure_nocorr_px[5-i]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_eta_closure_nocorr_px[5-i]->GetXaxis()->SetTitle("#eta");
    //h_eta_closure_nocorr_px[5-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_eta_closure_nocorr_px[5-i]->SetTitle(centVars[4-i].c_str());
    h_eta_closure_nocorr_px[5-i]->SetMarkerColor(kBlack);
    h_eta_closure_nocorr_px[5-i]->SetLineColor(kBlack);
    h_eta_closure_nocorr_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_nocorr_px[5-i]->Draw("e1 same");
    h_eta_closure_q_nocorr_px[5-i] = h_eta_closure_q_nocorr[5-i]->ProfileX();
    //h_eta_closure_q_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_q_nocorr_px[5-i]->SetMarkerColor(kBlue);
    h_eta_closure_q_nocorr_px[5-i]->SetLineColor(kBlue);
    h_eta_closure_q_nocorr_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_q_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_q_nocorr_px[5-i]->Draw("e1 same");
    h_eta_closure_g_nocorr_px[5-i] = h_eta_closure_g_nocorr[5-i]->ProfileX();
    //h_eta_closure_g_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_g_nocorr_px[5-i]->SetMarkerColor(kRed);
    h_eta_closure_g_nocorr_px[5-i]->SetLineColor(kRed);
    h_eta_closure_g_nocorr_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_g_nocorr_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_g_nocorr_px[5-i]->Draw("e1 same");
    tl6->Draw("same");
    tl7->Draw("same");
    tl8->Draw("same");

    if(i==1){
      TLegend *leg1 = new TLegend(0.2,0.15,0.85,0.4);
      leg1->SetLineColor(kWhite);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_eta_closure_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg1->AddEntry(h_eta_closure_q_nocorr_px[4], "Pre-correction Quark Jets", "lepf");
      leg1->AddEntry(h_eta_closure_g_nocorr_px[4], "Pre-correction Gluon Jets", "lepf");
      leg1->Draw("same");
    }
/*
    if(i==1){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(kWhite);
      leg2->AddEntry(h_eta_closure_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_eta_closure_q_nocorr_px[4], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_eta_closure_g_nocorr_px[4], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    }
*/
    TLegend *leg3 = new TLegend(0.4,0.75,0.85,0.89);
      leg3->SetLineColor(kWhite);
      leg3->SetFillColor(0);
      leg3->AddEntry((TObject*)0, "P+H", "");
      leg3->Draw("same");    
    }  
  }

//////////////////eta closures

  TCanvas *c_closure_eta = new TCanvas("c_closure_eta","",1200,240);
  c_closure_eta->Divide(5,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<5; i++){
  
  if(i==0){
    c_closure_eta->cd(i+1);
    
    h_eta_closure_ncs2_px[i] = h_eta_closure_ncs2[i]->ProfileX();
    //h_eta_closure_ncs2_px[i]->Rebin(2);
    h_eta_closure_ncs2_px[i]->GetYaxis()->SetRangeUser(0.85,1.12);
    h_eta_closure_ncs2_px[i]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_eta_closure_ncs2_px[i]->GetXaxis()->SetTitle("#eta");
    h_eta_closure_ncs2_px[i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_eta_closure_ncs2_px[i]->SetMarkerColor(kBlack);
    h_eta_closure_ncs2_px[i]->SetLineColor(kBlack);
    h_eta_closure_ncs2_px[i]->SetMarkerStyle(24);
    h_eta_closure_ncs2_px[i]->SetMarkerSize(0.4);
    h_eta_closure_ncs2_px[i]->Draw("e1 same");
        
    h_eta_closure_q_ncs2_px[i] = h_eta_closure_q_ncs2[i]->ProfileX();
    //h_eta_closure_q_ncs2_px[i]->Rebin(2);
    h_eta_closure_q_ncs2_px[i]->SetMarkerColor(kBlue);
    h_eta_closure_q_ncs2_px[i]->SetLineColor(kBlue);
    h_eta_closure_q_ncs2_px[i]->SetMarkerStyle(24);
    h_eta_closure_q_ncs2_px[i]->SetMarkerSize(0.4);
    h_eta_closure_q_ncs2_px[i]->Draw("e1 same");
    h_eta_closure_g_ncs2_px[i] = h_eta_closure_g_ncs2[i]->ProfileX();
    //h_eta_closure_g_ncs2_px[i]->Rebin(2);
    h_eta_closure_g_ncs2_px[i]->SetMarkerColor(kRed);
    h_eta_closure_g_ncs2_px[i]->SetLineColor(kRed);
    h_eta_closure_g_ncs2_px[i]->SetMarkerStyle(24);
    h_eta_closure_g_ncs2_px[i]->SetMarkerSize(0.4);
    h_eta_closure_g_ncs2_px[i]->Draw("e1 same");
    
    tl6->Draw("same");
    tl7->Draw("same");
    tl8->Draw("same");
    //tl4->Draw("same");  

    TLegend *leg0 = new TLegend(0.4,0.75,0.85,0.89);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry((TObject*)0, "Pythia", "");
      leg0->Draw("same");
      TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->AddEntry((TObject*)0, "Gen p_{T} > 120 GeV", "");
      leg->Draw("same"); 

  }

  else if(i>0){ 
    c_closure_eta->cd(i+1);

    h_eta_closure_ncs2_px[5-i] = h_eta_closure_ncs2[5-i]->ProfileX();
    //h_eta_closure_ncs2_px[5-i]->Rebin(2);
    h_eta_closure_ncs2_px[5-i]->GetYaxis()->SetRangeUser(0.85,1.12);
    h_eta_closure_ncs2_px[5-i]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_eta_closure_ncs2_px[5-i]->GetXaxis()->SetTitle("#eta");
    //h_eta_closure_ncs2_px[5-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_eta_closure_ncs2_px[5-i]->SetTitle(centVars[4-i].c_str());
    h_eta_closure_ncs2_px[5-i]->SetMarkerColor(kBlack);
    h_eta_closure_ncs2_px[5-i]->SetLineColor(kBlack);
    h_eta_closure_ncs2_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_ncs2_px[5-i]->Draw("e1 same");
    h_eta_closure_q_ncs2_px[5-i] = h_eta_closure_q_ncs2[5-i]->ProfileX();
    //h_eta_closure_q_ncs2_px[5-i]->Rebin(2);
    h_eta_closure_q_ncs2_px[5-i]->SetMarkerColor(kBlue);
    h_eta_closure_q_ncs2_px[5-i]->SetLineColor(kBlue);
    h_eta_closure_q_ncs2_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_q_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_q_ncs2_px[5-i]->Draw("e1 same");
    h_eta_closure_g_ncs2_px[5-i] = h_eta_closure_g_ncs2[5-i]->ProfileX();
    //h_eta_closure_g_ncs2_px[5-i]->Rebin(2);
    h_eta_closure_g_ncs2_px[5-i]->SetMarkerColor(kRed);
    h_eta_closure_g_ncs2_px[5-i]->SetLineColor(kRed);
    h_eta_closure_g_ncs2_px[5-i]->SetMarkerStyle(24);
    h_eta_closure_g_ncs2_px[5-i]->SetMarkerSize(0.4);
    h_eta_closure_g_ncs2_px[5-i]->Draw("e1 same");
    h_eta_closure_nocorr_px[5-i] = h_eta_closure_nocorr[5-i]->ProfileX();
    //h_eta_closure_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_nocorr_px[5-i]->SetLineColor(kBlack);
    h_eta_closure_nocorr_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    //h_eta_closure_nocorr_px[5-i]->Draw("e1 same");
    h_eta_closure_q_nocorr_px[5-i] = h_eta_closure_q_nocorr[5-i]->ProfileX();
    //h_eta_closure_q_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_q_nocorr_px[5-i]->SetLineColor(kBlue);
    //h_eta_closure_q_nocorr_px[5-i]->Draw("e1 same");
    h_eta_closure_g_nocorr_px[5-i] = h_eta_closure_g_nocorr[5-i]->ProfileX();
    //h_eta_closure_g_nocorr_px[5-i]->Rebin(2);
    h_eta_closure_g_nocorr_px[5-i]->SetLineColor(kRed);
    //h_eta_closure_g_nocorr_px[5-i]->Draw("e1 same");
    tl6->Draw("same");
    tl7->Draw("same");
    tl8->Draw("same");

    if(i==1){
      TLegend *leg1 = new TLegend(0.2,0.15,0.85,0.4);
      leg1->SetLineColor(kWhite);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_eta_closure_ncs2_px[4], "Corrected Incl. Jets", "lepf");
      leg1->AddEntry(h_eta_closure_q_ncs2_px[4], "Corrected Quark Jets", "lepf");
      leg1->AddEntry(h_eta_closure_g_ncs2_px[4], "Corrected Gluon Jets", "lepf");
      leg1->Draw("same");
    }
/*
    if(i==1){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(kWhite);
      leg2->AddEntry(h_eta_closure_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_eta_closure_q_nocorr_px[4], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_eta_closure_g_nocorr_px[4], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    }
*/
    TLegend *leg3 = new TLegend(0.4,0.75,0.85,0.89);
      leg3->SetLineColor(kWhite);
      leg3->SetFillColor(0);
      leg3->AddEntry((TObject*)0, "P+H", "");
      leg3->Draw("same");    
    }  
  }  

///reco closure
  TCanvas *c_closure_reco = new TCanvas("c_closure_reco","",1200,240);
  c_closure_reco->Divide(5,1);
  gStyle->SetOptStat(0);
    c_closure_reco->cd(1);

    h_jt_closure_reco_nocorr_px[0]=h_jt_closure_reco_nocorr[0]->ProfileX();
    h_jt_closure_reco_nocorr_px[0]->GetYaxis()->SetRangeUser(1.,1.3);
    h_jt_closure_reco_nocorr_px[0]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_reco_nocorr_px[0]->GetXaxis()->SetTitle("reco pT");
    h_jt_closure_reco_nocorr_px[0]->GetXaxis()->SetTitleSize(0.05);
    h_jt_closure_reco_nocorr_px[0]->GetXaxis()->CenterTitle();
    h_jt_closure_reco_nocorr_px[0]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_reco_nocorr_px[0]->SetMarkerColor(kBlack);
    h_jt_closure_reco_nocorr_px[0]->SetMarkerStyle(29);
    h_jt_closure_reco_nocorr_px[0]->Draw("e1 same");
    tl5->Draw("same");

      TLegend *leg0 = new TLegend(0.5,0.8,0.99,0.99);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry((TObject*)0, "Pythia", "");
      leg0->Draw("same");
      TLegend *leg = new TLegend(0.5,0.8,0.99,0.99);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(0);
      leg->AddEntry((TObject*)0, "| #eta | < 1.6", "");
      leg->Draw("same"); 

  for(int i=1; i<5; i++){
    c_closure_reco->cd(i+1);
    h_jt_closure_reco_nocorr_px[5-i]=h_jt_closure_reco_nocorr[5-i]->ProfileX();
    h_jt_closure_reco_nocorr_px[5-i]->GetYaxis()->SetRangeUser(1.,1.3);
    h_jt_closure_reco_nocorr_px[5-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_reco_nocorr_px[5-i]->GetXaxis()->SetTitle("reco pT");
    h_jt_closure_reco_nocorr_px[5-i]->GetXaxis()->SetTitleSize(0.05);
    h_jt_closure_reco_nocorr_px[5-i]->GetXaxis()->CenterTitle();
    h_jt_closure_reco_nocorr_px[5-i]->SetMarkerColor(kBlack);
    h_jt_closure_reco_nocorr_px[5-i]->SetMarkerStyle(29);
    h_jt_closure_reco_nocorr_px[5-i]->SetTitle(centVars[4-i].c_str());
    h_jt_closure_reco_nocorr_px[5-i]->Draw("e1 same");
    tl5->Draw("same");

    TLegend *leg3 = new TLegend(0.5,0.8,0.99,0.99);
      leg3->SetLineColor(kWhite);
      leg3->SetFillColor(0);
      leg3->AddEntry((TObject*)0, "P+H", "");
      leg3->Draw("same");
        if(i==1){
      TLegend *leg1 = new TLegend(0.5,0.8,0.99,0.99);
      leg1->SetLineColor(kWhite);
      leg1->SetFillColor(0);
      leg1->AddEntry(h_jt_closure_reco_nocorr_px[4], "Pre-correction Incl. Jets", "lepf");
      leg1->Draw("same");
    }

  }




}