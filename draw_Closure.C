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
#define nptBins 40
#define npfbins 21

char saythis[500];

using namespace std;

TString cent[4] = {"0","1","2","3"};
TString pt[41] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40"};

int jt_nbins = 41;
Double_t jt_bin_bounds[41] = {100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.};

void draw_Closure(){

  //TFile *closure_histos = TFile::Open("/home/dhanush/Documents/JFF_corrections/closure_histos_Jul6_header_id145.root");
  TFile *closure_histos = TFile::Open("/home/dhanush/Documents/JFF_corrections/closure_histos_Jul24_cymbal_header_id1_rebin.root");

// defining histos

  TH1D *h_reco[nCBins][nptBins];
  TH1D *h_reco_corr_iter[nCBins];
  TH1D *h_reco_corr_ncs[nCBins];
  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_corr[nCBins];
  TH1D *h_reco_iter[nCBins];
  TH1D *h_gen_full[nCBins];
  TH1D *h_reco_ratio[nCBins];
  TH1D *h_corr_ratio[nCBins];
  TH1D *h_iter_ratio[nCBins];
  TH1D *h_uncorr_ratio[nCBins];
  TH1D *h_reco_recogen2[nCBins];
  TH1D *h_reco_refpt20[nCBins];
  TH1D *h_reco_999[nCBins];
  TH2F *h_ncs_1[nCBins][nptBins];
  TH2F *h_ncs_2[nCBins][nptBins];
  TH2F *h_ncs_2_q[nCBins][nptBins];
  TH2F *h_ncs_2_g[nCBins][nptBins];
  TH2F *h_ncs_2_corr[nCBins][nptBins];
  TH2F *h_ncs_2_corr_q[nCBins][nptBins];
  TH2F *h_ncs_2_corr_g[nCBins][nptBins];
  TH2F *h_jt_closure[nCBins][nptBins];
  TH2F *h_jt_closure_ref_nocorr[nCBins];
  TH2F *h_jt_closure_ref_ncs1[nCBins];
  TH2F *h_jt_closure_ref_ncs2[nCBins];
  TH2F *h_jt_closure_reco_nocorr[nCBins];
  TH2F *h_jt_closure_reco_ncs1[nCBins];
  TH2F *h_jt_closure_reco_ncs2[nCBins];
  TH2F *h_jt_closure_q_nocorr[nCBins];
  TH2F *h_jt_closure_q_ncs1[nCBins];
  TH2F *h_jt_closure_q_ncs2[nCBins];
  TH2F *h_jt_closure_g_nocorr[nCBins];
  TH2F *h_jt_closure_g_ncs1[nCBins];
  TH2F *h_jt_closure_g_ncs2[nCBins];

  TProfile *h_ncs_1_closure[nCBins][nptBins];
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
  TH1D *h_mean_corr[nCBins];    
  TH1D *h_jt_closure_ref_nocorr_ppx[nCBins];

  TProfile *h_jt_closure_ref_nocorr_px[nCBins];
  TProfile *h_jt_closure_ref_ncs1_px[nCBins];
  TProfile *h_jt_closure_ref_ncs2_px[nCBins];
  TProfile *h_jt_closure_reco_nocorr_px[nCBins];
  TProfile *h_jt_closure_reco_ncs1_px[nCBins];
  TProfile *h_jt_closure_reco_ncs2_px[nCBins];
  TProfile *h_jt_closure_q_nocorr_px[nCBins];
  TProfile *h_jt_closure_q_ncs1_px[nCBins];
  TProfile *h_jt_closure_q_ncs2_px[nCBins];
  TProfile *h_jt_closure_g_nocorr_px[nCBins];
  TProfile *h_jt_closure_g_ncs1_px[nCBins];
  TProfile *h_jt_closure_g_ncs2_px[nCBins];
  TH1D *h_jt_closure_px[nCBins][nptBins];

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

    sprintf(saythis,"h_reco_ratio_cent%d",ibin);
    h_reco_ratio[ibin] = new TH1D("","",500,0.,500.);
    h_reco_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_corr_ratio_cent%d",ibin);
    h_corr_ratio[ibin] = new TH1D("","",45,50.,500.);
    h_corr_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_iter_ratio_cent%d",ibin);
    h_iter_ratio[ibin] = new TH1D("","",45,50.,500.);
    h_iter_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_uncorr_ratio_cent%d",ibin);
    h_uncorr_ratio[ibin] = new TH1D("","",45,50.,500.);
    h_uncorr_ratio[ibin]->Sumw2();

    sprintf(saythis,"h_res_cent%d",ibin);
    h_res[ibin] = new TH1D("","",45,50.,500.);
    h_res[ibin]->Sumw2();

    sprintf(saythis,"h_res_corr_cent%d",ibin);
    h_res_corr[ibin] = new TH1D("","",45,50.,500.);
    h_res_corr[ibin]->Sumw2();

    sprintf(saythis,"h_mean_cent%d",ibin);
    h_mean[ibin] = new TH1D("","",45,50.,500.);
    h_mean[ibin]->Sumw2();

    sprintf(saythis,"h_mean_corr_cent%d",ibin);
    h_mean_corr[ibin] = new TH1D("","",45,50.,500.);
    h_mean_corr[ibin]->Sumw2();

  }

  for(int ibin=0;ibin<nCBins;ibin++){

    h_jt_closure_ref_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_ref_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_nocorr_"+cent[ibin]));
    h_jt_closure_q_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_q_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_nocorr_"+cent[ibin]));
    h_jt_closure_g_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_g_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_nocorr_"+cent[ibin]));
    h_jt_closure_ref_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_ref_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_ref_ncs2_"+cent[ibin]));
    h_jt_closure_q_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_q_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_q_ncs2_"+cent[ibin]));
    h_jt_closure_g_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_g_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_g_ncs2_"+cent[ibin]));
    h_jt_closure_reco_nocorr[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_reco_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_reco_nocorr_"+cent[ibin]));
    h_jt_closure_reco_ncs2[ibin] = (TH2F*)closure_histos->Get((TString)("h_jt_closure_reco_ncs2_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_reco_ncs2_"+cent[ibin]));
    h_jt_closure_reco_nocorr_px[ibin] = h_jt_closure_reco_nocorr[ibin]->ProfileX();
    
    h_reco_full[ibin] = (TH1D*)closure_histos->Get((TString)("h_reco_full_cent"+cent[ibin]))->Clone((TString)("h_reco_full_"+cent[ibin]));
    h_reco_corr[ibin] = (TH1D*)closure_histos->Get((TString)("h_reco_corr_cent"+cent[ibin]))->Clone((TString)("h_reco_corr_"+cent[ibin]));
    h_gen_full[ibin] = (TH1D*)closure_histos->Get((TString)("h_gen_full_cent"+cent[ibin]))->Clone((TString)("h_gen_full_"+cent[ibin]));

    sprintf(saythis,"f_flatcorr_cent%d",ibin);
      f_closure[ibin] = new TF1(saythis, "pol1", 80., 500.);
      f_closure[ibin] ->SetParameter(0,0.995);
      f_closure[ibin] ->SetParameter(1,2e-05);

      sprintf(saythis,"f_flatcorr2_cent%d",ibin);
      f_closure_2[ibin] = new TF1(saythis, "pol0", 80., 500.);
      f_closure_2[ibin] ->SetParameter(0,1.);

    for(int ibin3=0;ibin3<nptBins;ibin3++){
      
      h_ncs_2[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_q[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_q_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_q_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_g[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_nocorr_g_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_g_"+cent[ibin]+"_"+pt[ibin3]));

      h_ncs_2_res[ibin][ibin3] = (TH1D*)closure_histos->Get((TString)("h_closure_nocorr_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_nocorr_res_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_res_corr[ibin][ibin3] = (TH1D*)closure_histos->Get((TString)("h_closure_ncs2_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_closure_ncs2_res_"+cent[ibin]+"_"+pt[ibin3]));
     
      h_ncs_2_corr[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_corr_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_corr_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_corr_q[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_corr_q_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_corr_q_"+cent[ibin]+"_"+pt[ibin3]));
      h_ncs_2_corr_g[ibin][ibin3] = (TH2F*)closure_histos->Get((TString)("h_ncs2_corr_g_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_ncs2_corr_g_"+cent[ibin]+"_"+pt[ibin3]));

      h_reco[ibin][ibin3] = (TH1D*)closure_histos->Get((TString)("h_reco_cent"+cent[ibin]+"_pt"+pt[ibin3]))->Clone((TString)("h_reco_"+cent[ibin]+"_"+pt[ibin3]));
    
      sprintf(saythis,"f_ncs_1_cent%d_pt%d",ibin,ibin3);
      f_ncs_1[ibin][ibin3] = new TF1(saythis, "1.+[1]*(x-[0])", 2., 30.);
      f_ncs_1[ibin][ibin3] ->SetParameter(0,10.);
      f_ncs_1[ibin][ibin3] ->SetParameter(1,-0.015);

      sprintf(saythis,"f_ncs_2_cent%d_pt%d",ibin,ibin3);
      f_ncs_2[ibin][ibin3] = new TF1(saythis, "1.+[1]*(x-[0])", 2., 30.);
      f_ncs_2[ibin][ibin3] ->SetParameter(0,10.);
      f_ncs_2[ibin][ibin3] ->SetParameter(1,-0.015);
      //f_ncs_2[ibin][ibin3] = new TF1(saythis, "pol1", 2., 30.);

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

  /// plotting parameters for ncs corrections
/*
  TFile f_corr;
  f_corr = new TFile("/home/dhanush/Documents/JEC/local/new_corr_files/ncscorrfactors_PbPb5TeV_Apr17.root", "RECREATE");   
  f_corr->cd();
*/
  Double_t x_mean[nptBins];
  Double_t x_mean_err[nptBins];

  Double_t par0_1[nptBins];
  Double_t par1_1[nptBins];
  Double_t par0_err_1[nptBins];
  Double_t par1_err_1[nptBins];

  Double_t par_res[nCBins][nptBins];
  Double_t par_res_err[nCBins][nptBins];

  Double_t par_res_corr[nCBins][nptBins];
  Double_t par_res_corr_err[nCBins][nptBins];

  Double_t par_mean[nCBins][nptBins];
  Double_t par_mean_err[nCBins][nptBins];

  Double_t par_mean_corr[nCBins][nptBins];
  Double_t par_mean_corr_err[nCBins][nptBins];

  Double_t par0_2[nptBins];
  Double_t par1_2[nptBins];
  Double_t par0_err_2[nptBins];
  Double_t par1_err_2[nptBins]; 

  Double_t ax_2[nptBins];
  Double_t ax_2_err[nptBins];
  Double_t ax_2_num[nptBins];
  Double_t ax_2_err_num[nptBins];
  Double_t ax_2_err_den[nptBins];

  // a0, a1 and ax parameter fits

  TF1 *f1_par0[nCBins];
  TF1 *f1_par1[nCBins];
  TF1 *f2_par0[nCBins];
  TF1 *f2_par1[nCBins];
  TF1 *f2_ax[nCBins];

  for(int ibin=0;ibin<nCBins;ibin++){

    for(int ibin3=0;ibin3<nptBins;ibin3++){ 

    //h_ncs_2_res[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProjectionY();
    h_ncs_2_res[ibin][ibin3] -> Rebin(5);
    h_ncs_2_res[ibin][ibin3] -> Scale(1./5.);
    Double_t integral_nocorr = h_ncs_2_res[ibin][ibin3]->Integral();
    h_ncs_2_res[ibin][ibin3]->Scale(1./integral_nocorr);

    h_ncs_2_res_q[ibin][ibin3] = h_ncs_2_q[ibin][ibin3]->ProjectionY();
    h_ncs_2_res_q[ibin][ibin3]->Scale(1./integral_nocorr);

    h_ncs_2_res_g[ibin][ibin3] = h_ncs_2_g[ibin][ibin3]->ProjectionY();
    h_ncs_2_res_g[ibin][ibin3]->Scale(1./integral_nocorr);

    //h_ncs_2_res_corr[ibin][ibin3] = h_ncs_2_corr[ibin][ibin3]->ProjectionY();
    h_ncs_2_res_corr[ibin][ibin3] -> Rebin(5);
    h_ncs_2_res_corr[ibin][ibin3] -> Scale(1./5.);
    Double_t integral_corr = h_ncs_2_res_corr[ibin][ibin3]->Integral();
    h_ncs_2_res_corr[ibin][ibin3]->Scale(1./integral_corr);

    h_ncs_2_res_corr_q[ibin][ibin3] = h_ncs_2_corr_q[ibin][ibin3]->ProjectionY();
    h_ncs_2_res_corr_q[ibin][ibin3]->Scale(1./integral_corr);

    h_ncs_2_res_corr_g[ibin][ibin3] = h_ncs_2_corr_g[ibin][ibin3]->ProjectionY();
    h_ncs_2_res_corr_g[ibin][ibin3]->Scale(1./integral_corr);

    //h_ncs_1_closure[ibin][ibin3] = h_ncs_1[ibin][ibin3]->ProfileX();
    //h_ncs_1_closure[ibin][ibin3]->Fit(((TString)("f_ncs_1_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q R M");
/*
    if(ibin3<5) h_ncs_1_closure[ibin][ibin3]->Fit(((TString)("f_ncs_1_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.); 

    else if(ibin3<10) h_ncs_1_closure[ibin][ibin3]->Fit(((TString)("f_ncs_1_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,30.);
    else if(ibin3<35) h_ncs_1_closure[ibin][ibin3]->Fit(((TString)("f_ncs_1_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",5.,35.);
   
    else h_ncs_1_closure[ibin][ibin3]->Fit(((TString)("f_ncs_1_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",8.,35.);
*/
    par0_1[ibin3] = f_ncs_1[ibin][ibin3]->GetParameter(0);
    par1_1[ibin3] = f_ncs_1[ibin][ibin3]->GetParameter(1);

    par0_err_1[ibin3] = f_ncs_1[ibin][ibin3]->GetParError(0);
    par1_err_1[ibin3] = f_ncs_1[ibin][ibin3]->GetParError(1);

    h_ncs_2_closure[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProfileX();
    h_ncs_2_closure_corr[ibin][ibin3] = h_ncs_2_corr[ibin][ibin3]->ProfileX();
 
    h_ncs_2_dist[ibin][ibin3] = h_ncs_2[ibin][ibin3]->ProjectionX();
    h_ncs_2_dist_q[ibin][ibin3] = h_ncs_2_q[ibin][ibin3]->ProjectionX();
    h_ncs_2_dist_q[ibin][ibin3] -> Scale (1./(h_ncs_2_dist_q[ibin][ibin3])->Integral());
    h_ncs_2_dist_q[ibin][ibin3] -> SetLineColor (kBlue);
    h_ncs_2_dist_g[ibin][ibin3] = h_ncs_2_g[ibin][ibin3]->ProjectionX();
    h_ncs_2_dist_g[ibin][ibin3] -> Scale (1./(h_ncs_2_dist_g[ibin][ibin3])->Integral());
    h_ncs_2_dist_g[ibin][ibin3] -> SetLineColor (kRed);
    //h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q R M");

    if(ibin3<5) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",0.,14.);


    else if(ibin3<10) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.);
    else if(ibin3<13) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,30.);
    else if(ibin3<20) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,35.);
    else if((ibin3==26) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,30.);
    else if(ibin3<30) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,38.);
    else if((ibin3==34) && (ibin==1)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,38.);    
    else if(ibin3<35 && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,35.);
    else if(ibin3<35) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,38.);
    else if((ibin3==37) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",5.,30.);
    else if(ibin3<40) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,35.);
    else if((ibin3==40) && ((ibin==0)||(ibin==1))) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",5.,30.);
    else if(ibin3<45) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",5.,35.);
    else h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",8.,25.);

/*
    else if(ibin3==7) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,18.);
    else if((ibin3==8) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,20.);
    else if(ibin3<9) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,17.);
    else if((ibin3==11) && (ibin==1)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.);
    else if((ibin3==11) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.);
    else if(ibin3<12) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);
    else if((ibin3==13) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",0.,21.);
    else if((ibin3==14) && (ibin==3)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,21.);
    else if((ibin3==27) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,18.);
    else if((ibin3==28) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);
    else if((ibin3==29) && (ibin==1)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,21.);
    else if((ibin3==31) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",0.,19.);
    else if((ibin3==32) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",1.,15.);
    else if((ibin3==30||ibin3==31||ibin3==32) && (ibin==1||ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);
    else if(ibin3==33) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,24.);
    else if((ibin3==34) && (ibin==3||ibin==1)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,23.);
    else if((ibin3==34) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,16.);
    else if(ibin3<35) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,28.);
    else if(ibin3==43 && ibin==1) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,22.);
    else if(ibin3==47 && ibin==2) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,24.);
    else if(ibin3==48 && ibin==0) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,12.);
    else if(ibin3==49 && ibin==1) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,19.);
    else if(ibin3>44 && ibin==0) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,15.);
*/
/*    
    else if(ibin3<9) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",1.,15.);
    else if(ibin3<15) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",1.,19.);
    else if(ibin3<25) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",1.,22.);
    else if(ibin3<30) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,23.);
    else if(ibin3<39) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,21.);
    else if(ibin3<44) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,19.);
    else if(ibin3<50) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,17.);
    else h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,16.); 
*/
    par0_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParameter(0);
    par1_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParameter(1);

    par0_err_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParError(0);
    par1_err_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParError(1);

    // mean pt in the pt bin

    x_mean[ibin3] = h_reco[ibin][ibin3]->GetMean();
    x_mean_err[ibin3] = h_reco[ibin][ibin3]->GetMeanError(); 

    }
  
    f1_par1[ibin] = new TF1((TString)("f1_par1_cent"+cent[ibin]),"([2]/x)*exp(([0]/((x-[3])^2))+([1]/x))",80.,550.);
    f1_par0[ibin] = new TF1((TString)("f1_par0_cent"+cent[ibin]), "exp(([0]/(x-[2]))+[1])", 80., 550.);

    f2_par1[ibin] = new TF1((TString)("f2_par1_cent"+cent[ibin]),"([2]/x)*exp(([0]/((x-[3])^2))+([1]/x))",60.,550.);
    f2_par0[ibin] = new TF1((TString)("f2_par0_cent"+cent[ibin]), "exp(([0]/(x-[2]))+[1])", 60., 550.);
        
    if(ibin==0){ 

      f2_par1[ibin]->SetParameter(0,5.23836e7);
      f2_par1[ibin]->SetParameter(1,-91.8229);
      f2_par1[ibin]->SetParameter(2,-0.675532);
      f2_par1[ibin]->SetParameter(3,5877.84);

      f2_par0[ibin]->SetParameter(0,-137.071);
      f2_par0[ibin]->SetParameter(1,2.97793);
      f2_par0[ibin]->SetParameter(2,-96.3265);
      
      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",60.,470.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",60.,470.);
    }

    else if(ibin==1){

      f2_par1[ibin]->SetParameter(0,170765);
      f2_par1[ibin]->SetParameter(1,-722.664);
      f2_par1[ibin]->SetParameter(2,-9.0319);
      f2_par1[ibin]->SetParameter(3,-73.4127);

      f2_par0[ibin]->SetParameter(0,-118.149);
      f2_par0[ibin]->SetParameter(1,2.94885);
      f2_par0[ibin]->SetParameter(2,-67.7417);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",60.,470.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",60.,470.); 
    }

    else if(ibin==2){

      f2_par1[ibin]->SetParameter(0,128199);
      f2_par1[ibin]->SetParameter(1,-580.536);
      f2_par1[ibin]->SetParameter(2,-7.48496);
      f2_par1[ibin]->SetParameter(3,-66.0678);

      f2_par0[ibin]->SetParameter(0,-115.096);
      f2_par0[ibin]->SetParameter(1,2.93954);
      f2_par0[ibin]->SetParameter(2,-60.6488);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",60.,480.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",60.,480.); 
    }

    else if(ibin==3){

      f2_par1[ibin]->SetParameter(0,81150.7);
      f2_par1[ibin]->SetParameter(1,-417.789);
      f2_par1[ibin]->SetParameter(2,-6.06929);
      f2_par1[ibin]->SetParameter(3,-57.8275);

      f2_par0[ibin]->SetParameter(0,-109.053);
      f2_par0[ibin]->SetParameter(1,2.93768);
      f2_par0[ibin]->SetParameter(2,-51.9022);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",60.,550.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",60.,550.);

    }
 
  }
  
/*
  f1_par0_cent0->Write();
  f1_par0_cent1->Write();
  f1_par0_cent2->Write();
  f1_par0_cent3->Write();
  f1_par1_cent0->Write();
  f1_par1_cent1->Write();
  f1_par1_cent2->Write();
  f1_par1_cent3->Write();

  f2_par0_cent0->Write();
  f2_par0_cent1->Write();
  f2_par0_cent2->Write();
  f2_par0_cent3->Write();
  f2_par1_cent0->Write();
  f2_par1_cent1->Write();
  f2_par1_cent2->Write();
  f2_par1_cent3->Write();
*/

  TLatex *l1[nCBins][nptBins];
  const string centVars[4] = {"Cent 0-10%", "Cent 10-30%", "Cent 30-50%",  "Cent 50-100%"};
  const string ptVars[55] = {"50<pT<60", "60<pT<70", "70<pT<80","80<pT<90", "90<pT<100", "100<pT<110", "110<pT<120","120<pT<130", "130<pT<140", "140<pT<150", "150<pT<160","160<pT<170", "170<pT<180", "180<pT<190", "190<pT<200", "200<pT<210","210<pT<220", "220<pT<230", "230<pT<240", "240<pT<250","250<pT<260", "260<pT<270", "270<pT<280", "280<pT<290","290<pT<300", "300<pT<310", "310<pT<320", "320<pT<330", "330<pT<340", "340<pT<350", "350<pT<360", "360<pT<370", "370<pT<380", "380<pT<390", "390<pT<400", "400<pT<410", "410<pT<420", "420<pT<430", "430<pT<440", "440<pT<450","450<pT<460", "460<pT<470", "470<pT<480", "480<pT<490", "490<pT<500", "500<pT<510", "510<pT<520", "520<pT<530", "530<pT<540","540<pT<550","550<pT<560", "560<pT<570", "570<pT<580", "580<pT<590","590<pT<600"};
  /*
  TCanvas *c_ncs_50_100 = new TCanvas("c_ncs_50_100","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_50_100->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_50_100->cd(i+(4*(j))+1);

      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+5]->GetXaxis()->SetRangeUser(0.,40.);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+5]->GetXaxis()->SetLabelSize(0.09);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetNdivisions(5);
      if(j%5==0)h_ncs_2_closure[3-i][j+5]->SetTitle(centVars[3-i].c_str());
      if(i==0) h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitle("<p_{T}^{reco}/p_{T}^{gen}>");
      h_ncs_2_closure[3-i][j+5]->GetXaxis()->SetTitle("nCS cand");
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+5]->GetXaxis()->SetTitleSize(0.07);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitleOffset(0.51);
      h_ncs_2_closure[3-i][j+5]->SetLineColor(kBlue);
      h_ncs_2_closure[3-i][j+5]->Draw("same");
      h_ncs_2_closure_corr[3-i][j+5]->SetLineColor(kRed);
      h_ncs_2_closure_corr[3-i][j+5]->Draw("same");
      //l1[3-i][j+5] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j+5]->SetTextSize(0.08);
      //l1[3-i][j+5]->Draw("same");
      TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+5].c_str(), "");
      legend ->Draw("same");
    

    }
  }
/*
  TCanvas *c_ncs_130_140 = new TCanvas("c_ncs_130_140","",1200,300);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  c_ncs_130_140->Divide(4,1);
  
    for(int i=0; i<4; i++){
      c_ncs_130_140->cd(i+1);

      h_ncs_2_closure[3-i][8]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][8]->GetXaxis()->SetRangeUser(0.,40.);
      h_ncs_2_closure[3-i][8]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][8]->GetXaxis()->SetLabelSize(0.09);
      h_ncs_2_closure[3-i][8]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][8]->SetTitle(centVars[3-i].c_str());
      if(i==0) h_ncs_2_closure[3-i][8]->GetYaxis()->SetTitle("<p_{T}^{reco}/p_{T}^{gen}>");
      h_ncs_2_closure[3-i][8]->GetXaxis()->SetTitle("nCS cand");
      h_ncs_2_closure[3-i][8]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][8]->GetXaxis()->SetTitleSize(0.07);
      h_ncs_2_closure[3-i][8]->GetYaxis()->SetTitleOffset(0.51);
      h_ncs_2_closure[3-i][8]->SetLineColor(kBlue);
      h_ncs_2_closure[3-i][8]->Draw("same");
      h_ncs_2_closure_corr[3-i][8]->SetLineColor(kRed);
      h_ncs_2_closure_corr[3-i][8]->Draw("same");
      //l1[3-i][8] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][8]->SetTextSize(0.08);
      //l1[3-i][8]->Draw("same");
      if(i==1){
        TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry((TObject*)0, ptVars[8].c_str(), "");
        legend ->Draw("same");
      }
      if(i==0){
        TLegend *legend1 = new TLegend(0.25,0.75,0.75,0.9);
        legend1 ->SetLineColor(kWhite);
        legend1 ->AddEntry(h_ncs_2_closure[3][8], "Pre-correction", "lepf");
        legend1 ->AddEntry(h_ncs_2_closure_corr[3][8], "Corrected", "lepf");
        legend1 ->Draw("same");
      }
      if(i==2){
        TLegend *legend2 = new TLegend(0.25,0.75,0.75,0.9);
        legend2 ->SetLineColor(kWhite);
        legend2 ->AddEntry((TObject*)0, "|#eta| < 1.6", "");
        legend2 ->Draw("same");
      }
    }
*/  

  TCanvas *c_ncs_50_100 = new TCanvas("c_ncs_50_100","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_50_100->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_50_100->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetTitle(ptVars[j].c_str());
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j]->Draw();
      //h_ncs_2_closure[3-i][j+5]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+5]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_100_150 = new TCanvas("c_ncs_100_150","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_100_150->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_100_150->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+5]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+5]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitle(ptVars[j+5].c_str());
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+5]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+5]->Draw();
      //h_ncs_2_closure[3-i][j+5]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+5]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_150_200 = new TCanvas("c_ncs_150_200","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_150_200->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_150_200->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+10]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+10]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetTitle(ptVars[j+10].c_str());
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+10]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+10]->Draw();
     // h_ncs_2_closure[3-i][j+10]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+10]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_200_250 = new TCanvas("c_ncs_200_250","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_200_250->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_200_250->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+15]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+15]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetTitle(ptVars[j+15].c_str());
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+15]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+15]->Draw();
      //h_ncs_2_closure[3-i][j+15]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+15]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_250_300 = new TCanvas("c_ncs_250_300","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_250_300->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_250_300->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+20]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+20]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetTitle(ptVars[j+20].c_str());
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+20]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+20]->Draw();
      //h_ncs_2_closure[3-i][j+20]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+20]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_300_350 = new TCanvas("c_ncs_300_350","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_300_350->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_300_350->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+25]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+25]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetTitle(ptVars[j+25].c_str());
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+25]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+25]->Draw();
      //h_ncs_2_closure[3-i][j+25]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+25]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_350_400 = new TCanvas("c_ncs_350_400","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_350_400->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_350_400->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+30]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+30]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetTitle(ptVars[j+30].c_str());
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+30]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+30]->Draw();
      //h_ncs_2_closure[3-i][j+30]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+30]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_400_450 = new TCanvas("c_ncs_400_450","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_400_450->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_400_450->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+35]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+35]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetTitle(ptVars[j+35].c_str());
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+35]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+35]->Draw();
     // h_ncs_2_closure[3-i][j+35]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+35]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }
/*
  TCanvas *c_ncs_450_500 = new TCanvas("c_ncs_450_500","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_450_500->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_450_500->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+40]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+40]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetTitle(ptVars[j+40].c_str());
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+40]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+40]->Draw();
     // h_ncs_2_closure[3-i][j+40]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+40]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_500_550 = new TCanvas("c_ncs_500_550","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_500_550->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_500_550->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+45]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+45]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetTitle(ptVars[j+45].c_str());
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+45]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+45]->Draw();
      //h_ncs_2_closure[3-i][j+45]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+45]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_ncs_550_600 = new TCanvas("c_ncs_550_600","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_ncs_550_600->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_ncs_550_600->cd(i+(4*j)+1);
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_closure[3-i][j+50]->GetXaxis()->SetLabelSize(0.12);
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_closure[3-i][j+50]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetTitle(ptVars[j+50].c_str());
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_closure[3-i][j+50]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_closure[3-i][j+50]->Draw();
      //h_ncs_2_closure[3-i][j+45]->SetLineColor(8);
      //h_ncs_2_closure[3-i][j+45]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }
*/

/*
  TCanvas *c_res_50_100 = new TCanvas("c_res_50_100","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_50_100->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_50_100->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j]->GetBinContent(h_ncs_2_res_corr[3-i][j]->GetMaximumBin()));
      //h_ncs_2_res[3-i][j]->GetXaxis()->SetRangeUser(0.,2.);
      h_ncs_2_res[3-i][j]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j]->Draw("same");
      
      h_ncs_2_res[3-i][j]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j]->Draw("same");
      h_ncs_2_res[3-i][j]->Fit(f_ncs_res[3-i][j],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j]->Fit(f_ncs_res_2[3-i][j],"Q M R","",(f_ncs_res[3-i][j]->GetParameter(1))*0.7,(f_ncs_res[3-i][j]->GetParameter(1))*1.3);
      par_res[3-i][j] = fabs(f_ncs_res_2[3-i][j]->GetParameter(2));
      par_res_err[3-i][j] = fabs(f_ncs_res_2[3-i][j]->GetParError(2));
      h_ncs_2_res_corr[3-i][j]->Fit(f_ncs_res_corr[3-i][j],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j]->Fit(f_ncs_res_corr_2[3-i][j],"Q M R","sames",(f_ncs_res_corr[3-i][j]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j]->GetParameter(1))*1.3);
      par_res_corr[3-i][j] = fabs(f_ncs_res_corr_2[3-i][j]->GetParameter(2));
      par_res_corr_err[3-i][j] = fabs(f_ncs_res_corr_2[3-i][j]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_100_150 = new TCanvas("c_res_100_150","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_100_150->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_100_150->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+5]->GetBinContent(h_ncs_2_res_corr[3-i][j+5]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+5]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+5]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+5]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+5]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+5]->Draw("same");
      
      h_ncs_2_res[3-i][j+5]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+5]->Draw("same");
      h_ncs_2_res[3-i][j+5]->Fit(f_ncs_res[3-i][j+5],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+5]->Fit(f_ncs_res_2[3-i][j+5],"Q M R","",(f_ncs_res[3-i][j+5]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+5]->GetParameter(1))*1.3);
      par_res[3-i][j+5] = fabs(f_ncs_res_2[3-i][j+5]->GetParameter(2));
      par_res_err[3-i][j+5] = fabs(f_ncs_res_2[3-i][j+5]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+5]->Fit(f_ncs_res_corr[3-i][j+5],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+5]->Fit(f_ncs_res_corr_2[3-i][j+5],"Q M R","sames",(f_ncs_res_corr[3-i][j+5]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+5]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+5] = fabs(f_ncs_res_corr_2[3-i][j+5]->GetParameter(2));
      par_res_corr_err[3-i][j+5] = fabs(f_ncs_res_corr_2[3-i][j+5]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+5].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_150_200 = new TCanvas("c_res_150_200","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_150_200->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_150_200->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+10]->GetBinContent(h_ncs_2_res_corr[3-i][j+10]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+10]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+10]->SetTitle(centVars[3-i].c_str());
      //h_ncs_2_res[3-i][j+10]->GetYaxis()->SetTitle(ptVars[j+10].c_str());
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+10]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res[3-i][j+10]->GetXaxis()->SetTitle("recopT/genpT");
      h_ncs_2_res_corr[3-i][j+10]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+10]->Draw("same");
      
      h_ncs_2_res[3-i][j+10]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+10]->Draw("same");
      h_ncs_2_res[3-i][j+10]->Fit(f_ncs_res[3-i][j+10],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+10]->Fit(f_ncs_res_2[3-i][j+10],"Q M R","",(f_ncs_res[3-i][j+10]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+10]->GetParameter(1))*1.3);
      par_res[3-i][j+10] = fabs(f_ncs_res_2[3-i][j+10]->GetParameter(2));
      par_res_err[3-i][j+10] = fabs(f_ncs_res_2[3-i][j+10]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+10]->Fit(f_ncs_res_corr[3-i][j+10],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+10]->Fit(f_ncs_res_corr_2[3-i][j+10],"Q M R","sames",(f_ncs_res_corr[3-i][j+10]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+10]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+10] = fabs(f_ncs_res_corr_2[3-i][j+10]->GetParameter(2));
      par_res_corr_err[3-i][j+10] = fabs(f_ncs_res_corr_2[3-i][j+10]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+10].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_200_250 = new TCanvas("c_res_200_250","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_200_250->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_200_250->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+15]->GetBinContent(h_ncs_2_res_corr[3-i][j+15]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+15]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+15]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+15]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+15]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+15]->Draw("same");
      
      h_ncs_2_res[3-i][j+15]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+15]->Draw("same");
      h_ncs_2_res[3-i][j+15]->Fit(f_ncs_res[3-i][j+15],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+15]->Fit(f_ncs_res_2[3-i][j+15],"Q M R","",(f_ncs_res[3-i][j+15]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+15]->GetParameter(1))*1.3);
      par_res[3-i][j+15] = fabs(f_ncs_res_2[3-i][j+15]->GetParameter(2));
      par_res_err[3-i][j+15] = fabs(f_ncs_res_2[3-i][j+15]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+15]->Fit(f_ncs_res_corr[3-i][j+15],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+15]->Fit(f_ncs_res_corr_2[3-i][j+15],"Q M R","sames",(f_ncs_res_corr[3-i][j+15]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+15]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+15] = fabs(f_ncs_res_corr_2[3-i][j+15]->GetParameter(2));
      par_res_corr_err[3-i][j+15] = fabs(f_ncs_res_corr_2[3-i][j+15]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+15].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_250_300 = new TCanvas("c_res_250_300","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_250_300->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_250_300->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+20]->GetBinContent(h_ncs_2_res_corr[3-i][j+20]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+20]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+20]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+20]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+20]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+20]->Draw("same");
      
      h_ncs_2_res[3-i][j+20]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+20]->Draw("same");
      h_ncs_2_res[3-i][j+20]->Fit(f_ncs_res[3-i][j+20],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+20]->Fit(f_ncs_res_2[3-i][j+20],"Q M R","",(f_ncs_res[3-i][j+20]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+20]->GetParameter(1))*1.3);
      par_res[3-i][j+20] = fabs(f_ncs_res_2[3-i][j+20]->GetParameter(2));
      par_res_err[3-i][j+20] = fabs(f_ncs_res_2[3-i][j+20]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+20]->Fit(f_ncs_res_corr[3-i][j+20],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+20]->Fit(f_ncs_res_corr_2[3-i][j+20],"Q M R","sames",(f_ncs_res_corr[3-i][j+20]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+20]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+20] = fabs(f_ncs_res_corr_2[3-i][j+20]->GetParameter(2));
      par_res_corr_err[3-i][j+20] = fabs(f_ncs_res_corr_2[3-i][j+20]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+20].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_300_350 = new TCanvas("c_res_300_350","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_300_350->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_300_350->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+25]->GetBinContent(h_ncs_2_res_corr[3-i][j+25]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+25]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+25]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+25]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+25]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+25]->Draw("same");
      
      h_ncs_2_res[3-i][j+25]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+25]->Draw("same");
      h_ncs_2_res[3-i][j+25]->Fit(f_ncs_res[3-i][j+25],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+25]->Fit(f_ncs_res_2[3-i][j+25],"Q M R","",(f_ncs_res[3-i][j+25]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+25]->GetParameter(1))*1.3);
      par_res[3-i][j+25] = fabs(f_ncs_res_2[3-i][j+25]->GetParameter(2));
      par_res_err[3-i][j+25] = fabs(f_ncs_res_2[3-i][j+25]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+25]->Fit(f_ncs_res_corr[3-i][j+25],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+25]->Fit(f_ncs_res_corr_2[3-i][j+25],"Q M R","sames",(f_ncs_res_corr[3-i][j+25]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+25]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+25] = fabs(f_ncs_res_corr_2[3-i][j+25]->GetParameter(2));
      par_res_corr_err[3-i][j+25] = fabs(f_ncs_res_corr_2[3-i][j+25]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+25].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_350_400 = new TCanvas("c_res_350_400","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_350_400->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_350_400->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+30]->GetBinContent(h_ncs_2_res_corr[3-i][j+30]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+30]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+30]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+30]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+30]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+30]->Draw("same");
      
      h_ncs_2_res[3-i][j+30]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+30]->Draw("same");
      h_ncs_2_res[3-i][j+30]->Fit(f_ncs_res[3-i][j+30],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+30]->Fit(f_ncs_res_2[3-i][j+30],"Q M R","",(f_ncs_res[3-i][j+30]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+30]->GetParameter(1))*1.3);
      par_res[3-i][j+30] = fabs(f_ncs_res_2[3-i][j+30]->GetParameter(2));
      par_res_err[3-i][j+30] = fabs(f_ncs_res_2[3-i][j+30]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+30]->Fit(f_ncs_res_corr[3-i][j+30],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+30]->Fit(f_ncs_res_corr_2[3-i][j+30],"Q M R","sames",(f_ncs_res_corr[3-i][j+30]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+30]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+30] = fabs(f_ncs_res_corr_2[3-i][j+30]->GetParameter(2));
      par_res_corr_err[3-i][j+30] = fabs(f_ncs_res_corr_2[3-i][j+30]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+30].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_400_450 = new TCanvas("c_res_400_450","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_400_450->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_400_450->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+35]->GetBinContent(h_ncs_2_res_corr[3-i][j+35]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+35]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+35]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+35]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+35]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+35]->Draw("same");
      
      h_ncs_2_res[3-i][j+35]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+35]->Draw("same");
      h_ncs_2_res[3-i][j+35]->Fit(f_ncs_res[3-i][j+35],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+35]->Fit(f_ncs_res_2[3-i][j+35],"Q M R","",(f_ncs_res[3-i][j+35]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+35]->GetParameter(1))*1.3);
      par_res[3-i][j+35] = fabs(f_ncs_res_2[3-i][j+35]->GetParameter(2));
      par_res_err[3-i][j+35] = fabs(f_ncs_res_2[3-i][j+35]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+35]->Fit(f_ncs_res_corr[3-i][j+35],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+35]->Fit(f_ncs_res_corr_2[3-i][j+35],"Q M R","sames",(f_ncs_res_corr[3-i][j+35]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+35]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+35] = fabs(f_ncs_res_corr_2[3-i][j+35]->GetParameter(2));
      par_res_corr_err[3-i][j+35] = fabs(f_ncs_res_corr_2[3-i][j+35]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+35].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_450_500 = new TCanvas("c_res_450_500","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_450_500->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_450_500->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+40]->GetBinContent(h_ncs_2_res_corr[3-i][j+40]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+40]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+40]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+40]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+40]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+40]->Draw("same");
      
      h_ncs_2_res[3-i][j+40]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+40]->Draw("same");
      h_ncs_2_res[3-i][j+40]->Fit(f_ncs_res[3-i][j+40],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+40]->Fit(f_ncs_res_2[3-i][j+40],"Q M R","",(f_ncs_res[3-i][j+40]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+40]->GetParameter(1))*1.3);
      par_res[3-i][j+40] = fabs(f_ncs_res_2[3-i][j+40]->GetParameter(2));
      par_res_err[3-i][j+40] = fabs(f_ncs_res_2[3-i][j+40]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+40]->Fit(f_ncs_res_corr[3-i][j+40],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+40]->Fit(f_ncs_res_corr_2[3-i][j+40],"Q M R","sames",(f_ncs_res_corr[3-i][j+40]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+40]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+40] = fabs(f_ncs_res_corr_2[3-i][j+40]->GetParameter(2));
      par_res_corr_err[3-i][j+40] = fabs(f_ncs_res_corr_2[3-i][j+40]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+40].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_500_550 = new TCanvas("c_res_500_550","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_500_550->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_500_550->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+45]->GetBinContent(h_ncs_2_res_corr[3-i][j+45]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+45]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+45]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+45]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+45]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+45]->Draw("same");
      
      h_ncs_2_res[3-i][j+45]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+45]->Draw("same");
      h_ncs_2_res[3-i][j+45]->Fit(f_ncs_res[3-i][j+45],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+45]->Fit(f_ncs_res_2[3-i][j+45],"Q M R","",(f_ncs_res[3-i][j+45]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+45]->GetParameter(1))*1.3);
      par_res[3-i][j+45] = fabs(f_ncs_res_2[3-i][j+45]->GetParameter(2));
      par_res_err[3-i][j+45] = fabs(f_ncs_res_2[3-i][j+45]->GetParError(2));
      h_ncs_2_res_corr[3-i][j+45]->Fit(f_ncs_res_corr[3-i][j+45],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+45]->Fit(f_ncs_res_corr_2[3-i][j+45],"Q M R","sames",(f_ncs_res_corr[3-i][j+45]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+45]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+45] = fabs(f_ncs_res_corr_2[3-i][j+45]->GetParameter(2));
      par_res_corr_err[3-i][j+45] = fabs(f_ncs_res_corr_2[3-i][j+45]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+45].c_str(), "");
      legend ->Draw("same");
    }
  }

  TCanvas *c_res_550_600 = new TCanvas("c_res_550_600","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_550_600->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_550_600->cd(i+(4*j)+1);
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetRangeUser(0.,1.2*h_ncs_2_res_corr[3-i][j+50]->GetBinContent(h_ncs_2_res_corr[3-i][j+50]->GetMaximumBin()));
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res[3-i][j+50]->GetXaxis()->SetLabelSize(0.08);
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res[3-i][j+50]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetTitle("Entries");
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetTitleSize(0.06);
      h_ncs_2_res[3-i][j+50]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+50]->SetLineColor(kRed);
      h_ncs_2_res_corr[3-i][j+50]->Draw("same");
      h_ncs_2_res[3-i][j+50]->SetLineColor(kBlue);
      h_ncs_2_res[3-i][j+50]->Draw("same");
      h_ncs_2_res[3-i][j+50]->Fit(f_ncs_res[3-i][j+50],"Q M R N","",0.,2.);
      h_ncs_2_res[3-i][j+50]->Fit(f_ncs_res_2[3-i][j+50],"Q M R","",(f_ncs_res[3-i][j+50]->GetParameter(1))*0.7,(f_ncs_res[3-i][j+50]->GetParameter(1))*1.3);
      par_res[3-i][j+50] = fabs(f_ncs_res_2[3-i][j+50]->GetParameter(2));
      par_res_err[3-i][j+50] = fabs(f_ncs_res_2[3-i][j+50]->GetParError(2));
      
      h_ncs_2_res_corr[3-i][j+50]->Fit(f_ncs_res_corr[3-i][j+50],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+50]->Fit(f_ncs_res_corr_2[3-i][j+50],"Q M R","sames",(f_ncs_res_corr[3-i][j+50]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+50]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+50] = fabs(f_ncs_res_corr_2[3-i][j+50]->GetParameter(2));
      par_res_corr_err[3-i][j+50] = fabs(f_ncs_res_corr_2[3-i][j+50]->GetParError(2));
      TLegend *legend = new TLegend(0.15,0.5,0.4,0.75);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry((TObject*)0, ptVars[j+50].c_str(), "");
      legend ->Draw("same");
    }
  }
*/
/*
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_res_corr_50_100 = new TCanvas("c_res_corr_50_100","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_50_100->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_50_100->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j]->GetYaxis()->SetRangeUser(0.5,1.5);
      //h_ncs_2_res[3-i][j]->GetXaxis()->SetRangeUser(0.,2.);
      h_ncs_2_res_corr[3-i][j]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j]->GetXaxis()->SetLabelSize(0.09);
      h_ncs_2_res_corr[3-i][j]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j]->GetYaxis()->SetTitle(ptVars[j].c_str());
      h_ncs_2_res_corr[3-i][j]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j]->Draw();
      h_ncs_2_res_corr[3-i][j]->Fit(f_ncs_res_corr[3-i][j],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j]->Fit(f_ncs_res_corr_2[3-i][j],"Q M R","",(f_ncs_res_corr[3-i][j]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j]->GetParameter(1))*1.3);
      par_res_corr[3-i][j] = fabs(f_ncs_res_corr_2[3-i][j]->GetParameter(2));
      par_res_corr_err[3-i][j] = fabs(f_ncs_res_corr_2[3-i][j]->GetParError(2));
      h_ncs_2_res_q[3-i][j]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j]->Draw("same");
      h_ncs_2_res_g[3-i][j]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j]->Draw("same");
      //h_ncs_2_res[3-i][j]->SetLineColor(8);
      //h_ncs_2_res[3-i][j]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_100_150 = new TCanvas("c_res_corr_100_150","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_100_150->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_100_150->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+5]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+5]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+5]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+5]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+5]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+5]->GetYaxis()->SetTitle(ptVars[j+5].c_str());
      h_ncs_2_res_corr[3-i][j+5]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+5]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+5]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+5]->Draw();
      h_ncs_2_res_corr[3-i][j+5]->Fit(f_ncs_res_corr[3-i][j+5],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+5]->Fit(f_ncs_res_corr_2[3-i][j+5],"Q M R","",(f_ncs_res_corr[3-i][j+5]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+5]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+5] = fabs(f_ncs_res_corr_2[3-i][j+5]->GetParameter(2));
      par_res_corr_err[3-i][j+5] = fabs(f_ncs_res_corr_2[3-i][j+5]->GetParError(2));
      h_ncs_2_res_q[3-i][j+5]->SetLineColor(kBlue);
      //h_ncs_2_res_q[3-i][j+5]->Draw("same");
      h_ncs_2_res_g[3-i][j+5]->SetLineColor(kRed);
      //h_ncs_2_res_g[3-i][j+5]->Draw("same");
      //h_ncs_2_res[3-i][j+5]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+5]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_150_200 = new TCanvas("c_res_corr_150_200","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_150_200->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_150_200->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+10]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+10]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+10]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+10]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+10]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+10]->GetYaxis()->SetTitle(ptVars[j+10].c_str());
      h_ncs_2_res_corr[3-i][j+10]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+10]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+10]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+10]->Draw();
      h_ncs_2_res_corr[3-i][j+10]->Fit(f_ncs_res_corr[3-i][j+10],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+10]->Fit(f_ncs_res_corr_2[3-i][j+10],"Q M R","",(f_ncs_res_corr[3-i][j+10]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+10]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+10] = fabs(f_ncs_res_corr_2[3-i][j+10]->GetParameter(2));
      par_res_corr_err[3-i][j+10] = fabs(f_ncs_res_corr_2[3-i][j+10]->GetParError(2));
      h_ncs_2_res_q[3-i][j+10]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+10]->Draw("same");
      h_ncs_2_res_g[3-i][j+10]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+10]->Draw("same");
      // h_ncs_2_res[3-i][j+10]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+10]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_200_250 = new TCanvas("c_res_corr_200_250","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_200_250->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_200_250->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+15]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+15]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+15]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+15]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+15]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+15]->GetYaxis()->SetTitle(ptVars[j+15].c_str());
      h_ncs_2_res_corr[3-i][j+15]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+15]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+15]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+15]->Draw();
      h_ncs_2_res_corr[3-i][j+15]->Fit(f_ncs_res_corr[3-i][j+15],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+15]->Fit(f_ncs_res_corr_2[3-i][j+15],"Q M R","",(f_ncs_res_corr[3-i][j+15]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+15]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+15] = fabs(f_ncs_res_corr_2[3-i][j+15]->GetParameter(2));
      par_res_corr_err[3-i][j+15] = fabs(f_ncs_res_corr_2[3-i][j+15]->GetParError(2));
      h_ncs_2_res_q[3-i][j+15]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+15]->Draw("same");
      h_ncs_2_res_g[3-i][j+15]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+15]->Draw("same");
      //h_ncs_2_res[3-i][j+15]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+15]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_250_300 = new TCanvas("c_res_corr_250_300","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_250_300->Divide(4,5);

  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_250_300->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+20]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+20]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+20]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+20]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+20]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+20]->GetYaxis()->SetTitle(ptVars[j+20].c_str());
      h_ncs_2_res_corr[3-i][j+20]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+20]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+20]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+20]->Draw();
      h_ncs_2_res_corr[3-i][j+20]->Fit(f_ncs_res_corr[3-i][j+20],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+20]->Fit(f_ncs_res_corr_2[3-i][j+20],"Q M R","",(f_ncs_res_corr[3-i][j+20]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+20]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+20] = fabs(f_ncs_res_corr_2[3-i][j+20]->GetParameter(2));
      par_res_corr_err[3-i][j+20] = fabs(f_ncs_res_corr_2[3-i][j+20]->GetParError(2));
      h_ncs_2_res_q[3-i][j+20]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+20]->Draw("same");
      h_ncs_2_res_g[3-i][j+20]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+20]->Draw("same");//h_ncs_2_res[3-i][j+20]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+20]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_300_350 = new TCanvas("c_res_corr_300_350","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_300_350->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_300_350->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+25]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+25]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+25]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+25]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+25]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+25]->GetYaxis()->SetTitle(ptVars[j+25].c_str());
      h_ncs_2_res_corr[3-i][j+25]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+25]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+25]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+25]->Draw();
      h_ncs_2_res_corr[3-i][j+25]->Fit(f_ncs_res_corr[3-i][j+25],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+25]->Fit(f_ncs_res_corr_2[3-i][j+25],"Q M R","",(f_ncs_res_corr[3-i][j+25]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+25]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+25] = fabs(f_ncs_res_corr_2[3-i][j+25]->GetParameter(2));
      par_res_corr_err[3-i][j+25] = fabs(f_ncs_res_corr_2[3-i][j+25]->GetParError(2));
      h_ncs_2_res_q[3-i][j+25]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+25]->Draw("same");
      h_ncs_2_res_g[3-i][j+25]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+25]->Draw("same");
      //h_ncs_2_res[3-i][j+25]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+25]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_350_400 = new TCanvas("c_res_corr_350_400","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_350_400->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_350_400->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+30]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+30]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+30]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+30]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+30]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+30]->GetYaxis()->SetTitle(ptVars[j+30].c_str());
      h_ncs_2_res_corr[3-i][j+30]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+30]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+30]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+30]->Draw();
      h_ncs_2_res_corr[3-i][j+30]->Fit(f_ncs_res_corr[3-i][j+30],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+30]->Fit(f_ncs_res_corr_2[3-i][j+30],"Q M R","",(f_ncs_res_corr[3-i][j+30]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+30]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+30] = fabs(f_ncs_res_corr_2[3-i][j+30]->GetParameter(2));
      par_res_corr_err[3-i][j+30] = fabs(f_ncs_res_corr_2[3-i][j+30]->GetParError(2));
      h_ncs_2_res_q[3-i][j+30]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+30]->Draw("same");
      h_ncs_2_res_g[3-i][j+30]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+30]->Draw("same");
      //h_ncs_2_res[3-i][j+30]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+30]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_400_450 = new TCanvas("c_res_corr_400_450","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_400_450->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_400_450->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+35]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+35]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+35]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+35]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+35]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+35]->GetYaxis()->SetTitle(ptVars[j+35].c_str());
      h_ncs_2_res_corr[3-i][j+35]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+35]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+35]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+35]->Draw();
      h_ncs_2_res_corr[3-i][j+35]->Fit(f_ncs_res_corr[3-i][j+35],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+35]->Fit(f_ncs_res_corr_2[3-i][j+35],"Q M R","",(f_ncs_res_corr[3-i][j+35]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+35]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+35] = fabs(f_ncs_res_corr_2[3-i][j+35]->GetParameter(2));
      par_res_corr_err[3-i][j+35] = fabs(f_ncs_res_corr_2[3-i][j+35]->GetParError(2));
      h_ncs_2_res_q[3-i][j+35]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+35]->Draw("same");
      h_ncs_2_res_g[3-i][j+35]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+35]->Draw("same");
      // h_ncs_2_res[3-i][j+35]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+35]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_450_500 = new TCanvas("c_res_corr_450_500","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_450_500->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_450_500->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+40]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+40]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+40]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+40]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+40]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+40]->GetYaxis()->SetTitle(ptVars[j+40].c_str());
      h_ncs_2_res_corr[3-i][j+40]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+40]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+40]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+40]->Draw();
      h_ncs_2_res_corr[3-i][j+40]->Fit(f_ncs_res_corr[3-i][j+40],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+40]->Fit(f_ncs_res_corr_2[3-i][j+40],"Q M R","",(f_ncs_res_corr[3-i][j+40]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+40]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+40] = fabs(f_ncs_res_corr_2[3-i][j+40]->GetParameter(2));
      par_res_corr_err[3-i][j+40] = fabs(f_ncs_res_corr_2[3-i][j+40]->GetParError(2));
      h_ncs_2_res_q[3-i][j+40]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+40]->Draw("same");
      h_ncs_2_res_g[3-i][j+40]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+40]->Draw("same");
      // h_ncs_2_res[3-i][j+40]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+40]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_500_550 = new TCanvas("c_res_corr_500_550","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_500_550->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_500_550->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+45]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+45]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+45]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+45]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+45]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+45]->GetYaxis()->SetTitle(ptVars[j+45].c_str());
      h_ncs_2_res_corr[3-i][j+45]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+45]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+45]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+45]->Draw();
      h_ncs_2_res_corr[3-i][j+45]->Fit(f_ncs_res_corr[3-i][j+45],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+45]->Fit(f_ncs_res_corr_2[3-i][j+45],"Q M R","",(f_ncs_res_corr[3-i][j+45]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+45]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+45] = fabs(f_ncs_res_corr_2[3-i][j+45]->GetParameter(2));
      par_res_corr_err[3-i][j+45] = fabs(f_ncs_res_corr_2[3-i][j+45]->GetParError(2));
      h_ncs_2_res_q[3-i][j+45]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+45]->Draw("same");
      h_ncs_2_res_g[3-i][j+45]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+45]->Draw("same");
      //h_ncs_2_res[3-i][j+45]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+45]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }

  TCanvas *c_res_corr_550_600 = new TCanvas("c_res_corr_550_600","",900,1200);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c_res_corr_550_600->Divide(4,5);
  
  for(int j=0; j<5; j++){  
    for(int i=0; i<4; i++){
      c_res_corr_550_600->cd(i+(4*j)+1);
      //h_ncs_2_res[3-i][j+50]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_ncs_2_res_corr[3-i][j+50]->GetYaxis()->SetLabelSize(0.07);
      h_ncs_2_res_corr[3-i][j+50]->GetXaxis()->SetLabelSize(0.10);
      h_ncs_2_res_corr[3-i][j+50]->GetYaxis()->SetNdivisions(5);
      h_ncs_2_res_corr[3-i][j+50]->SetTitle(centVars[3-i].c_str());
      h_ncs_2_res_corr[3-i][j+50]->GetYaxis()->SetTitle(ptVars[j+50].c_str());
      h_ncs_2_res_corr[3-i][j+50]->GetYaxis()->SetTitleSize(0.13);
      h_ncs_2_res_corr[3-i][j+50]->GetYaxis()->SetTitleOffset(0.30);
      h_ncs_2_res_corr[3-i][j+50]->SetLineColor(kBlack);
      h_ncs_2_res_corr[3-i][j+50]->Draw();
      h_ncs_2_res_corr[3-i][j+50]->Fit(f_ncs_res_corr[3-i][j+50],"Q M R N","",0.,2.);
      h_ncs_2_res_corr[3-i][j+50]->Fit(f_ncs_res_corr_2[3-i][j+50],"Q M R","",(f_ncs_res_corr[3-i][j+50]->GetParameter(1))*0.7,(f_ncs_res_corr[3-i][j+50]->GetParameter(1))*1.3);
      par_res_corr[3-i][j+50] = fabs(f_ncs_res_corr_2[3-i][j+50]->GetParameter(2));
      par_res_corr_err[3-i][j+50] = fabs(f_ncs_res_corr_2[3-i][j+50]->GetParError(2));
      h_ncs_2_res_q[3-i][j+50]->SetLineColor(kBlue);
      h_ncs_2_res_q[3-i][j+50]->Draw("same");
      h_ncs_2_res_g[3-i][j+50]->SetLineColor(kRed);
      h_ncs_2_res_g[3-i][j+50]->Draw("same");
      //h_ncs_2_res[3-i][j+45]->SetLineColor(8);
      //h_ncs_2_res[3-i][j+45]->Draw("same");   
      //l1[3-i][j] = new TLatex(130,1.3,centVars[3-i].c_str());
      //l1[3-i][j]->SetTextSize(0.08);
      //l1[3-i][j]->Draw("same");
    }
  }
  
*/

/////////////////////////////////////////////////////////////////////////////////

  ////////////////////  mean and resolution  //////////////////////

/////////////////////////////////////////////////////////////////////////////////

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin3=0; ibin3<nptBins; ibin3++){
      h_res[ibin]->SetBinContent(ibin3+1,par_res[ibin][ibin3]);
      h_res[ibin]->SetBinError(ibin3+1,par_res_err[ibin][ibin3]);

      h_res_corr[ibin]->SetBinContent(ibin3+1,par_res_corr[ibin][ibin3]);
      h_res_corr[ibin]->SetBinError(ibin3+1,par_res_corr_err[ibin][ibin3]);

      h_mean[ibin]->SetBinContent(ibin3+1,par_mean[ibin][ibin3]);
      h_mean[ibin]->SetBinError(ibin3+1,par_mean_err[ibin][ibin3]);

      h_mean_corr[ibin]->SetBinContent(ibin3+1,par_mean_corr[ibin][ibin3]);
      h_mean_corr[ibin]->SetBinError(ibin3+1,par_mean_corr_err[ibin][ibin3]);  
    }
  }

  TCanvas *c_res = new TCanvas("c_res","",1200,240);
  gStyle->SetOptFit(1);
  c_res->Divide(5,1);
  for(int i=0; i<4; i++){
    c_res->cd(i+2);
    h_res[3-i]->Rebin(4);
    h_res[3-i]->Scale(1./4.);
    h_res_corr[3-i]->Rebin(4);
    h_res_corr[3-i]->Scale(1./4.);
    h_res[3-i]->GetXaxis()->SetRangeUser(100.,480.);
    h_res[3-i]->GetYaxis()->SetRangeUser(0.05,0.25);
    h_res[3-i]->SetMarkerColor(4);
    h_res[3-i]->SetMarkerSize(0.5);
    h_res[3-i]->SetMarkerStyle(21);
    h_res_corr[3-i]->SetMarkerSize(0.5);
    h_res_corr[3-i]->SetMarkerStyle(21);
    h_res[3-i]->SetTitle(centVars[3-i].c_str());
    h_res[3-i]->GetYaxis()->SetTitle("#sigma (recopT/genpT)");
    h_res[3-i]->GetYaxis()->SetLabelSize(0.05);
    h_res[3-i]->GetXaxis()->SetTitle("gen pT");
    h_res[3-i]->SetMarkerColor(kBlack);
    h_res[3-i]->SetLineColor(kBlack);
    h_res[3-i]->Draw("same");
    h_res_corr[3-i]->SetMarkerColor(kRed);
    h_res_corr[3-i]->SetLineColor(kRed);
    h_res_corr[3-i]->Draw("same");
    
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

  TCanvas *c_mean = new TCanvas("c_mean","",1200,300);
  gStyle->SetOptFit(1);
  c_mean->Divide(4,1);
  for(int i=0; i<4; i++){
    c_mean->cd(i+1);
    h_mean[3-i]->Rebin(3);
    h_mean[3-i]->Scale(1./3.);
    h_mean_corr[3-i]->Rebin(3);
    h_mean_corr[3-i]->Scale(1./3.);
    h_mean[3-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_mean[3-i]->GetYaxis()->SetRangeUser(0.92,1.08);
    h_mean[3-i]->SetMarkerColor(4);
    h_mean[3-i]->SetMarkerSize(0.5);
    h_mean[3-i]->SetMarkerStyle(21);
    h_mean_corr[3-i]->SetMarkerSize(0.5);
    h_mean_corr[3-i]->SetMarkerStyle(21);
    h_mean[3-i]->SetTitle(centVars[3-i].c_str());
    h_mean[3-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_mean[3-i]->GetYaxis()->SetLabelSize(0.05);
    h_mean[3-i]->GetXaxis()->SetTitle("gen pT");
    h_mean[3-i]->SetMarkerColor(kBlack);
    h_mean[3-i]->SetLineColor(kBlack);
    h_mean[3-i]->Draw("same");
    h_mean_corr[3-i]->SetMarkerColor(kRed);
    h_mean_corr[3-i]->SetLineColor(kRed);
    h_mean_corr[3-i]->Draw("same");
    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");
  
    if(i==3){
      TLegend *leg_mean = new TLegend(0.5,0.8,0.99,0.99);
      leg_mean->SetLineColor(kWhite);
      leg_mean->AddEntry(h_mean[0][0], "uncorr reco jets", "lepf");
      leg_mean->AddEntry(h_mean_corr[0][0], "corr reco jets", "lepf");
      leg_mean->Draw("same");
    }
  }  

////// a0 and a1 for ncs pt > 2 //////////////

  TLegend *leg = new TLegend(0.5,0.8,0.99,0.99);
  leg->SetLineColor(kWhite);
  leg->AddEntry((TObject*)0, "nCS p_{T} > 2", "");
  //leg->AddEntry((TObject*)0, "|eta| < 1.6", "");

  TCanvas *c_param = new TCanvas("c_param","",1200,600);
  gStyle->SetOptFit(0);
  c_param->Divide(4,2);
  for(int i=0; i<4; i++){
    c_param->cd(i+1);
    gr_p0_param_2[3-i]->GetXaxis()->SetRangeUser(60.,500.);
    gr_p0_param_2[3-i]->GetYaxis()->SetRangeUser(7.,20.);
    //gr_p0_param_2[3-i]->GetYaxis()->SetRangeUser(4.,11.);
    gr_p0_param_2[3-i]->SetMarkerColor(4);
    gr_p0_param_2[3-i]->SetMarkerSize(0.5);
    gr_p0_param_2[3-i]->SetMarkerStyle(21);
    gr_p0_param_2[3-i]->SetTitle(centVars[3-i].c_str());
    gr_p0_param_2[3-i]->GetYaxis()->SetTitle("a0");
    gr_p0_param_2[3-i]->GetYaxis()->SetLabelSize(0.05);
    gr_p0_param_2[3-i]->GetXaxis()->SetTitle("reco pT");
    gr_p0_param_2[3-i]->Draw("ap same");

    c_param->cd(i+5);
    gr_p1_param_2[3-i]->GetYaxis()->SetRangeUser(-0.04,0.);
    gr_p1_param_2[3-i]->GetXaxis()->SetRangeUser(60.,500.);
    gr_p1_param_2[3-i]->SetMarkerColor(4);
    gr_p1_param_2[3-i]->SetMarkerSize(0.5);
    gr_p1_param_2[3-i]->SetMarkerStyle(21);
    gr_p1_param_2[3-i]->SetTitle("");
    gr_p1_param_2[3-i]->GetYaxis()->SetTitle("a1");
    gr_p1_param_2[3-i]->GetYaxis()->SetLabelSize(0.03);
    gr_p1_param_2[3-i]->GetXaxis()->SetTitle("reco pT");
    gr_p1_param_2[3-i]->Draw("ap same");
  }
  leg->Draw("same");
/*
////// a0 and a1 for ncs pt > 1 //////////////
 
  TLegend *legp = new TLegend(0.5,0.8,0.99,0.99);
  legp->SetLineColor(kWhite);
  legp->AddEntry((TObject*)0, "nCS p_{T} > 1", ""); 

  TCanvas *c_param_1 = new TCanvas("c_param_1","",1200,600);
  gStyle->SetOptFit(1);
  c_param_1->Divide(4,2);
  for(int i=0; i<4; i++){
    c_param_1->cd(i+1);
    gr_p0_param_1[3-i]->GetXaxis()->SetRangeUser(80.,550.);
    gr_p0_param_1[3-i]->GetYaxis()->SetRangeUser(5.,25.);
    gr_p0_param_1[3-i]->SetMarkerColor(4);
    gr_p0_param_1[3-i]->SetMarkerSize(0.5);
    gr_p0_param_1[3-i]->SetMarkerStyle(21);
    gr_p0_param_1[3-i]->SetTitle(centVars[3-i].c_str());
    gr_p0_param_1[3-i]->GetYaxis()->SetTitle("a0");
    gr_p0_param_1[3-i]->GetYaxis()->SetLabelSize(0.05);
    gr_p0_param_1[3-i]->GetXaxis()->SetTitle("reco pT");
    gr_p0_param_1[3-i]->Draw("ap same");

    c_param_1->cd(i+5);
    gr_p1_param_1[3-i]->GetYaxis()->SetRangeUser(-0.02,0.);
    gr_p1_param_1[3-i]->GetXaxis()->SetRangeUser(80.,550.);
    gr_p1_param_1[3-i]->SetMarkerColor(4);
    gr_p1_param_1[3-i]->SetMarkerSize(0.5);
    gr_p1_param_1[3-i]->SetMarkerStyle(21);
    gr_p1_param_1[3-i]->SetTitle("");
    gr_p1_param_1[3-i]->GetYaxis()->SetTitle("a1");
    gr_p1_param_1[3-i]->GetYaxis()->SetLabelSize(0.03);
    gr_p1_param_1[3-i]->GetXaxis()->SetTitle("reco pT");
    gr_p1_param_1[3-i]->Draw("ap same");
  } 
*/
  /////////////////////

  //Kurt's histos

  //////////////////////

  const bool ispp = false;
  
  int ncentbins=4;
  if(ispp) ncentbins=1;
  
  TProfile *fpx[3][4];
  TProfile *fpxPre[3][4];
  TProfile *fpxiter[3][4];
  TProfile *fpxpp;
  TProfile *fpxprepp;
  TH2D *jtclosures[4];
  TH2D *jtgclosures[4];
  TH2D *jtqclosures[4];
  
  TH2D *jtpreclosures[4];
  TH2D *jtqpreclosures[4];
  TH2D *jtgpreclosures[4];

  TH2D *jtiterclosures[4];
  TH2D *jtqiterclosures[4];
  TH2D *jtgiterclosures[4];
  
  TH2D *jtclosurepp;
  TH2D *jtpreclosurepp;
  TFile *fin, *fin2, *fin3;
  if(ispp){ 
    fin = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/fullIterativeJFFClosures_closureTest_pp.root");
    fin2 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/fullIterativeJFFClosures_closureTest_pp.root");
  }
  else{
    fin = new TFile("/home/dhanush/Documents/JEC/local/kurt/fullIterativeJFFClosures_closureTest_withEtaPhi.root");
    //fin2 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/fullIterativeJFFClosures_closureTest_pp.root");
    //fin3 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/fullIterativeJFFClosures_closureTest_withRXdown.root");
  }
  for(int i=0; i<ncentbins; i++){
    /*jtclosures[i] = (TH2D*)fin2->Get(Form("fullClosure_cent%d",i))->Clone(Form("h_jt_closure_cent%d",i));
    jtqclosures[i] = (TH2D*)fin2->Get(Form("quarkClosure_cent%d",i))->Clone(Form("h_jt_closure_q_cent%d",i)); 
    jtgclosures[i] = (TH2D*)fin2->Get(Form("gluonClosure_cent%d",i))->Clone(Form("h_jt_closure_g_cent%d",i));*/

    jtclosures[i] = (TH2D*)fin->Get(Form("s2refClosure_cent%d_iter0",i))->Clone(Form("h_jt_closure_cent%d",i));                          
    jtqclosures[i] = (TH2D*)fin->Get(Form("s2QrefClosure_cent%d_iter0",i))->Clone(Form("h_jt_closure_q_cent%d",i)); 
    jtgclosures[i] = (TH2D*)fin->Get(Form("s2GrefClosure_cent%d_iter0",i))->Clone(Form("h_jt_closure_g_cent%d",i));
    fpx[0][i] = jtclosures[i]->ProfileX();
    fpx[1][i] = jtqclosures[i]->ProfileX();
    fpx[2][i] = jtgclosures[i]->ProfileX();

	jtpreclosures[i] = (TH2D*)fin->Get(Form("s1refClosure_cent%d",i))->Clone(Form("h_jt_preclosure_cent%d",i));
	jtqpreclosures[i] = (TH2D*)fin->Get(Form("s1QrefClosure_cent%d",i))->Clone(Form("h_jt_preclosure_q_cent%d",i));
	jtgpreclosures[i] = (TH2D*)fin->Get(Form("s1GrefClosure_cent%d",i))->Clone(Form("h_jt_preclosure_g_cent%d",i));
	fpxPre[0][i] = jtpreclosures[i]->ProfileX();
	fpxPre[1][i] = jtqpreclosures[i]->ProfileX();
	fpxPre[2][i] = jtgpreclosures[i]->ProfileX();

    jtiterclosures[i] = (TH2D*)fin->Get(Form("s3refClosure_cent%d",i))->Clone(Form("h_jt_iterclosure_cent%d",i));                          
    jtqiterclosures[i] = (TH2D*)fin->Get(Form("s3QrefClosure_cent%d",i))->Clone(Form("h_jt_iterclosure_q_cent%d",i)); 
    jtgiterclosures[i] = (TH2D*)fin->Get(Form("s3GrefClosure_cent%d",i))->Clone(Form("h_jt_iterclosure_g_cent%d",i));
    fpxiter[0][i] = jtiterclosures[i]->ProfileX();
    fpxiter[1][i] = jtqiterclosures[i]->ProfileX();
    fpxiter[2][i] = jtgiterclosures[i]->ProfileX();

  }



  TCanvas *c_closure = new TCanvas("c_closure","",1200,300);
  c_closure->Divide(4,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c_closure->cd(i+1);
/*
    h_jt_closure_ref_ncs1_px[3-i] = h_jt_closure_ref_ncs1[3-i]->ProfileX();
    h_jt_closure_ref_ncs1_px[3-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_ncs1_px[3-i]->GetXaxis()->SetRangeUser(80.,500.);
    h_jt_closure_ref_ncs1_px[3-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_ncs1_px[3-i]->GetYaxis()->SetTitle("closure");
    h_jt_closure_ref_ncs1_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_ncs1_px[3-i]->SetMarkerStyle(29);
    h_jt_closure_ref_ncs1_px[3-i]->SetTitle(centVars[3-i].c_str());
    h_jt_closure_ref_ncs1_px[3-i]->Draw("e1 same");
    h_jt_closure_q_ncs1_px[3-i] = h_jt_closure_q_ncs1[3-i]->ProfileX();
    h_jt_closure_q_ncs1_px[3-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_ncs1_px[3-i]->SetMarkerStyle(29);
    h_jt_closure_q_ncs1_px[3-i]->Draw("e1 same");
    h_jt_closure_g_ncs1_px[3-i] = h_jt_closure_g_ncs1[3-i]->ProfileX();
    h_jt_closure_g_ncs1_px[3-i]->SetMarkerColor(kRed);
    h_jt_closure_g_ncs1_px[3-i]->SetMarkerStyle(29);    
    h_jt_closure_g_ncs1_px[3-i]->Draw("e1 same");
*/
    
    h_jt_closure_ref_ncs2_px[3-i] = h_jt_closure_ref_ncs2[3-i]->ProfileX();
    h_jt_closure_ref_ncs2_px[3-i]->Rebin(2);
    h_jt_closure_ref_ncs2_px[3-i]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_jt_closure_ref_ncs2_px[3-i]->GetXaxis()->SetRangeUser(100.,500.);
    h_jt_closure_ref_ncs2_px[3-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_ncs2_px[3-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_ncs2_px[3-i]->SetTitle(centVars[3-i].c_str());
    h_jt_closure_ref_ncs2_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_ncs2_px[3-i]->SetLineColor(kBlack);
    h_jt_closure_ref_ncs2_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_ref_ncs2_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_ref_ncs2_px[3-i]->Draw("e1 same");
    //h_jt_closure_ref_ncs2_px[3-i]->Fit(f_closure_2[3-i],"Q M","",80.,500.);
    h_jt_closure_q_ncs2_px[3-i] = h_jt_closure_q_ncs2[3-i]->ProfileX();
    h_jt_closure_q_ncs2_px[3-i]->Rebin(2);
    h_jt_closure_q_ncs2_px[3-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_ncs2_px[3-i]->SetLineColor(kBlue);
    h_jt_closure_q_ncs2_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_q_ncs2_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_q_ncs2_px[3-i]->Draw("e1 same");
    h_jt_closure_g_ncs2_px[3-i] = h_jt_closure_g_ncs2[3-i]->ProfileX();
    h_jt_closure_g_ncs2_px[3-i]->Rebin(2);
    h_jt_closure_g_ncs2_px[3-i]->SetMarkerColor(kRed);
    h_jt_closure_g_ncs2_px[3-i]->SetLineColor(kRed);
    h_jt_closure_g_ncs2_px[3-i]->SetMarkerStyle(24);
    h_jt_closure_g_ncs2_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_g_ncs2_px[3-i]->Draw("e1 same");
    h_jt_closure_ref_nocorr_px[3-i] = h_jt_closure_ref_nocorr[3-i]->ProfileX();
    /*
    h_jt_closure_ref_nocorr_px[3-i]->GetYaxis()->SetRangeUser(0.92,1.08);
    h_jt_closure_ref_nocorr_px[3-i]->GetXaxis()->SetTitle("gen pT");
    h_jt_closure_ref_nocorr_px[3-i]->GetYaxis()->SetTitle("#mu (recopT/genpT)");
    h_jt_closure_ref_nocorr_px[3-i]->SetTitle(centVars[3-i].c_str());
    */
    h_jt_closure_ref_nocorr_px[3-i]->Rebin(2);
    //h_jt_closure_ref_nocorr_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_ref_nocorr_px[3-i]->SetLineColor(kBlack);
    //h_jt_closure_ref_nocorr_px[3-i]->SetMarkerStyle(7);
    //h_jt_closure_ref_nocorr_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_ref_nocorr_px[3-i]->GetXaxis()->SetRangeUser(90.,500.);
    h_jt_closure_ref_nocorr_px[3-i]->Draw("e0 same");
    h_jt_closure_q_nocorr_px[3-i] = h_jt_closure_q_nocorr[3-i]->ProfileX();
    h_jt_closure_q_nocorr_px[3-i]->Rebin(2);
    //h_jt_closure_q_nocorr_px[3-i]->SetMarkerColor(kBlue);
    h_jt_closure_q_nocorr_px[3-i]->SetLineColor(kBlue);
    //h_jt_closure_q_nocorr_px[3-i]->SetMarkerStyle(7);
    //h_jt_closure_q_nocorr_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_q_nocorr_px[3-i]->Draw("e0 same");
    h_jt_closure_g_nocorr_px[3-i] = h_jt_closure_g_nocorr[3-i]->ProfileX();
    h_jt_closure_g_nocorr_px[3-i]->Rebin(2);
    //h_jt_closure_g_nocorr_px[3-i]->SetMarkerColor(kRed);
    h_jt_closure_g_nocorr_px[3-i]->SetLineColor(kRed);
    //h_jt_closure_g_nocorr_px[3-i]->SetMarkerStyle(7);
    //h_jt_closure_g_nocorr_px[3-i]->SetMarkerSize(0.7);
    h_jt_closure_g_nocorr_px[3-i]->Draw("e0 same");

    fpxiter[0][3-i]->Rebin(5);
    fpxiter[1][3-i]->Rebin(5);
    fpxiter[2][3-i]->Rebin(5);
    fpxiter[0][3-i]->SetLineColor(kBlack);
    fpxiter[0][3-i]->SetMarkerColor(kBlack);
    fpxiter[1][3-i]->SetLineColor(kBlue);
    fpxiter[1][3-i]->SetMarkerColor(kBlue);
    fpxiter[2][3-i]->SetLineColor(kRed);
    fpxiter[2][3-i]->SetMarkerColor(kRed);
    fpxiter[0][3-i]->SetMarkerStyle(31);
    fpxiter[1][3-i]->SetMarkerStyle(31);
    fpxiter[2][3-i]->SetMarkerStyle(31);
    fpxiter[0][3-i]->SetMarkerSize(1);
    fpxiter[1][3-i]->SetMarkerSize(1);
    fpxiter[2][3-i]->SetMarkerSize(1);
    //fpxiter[0][3-i]->Draw("e1 same");
    //fpxiter[1][3-i]->Draw("e1 same");
    //fpxiter[2][3-i]->Draw("e1 same");
/*
    fpxPre[0][3-i]->Rebin(5);
    fpxPre[1][3-i]->Rebin(5);
    fpxPre[2][3-i]->Rebin(5);
    fpxPre[0][3-i]->SetLineColor(kBlack);
    fpxPre[0][3-i]->SetMarkerColor(kBlack);
    fpxPre[1][3-i]->SetLineColor(kBlue);
    fpxPre[1][3-i]->SetMarkerColor(kBlue);
    fpxPre[2][3-i]->SetLineColor(kRed);
    fpxPre[2][3-i]->SetMarkerColor(kRed);
    fpxPre[0][3-i]->SetMarkerStyle(31);
    fpxPre[1][3-i]->SetMarkerStyle(31);
    fpxPre[2][3-i]->SetMarkerStyle(31);
    fpxPre[0][3-i]->SetMarkerSize(1);
    fpxPre[1][3-i]->SetMarkerSize(1);
    fpxPre[2][3-i]->SetMarkerSize(1);
    //fpxPre[0][3-i]->Draw("e1 same");
    //fpxPre[1][3-i]->Draw("e1 same");
    //fpxPre[2][3-i]->Draw("e1 same");
*/
    tl1->Draw("same");
    tl2->Draw("same");
    tl3->Draw("same");
    tl4->Draw("same");

    if(i==0){
      TLegend *leg3 = new TLegend(0.5,0.8,0.99,0.99);
      leg3->SetLineColor(0);
      leg3->SetFillColor(0);
      //leg1->AddEntry(h_jt_closure_ref_nocorr_px[0], "Pre-corr Inclusive", "lepf");
      leg3->AddEntry((TObject*)0, "nCS ID=1,4,5 cymbal", "");
      leg3->Draw("same");
    }

    if(i==1){
      TLegend *leg1 = new TLegend(0.5,0.8,0.99,0.99);
      leg1->SetLineColor(0);
      leg1->SetFillColor(0);
      //leg1->AddEntry(h_jt_closure_ref_nocorr_px[0], "Pre-corr Inclusive", "lepf");
      leg1->AddEntry(h_jt_closure_ref_ncs2_px[2], "Corrected Incl. Jets", "lepf");
      leg1->AddEntry(h_jt_closure_q_ncs2_px[2], "Corrected Quark Jets", "lepf");
      leg1->AddEntry(h_jt_closure_g_ncs2_px[2], "Corrected Gluon Jets", "lepf");
      leg1->Draw("same");
    }

    if(i==2){
      TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
      leg2->SetLineColor(0);
      leg2->SetFillColor(0);
      //leg1->AddEntry(h_jt_closure_ref_nocorr_px[0], "Pre-corr Inclusive", "lepf");
      leg2->AddEntry(h_jt_closure_ref_nocorr_px[1], "Pre-correction Incl. Jets", "lepf");
      leg2->AddEntry(h_jt_closure_q_nocorr_px[1], "Pre-correction Quark Jets", "lepf");
      leg2->AddEntry(h_jt_closure_g_nocorr_px[1], "Pre-correction Gluon Jets", "lepf");
      leg2->Draw("same");
    } 
  }

  
/*  
     TLegend *leg1 = new TLegend(0.5,0.8,0.99,0.99);
  leg1->SetLineColor(kWhite);
  //leg1->AddEntry(fpxiter[0][0], "Kurt", "lepf");
  //leg1->AddEntry(h_jt_closure_ref_nocorr_px[0], "Pre-corr Inclusive", "lepf");
  leg1->AddEntry(h_jt_closure_ref_ncs2_px[0], "Inclusive jets", "lepf");
  leg1->AddEntry(h_jt_closure_q_ncs2_px[0], "q jets", "lepf");
  leg1->AddEntry(h_jt_closure_g_ncs2_px[0], "g jets", "lepf");
   c_closure->cd(1);
   leg1->Draw("same");

   c_closure->cd(2);
   TLegend *leg2 = new TLegend(0.5,0.8,0.99,0.99);
   leg2->SetLineColor(kWhite);
   leg2->AddEntry(h_jt_closure_q_ncs2[0], "q jets", "lepf");
   leg2->AddEntry(h_jt_closure_g_ncs2[0], "g jets", "lepf");
   //leg2->Draw("same");
  /*
   TFile *f_pol1;
   f_pol1 = new TFile("/home/dhanush/Documents/JEC/local/new_corr_files/fitclosure_Apr17_pol1.root", "RECREATE");
   
   f_corr->cd();

   //f_corr->cd();

   f_flatcorr2_cent0->Write();
   f_flatcorr2_cent1->Write();
   f_flatcorr2_cent2->Write();
   f_flatcorr2_cent3->Write();
   f_corr->Close();
*/

////////////////////////////////////////////////////////

///////////////////  ncs cand dist /////////////////////

////////////////////////////////////////////////////////

  TCanvas *c_dist = new TCanvas("c_dist","",1200,300);
  c_dist->Divide(4,1);
  gStyle->SetOptStat(0);

  for(int i=0; i<4; i++){
    c_dist->cd(i+1);
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetRangeUser(0.,0.2);
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->SetRangeUser(0.,20.);
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetLabelSize(0.06);
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->SetLabelSize(0.07);
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetNdivisions(5);
    h_ncs_2_dist_q[3-i][3]->SetTitle(centVars[3-i].c_str());
    if(i==0) h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetTitle("fraction of inclusive jets");
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetTitleSize(0.09);
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->SetTitleOffset(0.8);
    h_ncs_2_dist_q[3-i][3]->GetYaxis()->CenterTitle();
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->SetTitle("nCS cand");
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->SetTitleSize(0.08);
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->SetTitleOffset(0.9);
    h_ncs_2_dist_q[3-i][3]->GetXaxis()->CenterTitle();
    h_ncs_2_dist_q[3-i][3]->Draw("same");
    h_ncs_2_dist_g[3-i][3]->Draw("same");
    //l1[3-i][j+5] = new TLatex(130,1.3,centVars[3-i].c_str());
    //l1[3-i][j+5]->SetTextSize(0.08);
    //l1[3-i][j+5]->Draw("same");
    if(i==1){
      TLegend *legend = new TLegend(0.5,0.75,0.9,0.9);
      legend ->SetLineColor(kWhite);
      legend->SetFillColor(0);
      legend ->AddEntry((TObject*)0, "130 < p_{T} < 140", "");
      legend ->Draw("same");
    }
    if(i==2){
      TLegend *legend = new TLegend(0.5,0.75,0.9,0.9);
      legend ->SetLineColor(kWhite);
      legend->SetFillColor(0);
      legend ->AddEntry((TObject*)0,"|#eta| < 1.6", "");
      legend ->Draw("same");
    }
    if(i==0){
      TLegend *leg0 = new TLegend(0.2,0.5,0.75,0.9);
      leg0->SetLineColor(kWhite);
      leg0->SetFillColor(0);
      leg0->AddEntry(h_ncs_2_dist_q[0][3],"q jets","lepf");
      leg0->AddEntry(h_ncs_2_dist_g[0][3],"g jets","lepf");
      leg0->Draw();
    }
  }


/////////////////////////////////////////////////////////

  TCanvas *c_closure_reco = new TCanvas("c_closure_reco","",1200,300);
  c_closure_reco->Divide(4,1);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c_closure_reco->cd(i+1);
    //h_jt_closure_reco_nocorr_px[3-i] = h_jt_closure_reco_nocorr[3-i]->ProfileX();
    h_jt_closure_reco_nocorr_px[3-i]->GetYaxis()->SetRangeUser(0.95,1.3);
    h_jt_closure_reco_nocorr_px[3-i]->GetXaxis()->SetRangeUser(90.,500.);
    h_jt_closure_reco_nocorr_px[3-i]->GetXaxis()->SetTitle("reco pT");
    h_jt_closure_reco_nocorr_px[3-i]->GetYaxis()->SetTitle("closure");
    h_jt_closure_reco_nocorr_px[3-i]->SetMarkerColor(kBlack);
    h_jt_closure_reco_nocorr_px[3-i]->SetMarkerStyle(29);
    h_jt_closure_reco_nocorr_px[3-i]->SetTitle(centVars[3-i].c_str());
    h_jt_closure_reco_nocorr_px[3-i]->Draw("e1 same");
    
    //tl1->Draw("same");
    //tl2->Draw("same");
    tl3->Draw("same");
    tl5->Draw("same");
  }


}