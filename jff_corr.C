#include <iostream>
#include "TFile.h"
#include "TRandom.h"
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
#include <TCut.h>
#include <vector>
#include "TCanvas.h"
#include "mixing_tree_new.h"
#include "nCScorr.h"

const int nCBins = 4;
const int nptBins = 55;
//const int nptBins = 10;
const int eta_nbins = 18;

using namespace std;

const bool derive_jff = false;
const bool derive_residual = false;
const bool is_pp = false;

int mypbin, mycbin, myptbin, myrefptbin;

char saythis[500];

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};
//TString pt[11] = {"0","1","2","3","4","5","6","7","8","9","10"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};
//int jt_nbins = 10;
//Double_t jt_bin_bounds[11] = {100., 110., 120., 130., 145., 160., 180., 220., 270., 350., 500.};

Double_t eta_bin_bounds[18] = {-1.5, -1.25,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.25, 1.5};

float CBins[5] = {0, 20, 60, 100, 200};

double calo_jtpt, calo_refpt, calo_jteta, calo_jtphi, closure_nocorr, closure_ncs2, closure_gen_pol1, closure_reco_pol1, hiBin, vz, pthat_weight, corr_factor2, corr_calo_jtpt_2, corr_calo_jtpt_pol1, corr_calo_jtpt_reco_pol1;

void jff_corr(){

  nCScorr *corrpt = new nCScorr(is_pp);

//// defining histos and profiles

  TH1D *h_reco[nCBins][nptBins];
  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_corr[nCBins];
  TH1D *h_gen_full[nCBins];
  TH1D *h_reco_ratio[nCBins]; 
  TH1F *h_closure_nocorr[nCBins][nptBins];
  TH1F *h_closure_ncs2[nCBins][nptBins];
  TH2F *h_ncs_pt[nCBins][nptBins];
  TH2F *h_ncs_pt_q[nCBins][nptBins];
  TH2F *h_ncs_pt_g[nCBins][nptBins];
  TH2F *h_ncs_2_nocorr[nCBins][nptBins];
  TH2F *h_ncs_2_nocorr_q[nCBins][nptBins];
  TH2F *h_ncs_2_nocorr_g[nCBins][nptBins];
  TH2F *h_ncs_2_corr[nCBins][nptBins];
  TH2F *h_ncs_2_corr_q[nCBins][nptBins];
  TH2F *h_ncs_2_corr_g[nCBins][nptBins];
  
  TH2F *h_jt_closure_ref_nocorr[nCBins];
  TH2F *h_jt_closure_ref_ncs2[nCBins];
  TH2F *h_jt_closure_ref_pol1[nCBins];
  TH2F *h_jt_closure_ref_pol1_gencorr[nCBins];
  TH2F *h_jt_closure_reco_nocorr[nCBins];
  TH2F *h_jt_closure_reco_ncs2[nCBins];
  TH2F *h_jt_closure_reco_pol1[nCBins];
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
  TProfile *h_jt_closure_ref_ncs2_px[nCBins];
  TProfile *h_jt_closure_ref_pol1_px[nCBins];
  TProfile *h_jt_closure_reco_px[nCBins];
  TProfile *h_jt_closure_reco_pol1_px[nCBins];

  TF1 *f_ncs_2[nCBins][nptBins];
  TF1 *f_closure[nCBins];
  TF1 *f_reco_ratio[nCBins];
  TF1 *f_flatcorr[nCBins];  
  TF1 *f2_par0[nCBins];
  TF1 *f2_par1[nCBins];

  TGraphErrors *gr_p1_param_2[nCBins];
  TGraphErrors *gr_p0_param_2[nCBins];

  Double_t x_mean[nptBins], x_mean_err[nptBins], par0_2[nptBins], par1_2[nptBins], par0_err_2[nptBins], par1_err_2[nptBins];

  for (int ibin=0;ibin<nCBins;ibin++){

    //closure histo vs refpt

    sprintf(saythis,"h_jt_closure_ref_nocorr_cent%d",ibin);
    h_jt_closure_ref_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_ref_ncs2_cent%d",ibin);
    h_jt_closure_ref_ncs2[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_ncs2[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_ref_pol1_cent%d",ibin);
    h_jt_closure_ref_pol1[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_pol1[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_ref_pol1_gencorr_cent%d",ibin);
    h_jt_closure_ref_pol1_gencorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_pol1_gencorr[ibin]->Sumw2();

    //closure histo vs recopt

    sprintf(saythis,"h_jt_closure_reco_nocorr_cent%d",ibin);
    h_jt_closure_reco_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_reco_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_reco_ncs2_cent%d",ibin);
    h_jt_closure_reco_ncs2[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_reco_ncs2[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_reco_pol1_cent%d",ibin);
    h_jt_closure_reco_pol1[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_reco_pol1[ibin]->Sumw2();

    //closure histo q vs refpt

    sprintf(saythis,"h_jt_closure_q_nocorr_cent%d",ibin);
    h_jt_closure_q_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);  
    h_jt_closure_q_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_q_ncs2_cent%d",ibin);
    h_jt_closure_q_ncs2[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);  
    h_jt_closure_q_ncs2[ibin]->Sumw2();

    //closure histo g vs refpt

    sprintf(saythis,"h_jt_closure_g_nocorr_cent%d",ibin);
    h_jt_closure_g_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_g_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_g_ncs2_cent%d",ibin);
    h_jt_closure_g_ncs2[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_g_ncs2[ibin]->Sumw2();

    //closure histo vs eta 
    sprintf(saythis,"h_eta_closure_nocorr_cent%d",ibin);
    h_eta_closure_nocorr[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);
    h_eta_closure_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_eta_closure_q_nocorr_cent%d",ibin);
    h_eta_closure_q_nocorr[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);  
    h_eta_closure_q_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_eta_closure_g_nocorr_cent%d",ibin);
    h_eta_closure_g_nocorr[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);
    h_eta_closure_g_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_eta_closure_ncs2_cent%d",ibin);
    h_eta_closure_ncs2[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);
    h_eta_closure_ncs2[ibin]->Sumw2();
    sprintf(saythis,"h_eta_closure_q_ncs2_cent%d",ibin);
    h_eta_closure_q_ncs2[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);  
    h_eta_closure_q_ncs2[ibin]->Sumw2();
    sprintf(saythis,"h_eta_closure_g_ncs2_cent%d",ibin);
    h_eta_closure_g_ncs2[ibin] = new TH2F(saythis, "", eta_nbins-1,eta_bin_bounds,200,0.,2.);
    h_eta_closure_g_ncs2[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_cent%d",ibin);
    h_reco_full[ibin] = new TH1D(saythis,"",500,0.,500.);
    h_reco_full[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_cent%d",ibin);
    h_reco_corr[ibin] = new TH1D(saythis,"",500,0.,500.);
    h_reco_corr[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_cent%d",ibin);
    h_gen_full[ibin] = new TH1D(saythis,"",500,0.,500.);
    h_gen_full[ibin]->Sumw2();

    sprintf(saythis,"h_reco_ratio_cent%d",ibin);
    h_reco_ratio[ibin] = new TH1D(saythis,"",55,50.,600.);
    h_reco_ratio[ibin]->Sumw2();

    for (int ibin2=0;ibin2<nptBins;ibin2++){

      sprintf(saythis,"h_reco_cent%d_pt%d",ibin,ibin2);
      h_reco[ibin][ibin2] = new TH1D(saythis,"",jt_bin_bounds[ibin2+1]-jt_bin_bounds[ibin2],jt_bin_bounds[ibin2],jt_bin_bounds[ibin2+1]);
      h_reco[ibin][ibin2]->Sumw2(); 

      sprintf(saythis,"h_closure_nocorr_cent%d_pt%d",ibin,ibin2);
      h_closure_nocorr[ibin][ibin2] = new TH1F(saythis,"",200,0.,2.);  
      h_closure_nocorr[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_closure_ncs2_cent%d_pt%d",ibin,ibin2);
      h_closure_ncs2[ibin][ibin2] = new TH1F(saythis,"",200,0.,2.);  
      h_closure_ncs2[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_nocorr_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_nocorr[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_nocorr[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_nocorr_q_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_nocorr_q[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_nocorr_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_nocorr_g_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_nocorr_g[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_nocorr_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_corr_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_corr[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_corr[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_corr_q_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_corr_q[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_corr_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs2_corr_g_cent%d_pt%d",ibin,ibin2);
      h_ncs_2_corr_g[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,200,0.,2.);  
      h_ncs_2_corr_g[ibin][ibin2]->Sumw2();
    
      sprintf(saythis,"h_ncs_pt_cent%d_pt%d",ibin,ibin2);
      h_ncs_pt[ibin][ibin2] = new TH2F(saythis,"",500,0.,500.,40,0.,40.);
      h_ncs_pt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs_pt_q_cent%d_pt%d",ibin,ibin2);
      h_ncs_pt_q[ibin][ibin2] = new TH2F(saythis,"",500,0.,500.,40,0.,40.);
      h_ncs_pt_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_ncs_pt_g_cent%d_pt%d",ibin,ibin2);
      h_ncs_pt_g[ibin][ibin2] = new TH2F(saythis,"",500,0.,500.,40,0.,40.);
      h_ncs_pt_g[ibin][ibin2]->Sumw2();
    }
  }

  ///// CUTS ///////
  const double etacut = 1.6;
  const double pTmincut = 50.;
  const double pTmaxcut = 600.;
  const double refpTmincut = 0.;
  //const double refpTmaxcut = 600.; 
  
  ///////////////// centrality reweighting ///////////////////////

  TFile *cen = TFile::Open("VertexHiNcollFits.root");

  TF1 *f_cen = (TF1*)cen->Get("xfit_hi")->Clone("f_cen");

  ///////////////// vz reweighting ///////////////////////

  TF1 *f_vz = (TF1*)cen->Get("xfit_vz")->Clone("f_vz"); 

  ////recopt correction file

  TH2D *h_jt_closure_reco_nocorr_corr_2D[4];
  TProfile *h_jt_closure_reco_nocorr_corr[4];
  
  TFile *f_reco;

  if(is_pp) f_reco = TFile::Open("/home/dhanush/Documents/JEC/local/ppreco_corr.root");
  else f_reco = TFile::Open("/home/dhanush/Documents/JEC/local/reco_corr.root");

  for (int ibin = 0; ibin < nCBins; ibin++){

    h_jt_closure_reco_nocorr_corr_2D[ibin] = (TH2D*)f_reco->Get((TString)("h_jt_closure_reco_nocorr_cent"+cent[ibin]))->Clone((TString)("h_jt_closure_reco_nocorr_cent"+cent[ibin]+"_pfx"));
    h_jt_closure_reco_nocorr_corr[ibin] = h_jt_closure_reco_nocorr_corr_2D[ibin]->ProfileX();

  }

  /////open main file //////

  TFile *my_file;

  if(is_pp) my_file = TFile::Open("unzippedSkim_5TeV_ppMC.root"); 
  else my_file = TFile::Open("unzippedSkim_PbPbMC_full.root");
  
  TTree *inp_tree = (TTree*)my_file->Get("unzipMixTree");
  mixing_tree_new *my_primary = new mixing_tree_new(inp_tree);
  std::cout << "Successfully retrieved tree from input file!" << std::endl;
  Long64_t n_jets = my_primary->fChain->GetEntriesFast();
  cout<<n_jets<<endl;

  //// Loop over all reco jets ////

  for (int jet = 0; jet < n_jets; jet++){ 
  //for (int jet = 0; jet < 200000; jet++){

    if (jet%1000000==0) cout<<jet<<endl;

    my_primary->fChain->GetEntry(jet);

    calo_jteta = my_primary->jteta;  
    if(fabs(calo_jteta) >= etacut) continue ;  

    calo_jtpt = my_primary->jtpt;
    if (calo_jtpt <= pTmincut || calo_jtpt >= pTmaxcut) continue;
            
    calo_refpt = my_primary->refpt;
    if (calo_refpt <= refpTmincut/* || calo_refpt >= refpTmaxcut*/) continue;

    int refparton_flavor = my_primary->refparton_flavor;

    int nCS_2 = my_primary->nCScandPt2;    
/*    
    //// vz and weight

    vz = my_primary->vz;
    double weight_vz = f_vz->Eval(vz);
*/
    //// centrality bin and weight 

    hiBin = my_primary->hiBin;
    if(is_pp) hiBin = 1;

    if (hiBin == 0 ) {continue; }
                 
    double weight_cen = f_cen->Eval(hiBin);

    if(is_pp) weight_cen = 1.;

    for (int cbin = 0; cbin < nCBins; cbin++){ 

      if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

        mycbin = cbin; 
      }
    }

    if(is_pp) mycbin = 0;

    ///// pthat weight

    pthat_weight = my_primary->weight;

    ///// pt bin

    for (int ptbin = 0; ptbin < nptBins; ptbin++){

      if (calo_jtpt > jt_bin_bounds[ptbin] && calo_jtpt <= jt_bin_bounds[ptbin+1]){
        
        myptbin = ptbin; 
      }
    }

    for (int ptbin = 0; ptbin < nptBins; ptbin++){

      if (calo_refpt > jt_bin_bounds[ptbin] && calo_refpt <= jt_bin_bounds[ptbin+1]){
        
        myrefptbin = ptbin; 
      }
    }

    //////apply ncs correction

    //// ncs pt>2
    corr_calo_jtpt_2 = corrpt->getCorrection(is_pp,nCS_2, hiBin, calo_jtpt, calo_jteta);

    /////// closure ///////// 

    closure_nocorr = calo_jtpt/calo_refpt;
    closure_ncs2 = corr_calo_jtpt_2/calo_refpt;

    /////filling histos

    h_ncs_pt[mycbin][myptbin]->Fill(calo_jtpt,nCS_2,pthat_weight*weight_cen);

    h_ncs_2_nocorr[mycbin][myptbin]->Fill(nCS_2,(closure_nocorr/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);    
    h_closure_nocorr[mycbin][myrefptbin]->Fill(closure_nocorr,pthat_weight*weight_cen);  

    h_ncs_2_corr[mycbin][myptbin]->Fill(nCS_2,(closure_ncs2/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);    
    h_closure_ncs2[mycbin][myrefptbin]->Fill(closure_ncs2,pthat_weight*weight_cen);

    h_reco[mycbin][myptbin]->Fill(calo_jtpt,pthat_weight*weight_cen);

    h_jt_closure_ref_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen);
    h_eta_closure_nocorr[mycbin]->Fill(calo_jteta,closure_nocorr,pthat_weight*weight_cen);
    h_jt_closure_ref_ncs2[mycbin]->Fill(calo_refpt,closure_ncs2,pthat_weight*weight_cen);
    h_eta_closure_ncs2[mycbin]->Fill(calo_jteta,closure_ncs2,pthat_weight*weight_cen);

    h_jt_closure_reco_nocorr[mycbin]->Fill(calo_jtpt,closure_nocorr,pthat_weight*weight_cen);
    h_jt_closure_reco_ncs2[mycbin]->Fill(calo_jtpt,closure_ncs2,pthat_weight*weight_cen);

    h_reco_full[mycbin]->Fill(calo_jtpt,pthat_weight*weight_cen);

    h_reco_corr[mycbin]->Fill(corr_calo_jtpt_2,pthat_weight*weight_cen);

    h_gen_full[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen);

    if(refparton_flavor == -999) continue;
               
    if (fabs(refparton_flavor) == 21){
      h_jt_closure_g_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen);
      h_eta_closure_g_nocorr[mycbin]->Fill(calo_jteta,closure_nocorr,pthat_weight*weight_cen);
      h_jt_closure_g_ncs2[mycbin]->Fill(calo_refpt,closure_ncs2,pthat_weight*weight_cen);
      h_eta_closure_g_ncs2[mycbin]->Fill(calo_jteta,closure_ncs2,pthat_weight*weight_cen);
      h_ncs_pt_g[mycbin][myptbin]->Fill(calo_jtpt,nCS_2,pthat_weight*weight_cen);
      h_ncs_2_nocorr_g[mycbin][myptbin]->Fill(nCS_2,(closure_nocorr/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);
      h_ncs_2_corr_g[mycbin][myptbin]->Fill(nCS_2,(closure_ncs2/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);
    }   
    else {
      h_jt_closure_q_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen);
      h_eta_closure_q_nocorr[mycbin]->Fill(calo_jteta,closure_nocorr,pthat_weight*weight_cen);
      h_jt_closure_q_ncs2[mycbin]->Fill(calo_refpt,closure_ncs2,pthat_weight*weight_cen);
      h_eta_closure_q_ncs2[mycbin]->Fill(calo_jteta,closure_ncs2,pthat_weight*weight_cen);
      h_ncs_pt_q[mycbin][myptbin]->Fill(calo_jtpt,nCS_2,pthat_weight*weight_cen);
      h_ncs_2_nocorr_q[mycbin][myptbin]->Fill(nCS_2,(closure_nocorr/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);
      h_ncs_2_corr_q[mycbin][myptbin]->Fill(nCS_2,(closure_ncs2/(h_jt_closure_reco_nocorr_corr[mycbin]->GetBinContent(myptbin+1))),pthat_weight*weight_cen);   
    }

  } //end of jet loop

  //jff correction parameters
  if(derive_jff){
  cout<<"deriving jff corrections"<<endl;

  for(int ibin=0;ibin<nCBins;ibin++){

    for(int ibin3=0;ibin3<nptBins;ibin3++){

      sprintf(saythis,"f_ncs_2_cent%d_pt%d",ibin,ibin3);
      f_ncs_2[ibin][ibin3] = new TF1(saythis, "1.+[1]*(x-[0])", 2., 30.);
      f_ncs_2[ibin][ibin3] ->SetParameter(0,10.);
      f_ncs_2[ibin][ibin3] ->SetParameter(1,-0.015);

      h_ncs_2_closure[ibin][ibin3] = h_ncs_2_nocorr[ibin][ibin3]->ProfileX();
      //h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q R M");
      
      // linear fits to nCS closures
      if(ibin3<5) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,18.); 
     
      else if(ibin3==7) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,18.);
      else if((ibin3==8) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,20.);
      else if(ibin3<9) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,17.);
      else if((ibin3==11) && (ibin==1)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.);
      else if((ibin3==11) && (ibin==2)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,25.);
      else if(ibin3<12) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);
      else if((ibin3==13) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",0.,21.);
      else if((ibin3==27) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,18.);
      else if((ibin3==28) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);
      else if((ibin3==31) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",0.,19.);
      else if((ibin3==32) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",1.,15.);
      else if(ibin3==33) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,24.);
      else if((ibin3==34) && (ibin==3)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,23.);
      else if((ibin3==34) && (ibin==0)) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,16.);
      else if(ibin3<35) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,28.);
      else if(ibin3==43 && ibin==1) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,22.);
      else if(ibin3==47 && ibin==2) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,24.);
      else if(ibin3==48 && ibin==0) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,12.);
      else if(ibin3==49 && ibin==1) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",3.,19.);
      else if(ibin3>44 && ibin==0) h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,15.);

      else h_ncs_2_closure[ibin][ibin3]->Fit(((TString)("f_ncs_2_cent"+cent[ibin]+"_pt"+pt[ibin3])),"Q M","",2.,22.);

      par0_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParameter(0);
      par1_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParameter(1);

      par0_err_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParError(0);
      par1_err_2[ibin3] = f_ncs_2[ibin][ibin3]->GetParError(1);

      // mean pt in the pt bin
      x_mean[ibin3] = h_reco[ibin][ibin3]->GetMean();
      x_mean_err[ibin3] = h_reco[ibin][ibin3]->GetMeanError();

    }
    
    f2_par0[ibin] = new TF1((TString)("f2_par0_cent"+cent[ibin]), "exp(([0]/(x-[2]))+[1])", 80., 550.);
    f2_par1[ibin] = new TF1((TString)("f2_par1_cent"+cent[ibin]),"([2]/x)*exp(([0]/((x-[3])^2))+([1]/x))",80.,550.);

    if(ibin==0){
      f2_par1[ibin]->SetParameter(0,-8.627e+06);
      f2_par1[ibin]->SetParameter(1,-199.9);
      f2_par1[ibin]->SetParameter(2,-15.14);
      f2_par1[ibin]->SetParameter(3,3236);

      f2_par0[ibin]->SetParameter(0,-339.6);
      f2_par0[ibin]->SetParameter(1,2.746);
      f2_par0[ibin]->SetParameter(2,-310.5);
    
      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->SetName((TString)("gr_par1_cent"+cent[ibin]));
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",80.,550.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->SetName((TString)("gr_par0_cent"+cent[ibin]));
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",80.,550.);  
    }

    if(ibin==1){
      f2_par1[ibin]->SetParameter(0,-331.2);
      f2_par1[ibin]->SetParameter(1,-122.9);
      f2_par1[ibin]->SetParameter(2,-4.509);
      f2_par1[ibin]->SetParameter(3,-0.1186);

      f2_par0[ibin]->SetParameter(0,-142.5);
      f2_par0[ibin]->SetParameter(1,2.673);
      f2_par0[ibin]->SetParameter(2,-113.5);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->SetName((TString)("gr_par1_cent"+cent[ibin]));
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",80.,550.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->SetName((TString)("gr_par0_cent"+cent[ibin]));
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",80.,550.);
      
    }
    
    if(ibin==2){
      f2_par1[ibin]->SetParameter(0,2.049e+05);
      f2_par1[ibin]->SetParameter(1,-859.7);
      f2_par1[ibin]->SetParameter(2,-13.44);
      f2_par1[ibin]->SetParameter(3,-73.72);

      f2_par0[ibin]->SetParameter(0,-112.8);
      f2_par0[ibin]->SetParameter(1,2.612);
      f2_par0[ibin]->SetParameter(2,-69.07);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->SetName((TString)("gr_par1_cent"+cent[ibin]));
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",80.,520.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->SetName((TString)("gr_par0_cent"+cent[ibin]));
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",80.,520.);
    }
    
    if(ibin==3){
      f2_par1[ibin]->SetParameter(0,9.643e+04);
      f2_par1[ibin]->SetParameter(1,-517.2);
      f2_par1[ibin]->SetParameter(2,-9.333);
      f2_par1[ibin]->SetParameter(3,-58.84);

      f2_par0[ibin]->SetParameter(0,-102.4);
      f2_par0[ibin]->SetParameter(1,2.539);
      f2_par0[ibin]->SetParameter(2,-70.84);

      gr_p1_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par1_2,x_mean_err,par1_err_2);
      gr_p1_param_2[ibin]->SetName((TString)("gr_par1_cent"+cent[ibin]));
      gr_p1_param_2[ibin]->Fit(f2_par1[ibin],"Q M","",80.,550.);

      gr_p0_param_2[ibin] = new TGraphErrors(nptBins,x_mean,par0_2,x_mean_err,par0_err_2);
      gr_p0_param_2[ibin]->SetName((TString)("gr_par0_cent"+cent[ibin]));
      gr_p0_param_2[ibin]->Fit(f2_par0[ibin],"Q M","",80.,550.);
    }
   }
  } 
  
  if(derive_residual){
    
    for(int ibin=0;ibin<nCBins;ibin++){
      sprintf(saythis,"f_closure_cent%d",ibin);
      f_closure[ibin] = new TF1(saythis, "pol1", 80., 500.);
      f_closure[ibin] ->SetParameter(0,0.995);
      f_closure[ibin] ->SetParameter(1,2e-05);
    
      h_jt_closure_ref_ncs2_px[ibin] = h_jt_closure_ref_ncs2[ibin]->ProfileX();
      h_jt_closure_ref_ncs2_px[ibin]->Fit(f_closure[ibin],"Q M","",80.,500.);
    
      //flat corrections
      sprintf(saythis,"f_flatcorr_cent%d",ibin);
      f_flatcorr[ibin] = new TF1(saythis, "pol0", 80., 600.);
      f_flatcorr[ibin] ->SetParameter(0,1.);

      h_jt_closure_ref_ncs2_px[ibin]->Fit(f_flatcorr[ibin],"Q M","",80.,600.);
    }

    //// Loop over all reco jets ////
    cout<<"deriving residual corrections"<<endl;

    for (int jet = 0; jet < n_jets; jet++){ 
    //for (int jet = 0; jet < 2000000; jet++){

      if (jet%1000000==0) cout<<jet<<endl;

      my_primary->fChain->GetEntry(jet);
 
      calo_jteta = my_primary->jteta;  
      if(fabs(calo_jteta) >= etacut) continue ;  

      calo_jtpt = my_primary->jtpt;
      if (calo_jtpt <= pTmincut || calo_jtpt >= pTmaxcut) continue;
              
      calo_refpt = my_primary->refpt;
      if (calo_refpt <= refpTmincut) continue;

      int nCS_2 = my_primary->nCScandPt2;
      
      //// centrality bin and weight 

      hiBin = my_primary->hiBin;
      if(is_pp) hiBin = 1;

      if (hiBin == 0 ) {continue; }
                 
      double weight_cen = f_cen->Eval(hiBin);

      if(is_pp) weight_cen = 1.;

      for (int cbin = 0; cbin < nCBins; cbin++){ 

        if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

          mycbin = cbin; 
        }
      }

      if(is_pp) mycbin = 0;

      ///// pthat weight

      pthat_weight = my_primary->weight;

      ///// pt bin

      for (int ptbin = 0; ptbin < nptBins; ptbin++){

        if (calo_jtpt > jt_bin_bounds[ptbin] && calo_jtpt <= jt_bin_bounds[ptbin+1]){
        
          myptbin = ptbin; 
        }
      }

      for (int ptbin = 0; ptbin < nptBins; ptbin++){

        if (calo_refpt > jt_bin_bounds[ptbin] && calo_refpt <= jt_bin_bounds[ptbin+1]){
        
          myrefptbin = ptbin; 
        }
      }

      ///after jff correction
      corr_calo_jtpt_2 = corrpt->getCorrection(is_pp,nCS_2, hiBin, calo_jtpt, calo_jteta);

      //////apply pol1 correction
      corr_calo_jtpt_pol1 = corr_calo_jtpt_2/(f_closure[mycbin]->Eval(calo_refpt));
    
      /////// closure ///////// 

      closure_gen_pol1 = corr_calo_jtpt_pol1/calo_refpt;

      /////filling histos

      h_jt_closure_ref_pol1_gencorr[mycbin]->Fill(calo_refpt,closure_gen_pol1,pthat_weight*weight_cen);
      h_jt_closure_reco_pol1[mycbin]->Fill(calo_jtpt,closure_gen_pol1,pthat_weight*weight_cen);
 
    } //end of jet loop

    for(int i=0; i<nCBins; i++){
      h_jt_closure_reco_px[i] = h_jt_closure_reco_ncs2[i]->ProfileX();
      h_jt_closure_reco_pol1_px[i] = h_jt_closure_reco_pol1[i]->ProfileX();
      
      sprintf(saythis,"f_reco_ratio_cent%d",i);
      f_reco_ratio[i] = new TF1(saythis, "pol1", 50., 600.);
      f_reco_ratio[i] ->SetParameter(0,0.99);
      f_reco_ratio[i] ->SetParameter(1,2e-05);
      
      for(int j=0; j<nptBins; j++){
        h_reco_ratio[i]->SetBinContent(j+1,(h_jt_closure_reco_px[i]->GetBinContent(j+1))/(h_jt_closure_reco_pol1_px[i]->GetBinContent(j+1)));
        h_reco_ratio[i]->SetBinError(j+1,sqrt((h_jt_closure_reco_px[i]->GetBinContent(j+1))/(h_jt_closure_reco_pol1_px[i]->GetBinContent(j+1)))*sqrt(pow((h_jt_closure_reco_px[i]->GetBinError(j+1))/(h_jt_closure_reco_px[i]->GetBinContent(j+1)),2)+pow((h_jt_closure_reco_pol1_px[i]->GetBinError(j+1))/(h_jt_closure_reco_pol1_px[i]->GetBinContent(j+1)),2)));
        h_reco_ratio[i]->Fit(f_reco_ratio[i],"Q M","",50.,600.);

      }
    }
/*
    cout<<"deriving flat corrections"<<endl;
    for (int jet = 0; jet < n_jets; jet++){ 
    //for (int jet = 0; jet < 2000000; jet++){

      if (jet%1000000==0) cout<<jet<<endl;

      my_primary->fChain->GetEntry(jet);

      calo_jteta = my_primary->jteta;  
      if(fabs(calo_jteta) >= etacut) continue ;  
  
      calo_jtpt = my_primary->jtpt;
      if (calo_jtpt <= pTmincut || calo_jtpt >= pTmaxcut) continue;
            
      calo_refpt = my_primary->refpt;
      if (calo_refpt <= refpTmincut) continue;
    
      int nCS_2 = my_primary->nCScandPt2;

      //// centrality bin and weight 

      hiBin = my_primary->hiBin;
      if(is_pp) hiBin = 1;

      if (hiBin == 0 ) {continue; }
                 
      double weight_cen = f_cen->Eval(hiBin);

      if(is_pp) weight_cen = 1.;

      for (int cbin = 0; cbin < nCBins; cbin++){ 

        if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

          mycbin = cbin; 
        }
      }

      if(is_pp) mycbin = 0;

      ///// pthat weight

      pthat_weight = my_primary->weight;

      ///// pt bin

      for (int ptbin = 0; ptbin < nptBins; ptbin++){

        if (calo_jtpt > jt_bin_bounds[ptbin] && calo_jtpt <= jt_bin_bounds[ptbin+1]){
        
          myptbin = ptbin; 
        }
      }

      for (int ptbin = 0; ptbin < nptBins; ptbin++){

        if (calo_refpt > jt_bin_bounds[ptbin] && calo_refpt <= jt_bin_bounds[ptbin+1]){
        
          myrefptbin = ptbin; 
        }
      }

      ///after jff correction
      corr_calo_jtpt_2 = corrpt->getCorrection(is_pp,nCS_2, hiBin, calo_jtpt, calo_jteta);

      //////apply reco pol1 correction

      //corr_calo_jtpt_reco_pol1 = corr_calo_jtpt_2/(f_reco_ratio[mycbin]->Eval(corr_calo_jtpt_2));
    
      /////// closure ///////// 

      closure_reco_pol1 = corr_calo_jtpt_reco_pol1/calo_refpt;

      /////filling histos

      h_jt_closure_ref_pol1[mycbin]->Fill(calo_refpt,closure_reco_pol1,pthat_weight*weight_cen);
 
    } //end of jet loop
*/ 
  }

/// recopt closure correction file ///
/*
  TFile *reco_corr;

  reco_corr = new TFile("/home/dhanush/Documents/JEC/local/ppreco_corr.root", "RECREATE");

  for(int ibin=0;ibin<nCBins;ibin++){

    h_jt_closure_reco_nocorr[ibin]->Write();  

  }
  reco_corr.Close();
*/

/// writing histos ///
  TFile *closure_histos;

  if(is_pp) closure_histos = new TFile("/home/dhanush/Documents/JEC/local/ppclosure_histos_Jun1_header.root", "RECREATE");
  else closure_histos = new TFile("/home/dhanush/Documents/JEC/local/closure_histos_Jun1_header.root", "RECREATE");
  closure_histos->cd();

  for(int ibin=0;ibin<nCBins;ibin++){
    h_jt_closure_ref_nocorr[ibin]->Write();
    h_eta_closure_nocorr[ibin]->Write();
    h_jt_closure_ref_ncs2[ibin]->Write();
    h_eta_closure_ncs2[ibin]->Write();

    h_jt_closure_reco_nocorr[ibin]->Write();
    h_jt_closure_reco_ncs2[ibin]->Write();
  
    h_jt_closure_q_nocorr[ibin]->Write();
    h_eta_closure_q_nocorr[ibin]->Write();
    h_jt_closure_q_ncs2[ibin]->Write();
    h_eta_closure_q_ncs2[ibin]->Write();

    h_jt_closure_g_nocorr[ibin]->Write();
    h_eta_closure_g_nocorr[ibin]->Write();
    h_jt_closure_g_ncs2[ibin]->Write();
    h_eta_closure_g_ncs2[ibin]->Write();

    if(derive_residual) h_jt_closure_ref_ncs2_px[ibin]->Write();
    h_jt_closure_reco_pol1[ibin]->Write();
    h_jt_closure_ref_pol1_gencorr[ibin]->Write();
    h_jt_closure_ref_pol1[ibin]->Write();

    h_reco_full[ibin]->Write();
    h_reco_corr[ibin]->Write();
    h_gen_full[ibin]->Write();
    h_reco_ratio[ibin]->Write();
    
    if(derive_jff){
      gr_p1_param_2[ibin]->Write();
      gr_p0_param_2[ibin]->Write();
    }

    for(int ibin3=0;ibin3<nptBins;ibin3++){
      h_ncs_2_nocorr[ibin][ibin3]->Write();
      h_ncs_2_nocorr_q[ibin][ibin3]->Write();
      h_ncs_2_nocorr_g[ibin][ibin3]->Write();
      h_ncs_2_corr[ibin][ibin3]->Write();
      h_ncs_2_corr_q[ibin][ibin3]->Write();
      h_ncs_2_corr_g[ibin][ibin3]->Write();
      h_reco[ibin][ibin3]->Write();
      h_ncs_pt[ibin][ibin3]->Write();
      h_ncs_pt_q[ibin][ibin3]->Write();
      h_ncs_pt_g[ibin][ibin3]->Write();
      h_closure_nocorr[ibin][ibin3]->Write();
      h_closure_ncs2[ibin][ibin3]->Write();
      //f_ncs_2[ibin][ibin3]->Write();
      if(derive_jff){
        h_ncs_2_closure[ibin][ibin3]->Write();
      }
    }
  }

  closure_histos->Close();

if(derive_jff && derive_residual){
  /// writing correction files ///
  TFile *corr_file;

  if(is_pp) corr_file = new TFile("/home/dhanush/Documents/JEC/local/corr_files_May29/ppcorr_file_May29.root", "RECREATE");
  else corr_file = new TFile("/home/dhanush/Documents/JEC/local/corr_files_May29/corr_file_May29.root", "RECREATE");
  corr_file->cd();

  for(int ibin=0;ibin<nCBins;ibin++){
    f2_par1[ibin]->Write();
    f2_par0[ibin]->Write();

    f_reco_ratio[ibin]->Write();
    f_flatcorr[ibin]->Write();
  }
  corr_file->Close();
}

}

