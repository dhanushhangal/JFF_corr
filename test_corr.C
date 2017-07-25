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
#include "mixing_tree_cymbal.h"
#include "nCScorr.h"

const int nCBins = 4;
const int nptBins = 55;
const int npt_histoBins = 15;
const int eta_nbins = 18;
const int phi_nbins = 20;

using namespace std;

const bool is_cymbal = true;
const bool ispp = false;
const bool isdata = false;

TString dataset_type_file_names[4] = {"PbPb_Data_skim.txt","pp_Data_skim.txt","PbPb_MC_skim_new.txt","pp_MC_skim.txt"};

int mypbin, mycbin, myptbin, myrefptbin;

char saythis[500];

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};
Double_t pt_bounds[15] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 170., 220., 280., 370., 500.};

Double_t eta_bin_bounds[18] = {-1.5, -1.25,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.25, 1.5};
Double_t phi_bin_bounds[20] = {-1.50796, -1.256635, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531, 1.256635, 1.50796};

float CBins[5] = {0, 20, 60, 100, 200};

double vz, calo_jtpt, calo_corrpt, calo_refpt, calo_jteta, calo_jtphi, closure_nocorr, closure_corr, pt_size, hiBin, pthat_weight;

//Auxillary functions defined below
void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);

//cymbal tune centrality reweighting
double fcent(double centrality, TF1* fcent1){ 
  return (centrality < 194) ? fcent1->Eval(centrality) : 1;
}

void test_corr(){

  nCScorr *corrpt = new nCScorr(ispp,is_cymbal);

//// defining histos and profiles

  TH1D *h_reco[nCBins][nptBins];

  TH1D *h_hibin = new TH1D("h_hibin","",200,0.,200.);
  TH1D *h_hibin_nrw = new TH1D("h_hibin_nrw","",200,0.,200.);
  TH1D *h_vz[nCBins];
  TH1D *h_vz_nrw[nCBins];

  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_full_q[nCBins];
  TH1D *h_reco_full_g[nCBins];
  TH1D *h_reco_corr[nCBins];
  TH1D *h_reco_corr_q[nCBins];
  TH1D *h_reco_corr_g[nCBins];

  TH1D *h_eta_full[nCBins];
  TH1D *h_eta_full_nrw[nCBins];
  TH1D *h_eta_full_q[nCBins];
  TH1D *h_eta_full_g[nCBins];
  
  TH1D *h_phi_full[nCBins];
  TH1D *h_phi_full_nrw[nCBins];
  TH1D *h_phi_full_q[nCBins];
  TH1D *h_phi_full_g[nCBins];

  TH1D *h_gen_full[nCBins];
  TH2F *h_ncs_pt[nCBins];
  TH2F *h_npf_pt[nCBins];  

  TH2F *h_jt_closure_ref_nocorr[nCBins];
  TH2F *h_jt_closure_ref_corr[nCBins];
  TH2F *h_jt_closure_reco_nocorr[nCBins];
  TH2F *h_jt_closure_reco_corr[nCBins];
  TH2F *h_jt_closure_q_nocorr[nCBins];
  TH2F *h_jt_closure_q_corr[nCBins];
  TH2F *h_jt_closure_g_nocorr[nCBins];
  TH2F *h_jt_closure_g_corr[nCBins];

  for (int ibin=0;ibin<nCBins;ibin++){

    //closure histo vs refpt

    sprintf(saythis,"h_jt_closure_ref_nocorr_cent%d",ibin);
    h_jt_closure_ref_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_ref_corr_cent%d",ibin);
    h_jt_closure_ref_corr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_ref_corr[ibin]->Sumw2();

    //closure histo vs recopt

    sprintf(saythis,"h_jt_closure_reco_nocorr_cent%d",ibin);
    h_jt_closure_reco_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_reco_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_reco_corr_cent%d",ibin);
    h_jt_closure_reco_corr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_reco_corr[ibin]->Sumw2();

    //closure histo q vs refpt

    sprintf(saythis,"h_jt_closure_q_nocorr_cent%d",ibin);
    h_jt_closure_q_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);  
    h_jt_closure_q_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_q_corr_cent%d",ibin);
    h_jt_closure_q_corr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);  
    h_jt_closure_q_corr[ibin]->Sumw2();

    //closure histo g vs refpt

    sprintf(saythis,"h_jt_closure_g_nocorr_cent%d",ibin);
    h_jt_closure_g_nocorr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_g_nocorr[ibin]->Sumw2();
    sprintf(saythis,"h_jt_closure_g_corr_cent%d",ibin);
    h_jt_closure_g_corr[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,200,0.,2.);
    h_jt_closure_g_corr[ibin]->Sumw2();
/*
    sprintf(saythis,"h_gen_full_cent%d",ibin);
    h_gen_full[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_gen_full[ibin]->Sumw2();
*/
    sprintf(saythis,"h_gen_full_cent%d",ibin);
    h_gen_full[ibin] = new TH1D(saythis,"",npt_histoBins-1,pt_bounds);
    h_gen_full[ibin]->Sumw2();

/*
    sprintf(saythis,"h_reco_full_cent%d",ibin);
    h_reco_full[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_full[ibin]->Sumw2();
*/
    sprintf(saythis,"h_vz_cent%d",ibin);
    h_vz[ibin] = new TH1D(saythis, "", 30,-15.,15.);
    h_vz[ibin]->Sumw2();

    sprintf(saythis,"h_vz_nrw_cent%d",ibin);
    h_vz_nrw[ibin] = new TH1D(saythis, "", 30,-15.,15.);
    h_vz_nrw[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_cent%d",ibin);
    h_reco_full[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_full[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_q_cent%d",ibin);
    h_reco_full_q[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_g_cent%d",ibin);
    h_reco_full_g[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_cent%d",ibin);
    h_eta_full[ibin] = new TH1D(saythis, "", eta_nbins-1,eta_bin_bounds);
    h_eta_full[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_nrw_cent%d",ibin);
    h_eta_full_nrw[ibin] = new TH1D(saythis, "", eta_nbins-1,eta_bin_bounds);
    h_eta_full_nrw[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_q_cent%d",ibin);
    h_eta_full_q[ibin] = new TH1D(saythis,"",eta_nbins-1,eta_bin_bounds);
    h_eta_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_g_cent%d",ibin);
    h_eta_full_g[ibin] = new TH1D(saythis,"",eta_nbins-1,eta_bin_bounds);
    h_eta_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_phi_full_cent%d",ibin);
    h_phi_full[ibin] = new TH1D(saythis, "", phi_nbins-1,phi_bin_bounds);
    h_phi_full[ibin]->Sumw2();

    sprintf(saythis,"h_phi_full_nrw_cent%d",ibin);
    h_phi_full_nrw[ibin] = new TH1D(saythis, "", phi_nbins-1,phi_bin_bounds);
    h_phi_full_nrw[ibin]->Sumw2();

    sprintf(saythis,"h_phi_full_q_cent%d",ibin);
    h_phi_full_q[ibin] = new TH1D(saythis,"",phi_nbins-1,phi_bin_bounds);
    h_phi_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_phi_full_g_cent%d",ibin);
    h_phi_full_g[ibin] = new TH1D(saythis,"",phi_nbins-1,phi_bin_bounds);
    h_phi_full_g[ibin]->Sumw2();            

/*
    sprintf(saythis,"h_reco_corr_cent%d",ibin);
    h_reco_corr[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_corr[ibin]->Sumw2();
*/
    sprintf(saythis,"h_reco_corr_cent%d",ibin);
    h_reco_corr[ibin] = new TH1D(saythis, "", npt_histoBins-1,pt_bounds);
    h_reco_corr[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_q_cent%d",ibin);
    h_reco_corr_q[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_corr_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_g_cent%d",ibin);
    h_reco_corr_g[ibin] = new TH1D(saythis,"",450,50.,500.);
    h_reco_corr_g[ibin]->Sumw2();

    //for (int ibin2=0;ibin2<nptBins;ibin2++){
      sprintf(saythis,"h_ncs_pt_cent%d",ibin);
      h_ncs_pt[ibin] = new TH2F(saythis,"",500,0.,500.,40,0.,40.);
      h_ncs_pt[ibin]->Sumw2();

      sprintf(saythis,"h_npf_pt_cent%d",ibin);
      h_npf_pt[ibin] = new TH2F(saythis,"",500,0.,500.,40,0.,40.);
      h_npf_pt[ibin]->Sumw2();
    //}

    for (int ibin2=0;ibin2<nptBins;ibin2++){
      sprintf(saythis,"h_reco_cent%d_pt%d",ibin,ibin2);
      h_reco[ibin][ibin2] = new TH1D(saythis,"",jt_bin_bounds[ibin2+1]-jt_bin_bounds[ibin2],jt_bin_bounds[ibin2],jt_bin_bounds[ibin2+1]);
      h_reco[ibin][ibin2]->Sumw2();  
    }
  }

  int total_n_jets = 0;

  ///// CUTS ///////
  const double etacut = 1.6;
  const double pTmincut = 50.;
  const double pTmaxcut = 600.;
  const double refpTmincut = 0.;
  
  ///////////////// centrality reweighting ///////////////////////

  TFile *cen = TFile::Open("VertexHiNcollFits.root");

  TF1 *f_cen = (TF1*)cen->Get("xfit_hi")->Clone("f_cen"); 

  TF1* f_cent= new TF1("f_cent","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);  
  f_cent->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);

  ///////////////// vz reweighting ///////////////////////

  TF1 *f_vz = (TF1*)cen->Get("xfit_vz")->Clone("f_vz");

  TF1 *fWeight = new TF1("fWeight","gaus(0)/(gaus(3))",-30.,30.);
  fWeight->SetParameters(0.08,0.44,5.12,0.08,3.25,5.23);

  /////open files //////
  
  std::vector<TString> file_names;   file_names.clear();

  ReadFileList(file_names, dataset_type_file_names[2], true);

  cout<<"got file"<<endl;

  for(int fi = 0; fi < (int) file_names.size(); fi++) {
  //for(int fi = 0; fi < 10; fi++) {
    
    TFile *my_file = TFile::Open(file_names.at(fi));
    //if(fi%5==0) std::cout << "Current file: " << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
    if(my_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }

/*
  TFile *my_file;

  if(ispp) my_file = TFile::Open("unzippedSkim_5TeV_ppMC.root"); 
  else my_file = TFile::Open("unzippedSkim_PbPbMC_full.root");
*/

  TTree *inp_tree = (TTree*)my_file->Get("unzipMixTree");
  mixing_tree_cymbal *my_primary = new mixing_tree_cymbal(inp_tree);
  std::cout << "Successfully retrieved tree from input file!" << std::endl;
  Long64_t n_jets = my_primary->fChain->GetEntriesFast();
  total_n_jets += n_jets;
  //cout<<total_n_jets<<endl;

  //// Loop over all reco jets ////

  for (int jet = 0; jet < n_jets; jet++){ 
  //for (int jet = 0; jet < 200000; jet++){
    
    if (jet%1000000==0) cout<<jet<<endl;

    my_primary->fChain->GetEntry(jet);

    calo_jteta = my_primary->jteta;  
    if(fabs(calo_jteta) >= etacut) continue ;  

    calo_jtpt = my_primary->jtpt;
    if (calo_jtpt <= pTmincut || calo_jtpt >= pTmaxcut) continue;
    
    calo_jtphi = my_primary->jtphi;
    //calo_corrpt = my_primary->corrpt;
 
    calo_refpt = my_primary->refpt;
    if (calo_refpt <= refpTmincut/* || calo_refpt >= refpTmaxcut*/) continue;

    int refparton_flavor = my_primary->refparton_flavor;

    int nCS_2 = my_primary->nCScandPt2_id145;
    int nPF =0;
    //int nCS_2 = my_primary->nCScand;
    //int nCS_2 = my_primary->nCScandPt2;  
    //cout<<nCS_2<<endl;  
    
    //// centrality bin and weight 

    if(ispp) hiBin = 1;
    else hiBin = my_primary->hiBin;

    if (hiBin == 0 ) {continue; }
                 
    //double weight_cen = f_cen->Eval(hiBin);
    double weight_cen = fcent(hiBin,f_cent);

    if(ispp || isdata) weight_cen = 1.;

    //// vz and weight

    vz = my_primary->vz;
    if (fabs(vz) > 15.) continue;
    //double weight_vz = f_vz->Eval(vz);
    double weight_vz = fWeight->Eval(vz);

    if(isdata) weight_vz = 1.;
    if(ispp) weight_vz = 1.;  

    for (int cbin = 0; cbin < nCBins; cbin++){ 

      if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

        mycbin = cbin; 
      }
    }

    if(ispp) mycbin = 0;

    ///// pthat weight

    pthat_weight = my_primary->weight;

    if(isdata) pthat_weight = 1.;

    ///// pt bin

    for (int ptbin = 0; ptbin < nptBins; ptbin++){

      if (calo_jtpt > jt_bin_bounds[ptbin] && calo_jtpt <= jt_bin_bounds[ptbin+1]){
        
        myptbin = ptbin; 
      }
    }

    if(!isdata){
      for (int ptbin = 0; ptbin < nptBins; ptbin++){

        if (calo_refpt > jt_bin_bounds[ptbin] && calo_refpt <= jt_bin_bounds[ptbin+1]){
        
          myrefptbin = ptbin; 
        }
      }
    }

    calo_corrpt = corrpt->getCorrection(ispp, nCS_2, hiBin, calo_jtpt, calo_jteta);

    if(!isdata){
    
    /////// closure ///////// 

      closure_nocorr = calo_jtpt/calo_refpt;
      closure_corr = calo_corrpt/calo_refpt;
    }

    /////filling histos
    h_reco[mycbin][myptbin]->Fill(calo_jtpt,pthat_weight*weight_cen*weight_vz);

    h_ncs_pt[mycbin]->Fill(calo_jtpt,nCS_2,pthat_weight*weight_cen*weight_vz);
    h_npf_pt[mycbin]->Fill(calo_jtpt,nPF,pthat_weight*weight_cen*weight_vz);

    if(!isdata){
      h_jt_closure_ref_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen*weight_vz);
      h_jt_closure_ref_corr[mycbin]->Fill(calo_refpt,closure_corr,pthat_weight*weight_cen*weight_vz);

      h_jt_closure_reco_nocorr[mycbin]->Fill(calo_jtpt,closure_nocorr,pthat_weight*weight_cen*weight_vz);
      h_jt_closure_reco_corr[mycbin]->Fill(calo_jtpt,closure_corr,pthat_weight*weight_cen*weight_vz);
       
      h_gen_full[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
    
      if(refparton_flavor == -999) continue;
               
      if (fabs(refparton_flavor) == 21){
        h_jt_closure_g_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen*weight_vz);
        h_jt_closure_g_corr[mycbin]->Fill(calo_refpt,closure_corr,pthat_weight*weight_cen*weight_vz);
        h_reco_full_g[mycbin]->Fill(calo_jtpt,pthat_weight*weight_cen*weight_vz);
        h_reco_corr_g[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
        h_eta_full_g[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        h_phi_full_g[mycbin]->Fill(calo_jtphi,pthat_weight*weight_cen*weight_vz);
      }   
      else {
        h_jt_closure_q_nocorr[mycbin]->Fill(calo_refpt,closure_nocorr,pthat_weight*weight_cen*weight_vz);
        h_jt_closure_q_corr[mycbin]->Fill(calo_refpt,closure_corr,pthat_weight*weight_cen*weight_vz);
        h_reco_full_q[mycbin]->Fill(calo_jtpt,pthat_weight*weight_cen*weight_vz);
        h_reco_corr_q[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
        h_eta_full_q[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        h_phi_full_q[mycbin]->Fill(calo_jtphi,pthat_weight*weight_cen*weight_vz);
      }
    }
    
    h_vz[mycbin]->Fill(vz,pthat_weight*weight_cen*weight_vz);
    h_vz_nrw[mycbin]->Fill(vz,pthat_weight);
    h_hibin->Fill(hiBin,pthat_weight*weight_cen*weight_vz);
    h_hibin_nrw->Fill(hiBin,pthat_weight);
    h_reco_full[mycbin]->Fill(calo_jtpt,pthat_weight*weight_cen*weight_vz);
    h_reco_corr[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
    h_eta_full[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
    h_phi_full[mycbin]->Fill(calo_jtphi,pthat_weight*weight_cen*weight_vz);
    h_eta_full_nrw[mycbin]->Fill(calo_jteta,pthat_weight);
    h_phi_full_nrw[mycbin]->Fill(calo_jtphi,pthat_weight);

    }//end of jet loop
  }//end of file loop

  TFile *closure_histos;

  if(isdata){
    if(ispp) closure_histos = new TFile("/home/dhanush/Documents/JFF_corrections/pptest_histos_data_Jun26.root", "RECREATE");
    else closure_histos = new TFile("/home/dhanush/Documents/JFF_corrections/test_histos_data_Jun26.root", "RECREATE");
  }
  else{
    if(ispp) closure_histos = new TFile("/home/dhanush/Documents/JFF_corrections/pptest_histos_MC_Jun30.root", "RECREATE");
    else closure_histos = new TFile("/home/dhanush/Documents/JFF_corrections/test_histos_MC_Jul25.root", "RECREATE");
  }  

  closure_histos->cd();

  for(int ibin=0;ibin<nCBins;ibin++){
 
  if(!isdata){
    h_jt_closure_ref_nocorr[ibin]->Write();
    h_jt_closure_ref_corr[ibin]->Write();

    h_jt_closure_reco_nocorr[ibin]->Write();
    h_jt_closure_reco_corr[ibin]->Write();
  
    h_jt_closure_q_nocorr[ibin]->Write();
    h_jt_closure_q_corr[ibin]->Write();

    h_jt_closure_g_nocorr[ibin]->Write();
    h_jt_closure_g_corr[ibin]->Write();
    h_gen_full[ibin]->Write();
  }

    h_reco_full[ibin]->Write();
    h_reco_corr[ibin]->Write();
    h_reco_full_q[ibin]->Write();
    h_reco_corr_q[ibin]->Write();
    h_reco_full_g[ibin]->Write();
    h_reco_corr_g[ibin]->Write();

    h_eta_full[ibin]->Write();
    h_eta_full_nrw[ibin]->Write();
    h_eta_full_q[ibin]->Write();
    h_eta_full_g[ibin]->Write();
    h_phi_full[ibin]->Write();
    h_phi_full_nrw[ibin]->Write();
    h_phi_full_q[ibin]->Write();
    h_phi_full_g[ibin]->Write();

    h_vz[ibin]->Write();
    h_vz_nrw[ibin]->Write();
    h_hibin->Write();
    h_hibin_nrw->Write();

    for (int ibin2=0;ibin2<nptBins;ibin2++){
      h_reco[ibin][ibin2]->Write();
    }
    h_ncs_pt[ibin]->Write();
    h_npf_pt[ibin]->Write();
    
  }

  closure_histos->Close();

}

void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
{
  ifstream file_stream(file_of_names);
  std::string line;
  my_file_names.clear();
  if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
  if( file_stream.is_open() ) {
    if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
    int line_num = 0;
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug ) std::cout << line_num << ": " << line << std::endl;
      TString tstring_line(line);
      if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
      line_num++;
    }
  } else {
    std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
    assert(0);
  }
}
