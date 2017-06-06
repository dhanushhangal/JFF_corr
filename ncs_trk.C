#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

using namespace std;

const int nCBins = 4;
const int nptBins = 55;

int mypbin, mycbin, myptbin, myrefptbin;

char saythis[500];

TString cent[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};

float CBins[5] = {0, 20, 60, 100, 200};

enum enum_dataset_types {e_Data2015,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};

int dataset_type_code = -999;

//arg 1 = which data set, arg 2 = output file number
void ncs_trk(bool doCrab=0, int jobID=0, int endfile = 500, int dataset_type_code = 2, int output_file_num = 1)
{

    //defining histos
    TH2F *h_cs_trk[nCBins];
    TH2F *h_pf_trk[nCBins];
    TH2F *h_pf_cs[nCBins];
    TH1D *h_cs_dist;
    TH1D *h_cs_id;
    TH1D *h_jt_pt[nCBins];

    for (int ibin=0;ibin<nCBins;ibin++){

      //cs_trk
      sprintf(saythis,"h_cs_trk_cent%d",ibin);
      h_cs_trk[ibin] = new TH2F(saythis, "",30,0,30,30,0,30);
      h_cs_trk[ibin]->Sumw2();

      //cs_trk
      sprintf(saythis,"h_pf_trk_cent%d",ibin);
      h_pf_trk[ibin] = new TH2F(saythis, "",30,0,30,30,0,30);
      h_pf_trk[ibin]->Sumw2();

      //cs_trk
      sprintf(saythis,"h_pf_cs_cent%d",ibin);
      h_pf_cs[ibin] = new TH2F(saythis, "",30,0,30,30,0,30);
      h_pf_cs[ibin]->Sumw2();

      //jtpt
      sprintf(saythis,"h_jt_pt_cent%d",ibin);
      h_jt_pt[ibin] = new TH1D(saythis, "",30,120.,150.);
      h_jt_pt[ibin]->Sumw2();      

    }

    //csdist
    h_cs_dist = new TH1D("h_cs_dist", "",40,0.,40.);
    h_cs_dist->Sumw2();

    //pfId
    h_cs_id = new TH1D("h_cs_id", "",10,0,10);
    h_cs_id->Sumw2();

	bool is_data = false;

	if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;

	// assert(!is_data); //for now I'm interested in MC

	//-----------------------------
	// Set JFF-dependent corrections
	//-----------------------------

	float reco_eta, reco_phi, reco_pt;
	
	bool do_PbPb=1, do_pp_tracking=0;

	int radius = 4;
	
	if(dataset_type_code== 1 || dataset_type_code > 10){do_PbPb = 0;   do_pp_tracking = 1;}

	cout<<"do_PbPb = "<<do_PbPb<<endl;
	cout<<"do_pp_tracking = "<<do_pp_tracking<<endl;


	//////////###### centrality Bins ###########///////////////

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree7=0;
	TTree *pftree;
	TTree *inp_tree_CS;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "ppdata_filelist.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPbdata_filelist_small.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "ppMC_pthat80_filelist.txt";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		in_file_name = "PbPb_MC_filelist_withCSfix_small.txt";
	}else{
		cerr<<"need to set up to run on that sample..."<<endl;
	}

	cout << "trying a filelist named "<< in_file_name << endl;

	TFile *output_file = new TFile("histos_data.root", "RECREATE");
	//TTree *mixing_tree = new TTree("mixing_tree", "");

	const int MAXPARTICLES = 100000;

	Int_t HBHENoiseFilter;
	Int_t eventSelection, pvFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_ps, HLT_Jet100_ps;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t hiBin = -999;
	Float_t evtPlane_HF2 = -999;
	Float_t pthat = -999;
	Float_t vz = -999;//, sumpt[15];
	Int_t nPFpart, nCSpart;

    std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;
	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS];
	Float_t t_calo_refpt[MAXJETS];

	Int_t t_calo_parton_flavor[MAXJETS];	
	Float_t t_trkPt[MAXPARTICLES], t_trkEta[MAXPARTICLES], t_trkPhi[MAXPARTICLES], t_trkDxy1[MAXPARTICLES], t_trkDxyError1[MAXPARTICLES], t_trkDz1[MAXPARTICLES], t_trkDzError1[MAXPARTICLES], t_trkPtError[MAXPARTICLES], t_pfEcal[MAXPARTICLES], t_pfHcal[MAXPARTICLES], t_trkChi2[MAXPARTICLES];
	Bool_t t_trkMVALoose[MAXPARTICLES], t_trkMVATight[MAXPARTICLES];
	UChar_t t_trkAlgo[MAXPARTICLES], t_trkNHit[MAXPARTICLES], t_trkNlayer[MAXPARTICLES], t_trkNdof[MAXPARTICLES]; //Run2
	//Float_t t_trkAlgo[MAXPARTICLES], t_trkNdof[MAXPARTICLES]; //Run1
	//Int_t t_trkNHit[MAXPARTICLES], t_trkNlayer[MAXPARTICLES]; //Run1
	Bool_t t_highPurity[MAXPARTICLES];

	Float_t t_hiEvtPlane[50];

	vector<float> *csCandEta=0, *csCandPhi=0, *csCandPt=0;
	vector<int> *csCandId=0;

	vector<float> *pfCandEta=0, *pfCandPhi=0, *pfCandPt=0;
        vector<int> *pfCandId=0;

	Int_t nTrk, ngen, mult, calo_nref, pf_nref;

    int total_n_jets = 0;

	while(instr>>filename && ifile<endfile){
		filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
		cout<<"File name is "<< filename <<endl;
		ifile++;

		//TFile *my_file = TFile::Open(filename.c_str());

		//if(!my_file){
			int pos = filename.find_first_of('s');
			string reducedfn = filename.substr(pos-1);
			string xrdPrefix = "root://cmsxrootd.fnal.gov//";
			cout << "local file not detected. Trying " << xrdPrefix+reducedfn << endl;
			TFile *my_file = TFile::Open((xrdPrefix+reducedfn).c_str());
			//TFile::Open((xrdPrefix+reducedfn).c_str());
		//}

        //TFile *my_file = TFile::Open("/data/forests/temp/HiForestAOD_80.root");

/*		
		if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	

		if(my_file->IsZombie()) { 
			std::cout << "Is zombie" << std::endl;
		}    
*/
		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
			pftree = (TTree*) my_file->Get(Form("akPu%dPFJetAnalyzer/t",radius));
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
			pftree = (TTree*) my_file->Get(Form("ak%dPFJetAnalyzer/t",radius));
		}

		inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
		if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree2);

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree3);

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);

		if(do_PbPb){
			inp_tree5 = (TTree*) my_file->Get("anaTrack/trackTree");
		}else{
			inp_tree5 = (TTree*) my_file->Get("ppTrack/trackTree");
		}
		if(!inp_tree5){ cout << "track Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree5);

		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(is_data && !inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		if(!is_data){ 
			inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
			if(!inp_tree7){ cout << "GenPart Tree not found!! Exiting..." << endl; exit(1); }
			else inp_tree->AddFriend(inp_tree7);
		}

		cout << "trees loaded" << endl;

		inp_tree->SetBranchAddress("vz",&vz);
		if(do_PbPb){
			inp_tree->SetBranchAddress("hiBin",&hiBin);
			inp_tree->SetBranchAddress("hiEvtPlanes",t_hiEvtPlane);
		}
	    inp_tree->SetBranchAddress("nTrk",&nTrk);
		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);

		inp_tree->SetBranchAddress("trkPt",t_trkPt);
		inp_tree->SetBranchAddress("trkEta",t_trkEta);
		inp_tree->SetBranchAddress("trkPhi",t_trkPhi);
		inp_tree->SetBranchAddress("trkNlayer",t_trkNlayer);
		inp_tree->SetBranchAddress("trkChi2",t_trkChi2);
		inp_tree->SetBranchAddress("trkNdof",t_trkNdof);

		inp_tree->SetBranchAddress("trkAlgo",t_trkAlgo);
		inp_tree->SetBranchAddress("highPurity",t_highPurity);
		
		if(!is_data) pftree->SetBranchAddress("pthat",&pthat);

		//Run2
		inp_tree->SetBranchAddress("nPFpart",&nPFpart);
		inp_tree->SetBranchAddress("pfPt",&pfCandPt);
                inp_tree->SetBranchAddress("pfEta",&pfCandEta);
                inp_tree->SetBranchAddress("pfPhi",&pfCandPhi);
                inp_tree->SetBranchAddress("pfId",&pfCandId);
                
		if(do_PbPb){
			inp_tree_CS->SetBranchAddress("nPFpart",&nCSpart);
			inp_tree_CS->SetBranchAddress("pfPt",&csCandPt);
			inp_tree_CS->SetBranchAddress("pfEta",&csCandEta);
			inp_tree_CS->SetBranchAddress("pfPhi",&csCandPhi);
			inp_tree_CS->SetBranchAddress("pfId",&csCandId);
		}


		inp_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilter);
		
		if(!do_PbPb){
			//inp_tree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter); //Run2
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		}

		if(is_data){
			if(!do_PbPb){ 
				inp_tree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_Jet100);
			}
			else{
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
			}
		}


		int n_evt = inp_tree->GetEntriesFast();

		//n_evt=2;
		cout << "Entries: "<< n_evt << endl;
		for(int evi = 0; evi < n_evt; evi++) {

            if(evi%20==0) cout<<evi<<endl;
			inp_tree->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);
			pftree->GetEntry(evi);

			if(is_data && !HLT_Jet80 && !HLT_Jet100) continue;

			if(!do_PbPb && (!pvFilter)) continue;

			if(do_PbPb) evtPlane_HF2 = t_hiEvtPlane[8]; 
			else evtPlane_HF2 = -999;

			if(do_PbPb && (!pvFilter || !HBHENoiseFilter || !eventSelection)) continue;
            
            if (hiBin == 0 ) {continue; }
                 
            for (int cbin = 0; cbin < nCBins; cbin++){ 
              if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

                mycbin = cbin; 
              }
            }

			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if( fabs(t_calo_jteta[j4i]) > 1.6 ) continue;

                if( t_calo_jtpt[j4i] < 120. || t_calo_jtpt[j4i] > 150.) continue;

                total_n_jets++; 

                reco_pt = t_calo_jtpt[j4i];
				reco_phi = t_calo_jtphi[j4i];
				reco_eta = t_calo_jteta[j4i];
                
                h_jt_pt[mycbin]->Fill(reco_pt);
				
				int ncs2_id1 = 0;
				if(do_PbPb){
					for(int icand=0; icand<nCSpart; icand++){
						if(abs(csCandEta->at(icand))>2.4) continue;
						double dr = sqrt(pow(csCandEta->at(icand) - reco_eta, 2) + pow(acos(cos(reco_phi - csCandPhi->at(icand))),2));
						if(dr < 0.4){
                          h_cs_id->Fill(csCandId->at(icand));
						  if(csCandPt->at(icand)>2 && csCandId->at(icand)==1){
						    ncs2_id1++;
						  } 
						}
					}
				}
				h_cs_dist->Fill(ncs2_id1);


                //npf cand for PbPb and P+H
				int npf2_id1 = 0;
				if(do_PbPb){
					for(int icand=0; icand<nPFpart; icand++){
						if(abs(pfCandEta->at(icand))>2.4) continue;
						double dr = sqrt(pow(pfCandEta->at(icand) - reco_eta, 2) + pow(acos(cos(reco_phi - pfCandPhi->at(icand))),2));
						if(dr < 0.4){
						  if(pfCandPt->at(icand)>2 && pfCandId->at(icand)==1){
						    npf2_id1++;
						  } 
						}
					}
				}				

                // reco tracks
                int ntracksInCone = 0; 

                for(int itrk=0; itrk < nTrk; itrk++){

                	if(!t_highPurity[itrk] || (float)t_trkChi2[itrk]/(float)t_trkNdof[itrk]/(float)t_trkNlayer[itrk] > 0.15 ) continue;
                    double dr = sqrt(pow(t_trkEta[itrk] - reco_eta,2) + pow(acos(cos(reco_phi - t_trkPhi[itrk])),2));
                    
                    if(dr < 0.4 && (t_trkPt[itrk])>2.){

                    	ntracksInCone++;
                    }
                }  

                h_cs_trk[mycbin]->Fill(ncs2_id1,ntracksInCone);
                h_pf_trk[mycbin]->Fill(npf2_id1,ntracksInCone);
                h_pf_cs[mycbin]->Fill(npf2_id1,ncs2_id1);

			} /// calo jet loop
		
		}  ///event loop

		my_file->Close();
	}  //file loop

    h_cs_dist->Scale(1./total_n_jets);
    h_cs_id->Scale(1./total_n_jets);

	cout<<"writing"<<endl;

	//mixing_tree->Write();
	output_file->cd();
	for(int ibin=0;ibin<nCBins;ibin++){
	  h_cs_trk[ibin]->Write();
      h_pf_trk[ibin]->Write();
      h_pf_cs[ibin]->Write();
      h_jt_pt[ibin]->Write();
	}

	h_cs_dist->Write();
	h_cs_id->Write();

	output_file->Close();

	cout<<"done"<<endl;

}
