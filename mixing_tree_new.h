//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  2 22:09:28 2017 by ROOT version 5.34/14
// from TTree unzipMixTree/
// found on file: unzippedSkim_PbPbMC_full.root
//////////////////////////////////////////////////////////

#ifndef mixing_tree_new_h
#define mixing_tree_new_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class mixing_tree_new {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        jtpt;
   Double_t        corrpt;
   Double_t        refpt;
   Double_t        jteta;
   Double_t        jtphi;
   Int_t           nCScandPt2;
   Int_t           nCScandPt1;
   Int_t           nPFcand;
   Int_t           hiBin;
   Float_t         rxPlane;
   Int_t           refparton_flavor;
   Double_t        weight;

   // List of branches
   TBranch        *b_jtpt;   //!
   TBranch        *b_corrpt;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_nCScandPt2;   //!
   TBranch        *b_nCScandPt1;   //!
   TBranch        *b_nPFcand;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_rxPlane;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_weight;   //!

   mixing_tree_new(TTree *tree=0);
   virtual ~mixing_tree_new();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef mixing_tree_new_cxx
mixing_tree_new::mixing_tree_new(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("unzippedSkim_PbPbMC_full.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("unzippedSkim_PbPbMC_full.root");
      }
      f->GetObject("unzipMixTree",tree);

   }
   Init(tree);
}

mixing_tree_new::~mixing_tree_new()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mixing_tree_new::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mixing_tree_new::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mixing_tree_new::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("corrpt", &corrpt, &b_corrpt);
   fChain->SetBranchAddress("refpt", &refpt, &b_refpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("nCScandPt2", &nCScandPt2, &b_nCScandPt2);
   fChain->SetBranchAddress("nCScandPt1", &nCScandPt1, &b_nCScandPt1);
   fChain->SetBranchAddress("nPFcand", &nPFcand, &b_nPFcand);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("rxPlane", &rxPlane, &b_rxPlane);
   fChain->SetBranchAddress("refparton_flavor", &refparton_flavor, &b_refparton_flavor);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t mixing_tree_new::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mixing_tree_new::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mixing_tree_new::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mixing_tree_new_cxx
