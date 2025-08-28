//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 18 04:54:41 2024 by ROOT version 6.32.06
// from TTree saveTree/amstreea
// found on file: root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Lit.root
//////////////////////////////////////////////////////////

#ifndef analysTree_h
#define analysTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class analysTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   Float_t         tof_pos[4][3];
   Float_t         tof_chist;
   Float_t         tof_chisc;
   Float_t         tof_ql[4];
   Float_t         rich_pb;
   Int_t           rich_hit;
   Int_t           rich_usedm;
   Float_t         rich_q[2];
   Float_t         rich_npe[3];
   Bool_t          rich_good;
   Bool_t          rich_clean;
   Bool_t          rich_NaF;
   Int_t           rich_pmt;
   Float_t         rich_pos[3];
   Float_t         rich_theta;
   Float_t         rich_phi;
   UInt_t          cutStatus;
   Double_t        InnerRig;
   Double_t        richBeta;
   Double_t        tofBeta;
   Double_t        cutOffRig;
   Double_t        modiRichX;
   Double_t        modiRichY;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_tof_pos;   //!
   TBranch        *b_tof_chist;   //!
   TBranch        *b_tof_chisc;   //!
   TBranch        *b_tof_ql;   //!
   TBranch        *b_rich_pb;   //!
   TBranch        *b_rich_hit;   //!
   TBranch        *b_rich_usedm;   //!
   TBranch        *b_rich_q;   //!
   TBranch        *b_rich_npe;   //!
   TBranch        *b_rich_good;   //!
   TBranch        *b_rich_clean;   //!
   TBranch        *b_rich_NaF;   //!
   TBranch        *b_rich_pmt;   //!
   TBranch        *b_rich_pos;   //!
   TBranch        *b_rich_theta;   //!
   TBranch        *b_rich_phi;   //!
   TBranch        *b_cutStatus;   //!
   TBranch        *b_InnerRig;   //!
   TBranch        *b_richBeta;   //!
   TBranch        *b_tofBeta;   //!
   TBranch        *b_cutOffRig;   //!
   TBranch        *b_modiRichX;   //!
   TBranch        *b_modiRichY;   //!

   analysTree(TTree *tree=0);
   virtual ~analysTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysTree_cxx
analysTree::analysTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Lit.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Lit.root");
      }
      f->GetObject("saveTree",tree);

   }
   Init(tree);
}

analysTree::~analysTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysTree::LoadTree(Long64_t entry)
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

void analysTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("tof_pos", tof_pos, &b_tof_pos);
   fChain->SetBranchAddress("tof_chist", &tof_chist, &b_tof_chist);
   fChain->SetBranchAddress("tof_chisc", &tof_chisc, &b_tof_chisc);
   fChain->SetBranchAddress("tof_ql", tof_ql, &b_tof_ql);
   fChain->SetBranchAddress("rich_pb", &rich_pb, &b_rich_pb);
   fChain->SetBranchAddress("rich_hit", &rich_hit, &b_rich_hit);
   fChain->SetBranchAddress("rich_usedm", &rich_usedm, &b_rich_usedm);
   fChain->SetBranchAddress("rich_q", rich_q, &b_rich_q);
   fChain->SetBranchAddress("rich_npe", rich_npe, &b_rich_npe);
   fChain->SetBranchAddress("rich_good", &rich_good, &b_rich_good);
   fChain->SetBranchAddress("rich_clean", &rich_clean, &b_rich_clean);
   fChain->SetBranchAddress("rich_NaF", &rich_NaF, &b_rich_NaF);
   fChain->SetBranchAddress("rich_pmt", &rich_pmt, &b_rich_pmt);
   fChain->SetBranchAddress("rich_pos", rich_pos, &b_rich_pos);
   fChain->SetBranchAddress("rich_theta", &rich_theta, &b_rich_theta);
   fChain->SetBranchAddress("rich_phi", &rich_phi, &b_rich_phi);
   fChain->SetBranchAddress("cutStatus", &cutStatus, &b_cutStatus);
   fChain->SetBranchAddress("InnerRig", &InnerRig, &b_InnerRig);
   fChain->SetBranchAddress("richBeta", &richBeta, &b_richBeta);
   fChain->SetBranchAddress("tofBeta", &tofBeta, &b_tofBeta);
   fChain->SetBranchAddress("cutOffRig", &cutOffRig, &b_cutOffRig);
   fChain->SetBranchAddress("modiRichX", &modiRichX, &b_modiRichX);
   fChain->SetBranchAddress("modiRichY", &modiRichY, &b_modiRichY);
   Notify();
}

bool analysTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void analysTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysTree_cxx
