void cpSaveTree() {
    TFile *f1 = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MCBeryllium10.root");
    TFile *f2 = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be10.root");
    
    TTree *t1 = (TTree*)f1->Get("saveTree");
    TTree *t2 = (TTree*)f2->Get("saveTree");
    
    Int_t cutStatus1;
    UInt_t cutStatus2;
    
    // 只启用需要的分支
    t1->SetBranchStatus("*", 0);
    t2->SetBranchStatus("*", 0);
    t1->SetBranchStatus("cutStatus", 1);
    t2->SetBranchStatus("cutStatus", 1);
    
    t1->SetBranchAddress("cutStatus", &cutStatus1);
    t2->SetBranchAddress("cutStatus", &cutStatus2);
    
    Long64_t nentries = t1->GetEntries();
    
    for(Long64_t i = 0; i < nentries; i++) {
        t1->GetEntry(i);
        t2->GetEntry(i);
        
        Int_t cuts1 = cutStatus1 & 0x7;
        Int_t cuts2 = cutStatus2 & 0x7;
        
        if(cuts1 != cuts2) {
            cout << "Entry " << i << ": CutStatus mismatch!" << endl;
            cout << "File1: cutStatus(first 3 bits)=" << std::bitset<3>(cuts1) << endl;
            cout << "File2: cutStatus(first 3 bits)=" << std::bitset<3>(cuts2) << endl;
        }
        
        if(i % 1000000 == 0) {
            cout << "Processed " << i << " entries..." << endl;
        }
    }
    
    f1->Close();
    f2->Close();
    delete f1;
    delete f2;
}