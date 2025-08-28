void checkrun_cf(int nmax=20) {
    TFile *f = new TFile("/eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST/1713717018_10.root","READ");
    TTree *t = (TTree*)f->Get("amstreea");
    Long64_t n = t->GetEntries();
    n = (nmax>0 && nmax<n) ? nmax : n;
    float cutoff[4][2], thetam, phim;
    UInt_t run;

    t->SetBranchAddress("mcutoffi", &cutoff);
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("thetam", &thetam);
    t->SetBranchAddress("phim", &phim);
    printf("%-10s %-15s %-10s\n", "Row", "run", "cutoff[1][1]");
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        printf("%-10lld %u %f %-10.6f %-10.6f\n", i, run, cutoff[1][1], thetam, phim);
    }
}