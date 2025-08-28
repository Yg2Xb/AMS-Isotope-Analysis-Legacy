void checkBT() {
    TFile *f = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7_ntrk.root");
    TTree *tree = (TTree*)f->Get("saveTree");

    // 打印总事例数
    Long64_t totalEntries = tree->GetEntries();
    cout << "Total entries: " << totalEntries << endl;

    // 方法1：使用GetEntries统计
    cout << "\nMethod 1 - Using GetEntries:" << endl;
    for(int i = 0; i <= 3; i++) {
        Long64_t count = tree->GetEntries(Form("btstat==%d", i));
        cout << "btstat=" << i << ": " << count << endl;
    }

    // 方法2：使用Scan查看前100个事例的分布
    cout << "\nMethod 2 - First 100 entries:" << endl;
    tree->Scan("btstat", "", "length=100");

    // 方法3：使用直方图统计
    cout << "\nMethod 3 - Using histogram:" << endl;
    TH1F *h = new TH1F("h", "btstat distribution", 4, -0.5, 3.5);
    tree->Draw("btstat>>h", "", "goff");
    for(int i = 0; i <= 3; i++) {
        cout << "btstat=" << i << ": " << h->GetBinContent(i+1) << endl;
    }

    // 方法4：手动遍历
    cout << "\nMethod 4 - Manual counting:" << endl;
    Int_t btstat;
    tree->SetBranchAddress("btstat", &btstat);
    int counts[4] = {0};
    
    // 只遍历前10000个事例作为样本
    Long64_t nEntries = min(totalEntries, (Long64_t)1000000);
    for(Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if(btstat >= 0 && btstat <= 3) {
            counts[btstat]++;
        } else {
            cout << "Found invalid btstat value: " << btstat << " at entry " << i << endl;
        }
    }
    
    for(int i = 0; i <= 3; i++) {
        cout << "btstat=" << i << ": " << counts[i] << " (in first " << nEntries << " entries)" << endl;
    }

    // 检查是否有任何其他的selection cuts
    cout << "\nCurrent tree selection (if any):" << endl;

    delete h;
    f->Close();
}