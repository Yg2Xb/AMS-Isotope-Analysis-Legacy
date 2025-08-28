void cpe() {
    // 打开两个文件
    TFile *f1 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10/1451905820_10.root");
    TFile *f2 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10/1451911372_5.root");
    
    // 获取树
    TTree *t1 = (TTree*)f1->Get("saveTree");
    TTree *t2 = (TTree*)f2->Get("saveTree");
    
    // 定义变量
    unsigned int cutStatus1, cutStatus2;
    double AGLBeta1, AGLBeta2;
    unsigned int run1, event1, run2, event2;
    
    // 设置分支地址
    t1->SetBranchAddress("cutStatus", &cutStatus1);
    t1->SetBranchAddress("AGLBeta", &AGLBeta1);
    t1->SetBranchAddress("run", &run1);
    t1->SetBranchAddress("event", &event1);
    
    t2->SetBranchAddress("cutStatus", &cutStatus2);
    t2->SetBranchAddress("AGLBeta", &AGLBeta2);
    t2->SetBranchAddress("run", &run2);
    t2->SetBranchAddress("event", &event2);
    
    // 创建run条件的TCut
    TCut runCut = "run>=1451911372 && run<=1451916924";
    
    // 创建map来存储第一个文件的数据
    std::map<std::pair<int,int>, std::pair<unsigned int, double>> file1Data;
    
    // 处理第一个文件
    Long64_t nentries1 = t1->GetEntries();
    for(Long64_t i = 0; i < nentries1; i++) {
        t1->GetEntry(i);
        if(run1 >= 1451911372 && run1 <= 1451916924) {
            file1Data[{run1,event1}] = {cutStatus1, AGLBeta1};
        }
    }
    
    // 比较两个文件
    Long64_t nentries2 = t2->GetEntries();
    int mismatchCount = 0;
    
    for(Long64_t i = 0; i < nentries2; i++) {
        t2->GetEntry(i);
        if(run2 >= 1451911372 && run2 <= 1451916924) {
            auto it = file1Data.find({run2,event2});
            if(it != file1Data.end()) {
                // 检查cutStatus的前三位
                bool cutStatusMatch = ((it->second.first & 0x7) == (cutStatus2 & 0x7));
                bool AGLBetaMatch = (it->second.second == AGLBeta2);
                
                if((!cutStatusMatch || !AGLBetaMatch) && AGLBeta2 > 0) {
                    cout << "Mismatch found for run " << run2 << " event " << event2 << endl;
                    cout << "File1 cutStatus (first 3 bits): " << (it->second.first & 0x7) 
                         << " AGLBeta: " << it->second.second << endl;
                    cout << "File2 cutStatus (first 3 bits): " << (cutStatus2 & 0x7) 
                         << " AGLBeta: " << AGLBeta2 << endl;
                        cout<<it->second.second-AGLBeta2<<endl; 
                    cout << "-------------------" << endl;
                    mismatchCount++;
                }
            }
        }
    }
    
    cout << "Total events checked in run range: " << file1Data.size() << endl;
    cout << "Number of mismatches found: " << mismatchCount << endl;
    
    // 关闭文件
    f1->Close();
    f2->Close();
    delete f1;
    delete f2;
}