void cp_bt() {
    auto f = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root");
    auto t = (TTree*)f->Get("saveTree");
    auto c = new TCanvas();
    c->SetGrid();
    
    vector<int> colors = {kRed, kBlue, kGreen+2};
    string path = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/";
    
    // Create all histograms
    for(int i = 1; i <= 3; i++) {
        t->Draw(Form("InnerRig/cutoffpi[1]>>h%d(100,0.1,2.1)", i), Form("btstat_new==%d", i));
        t->Draw(Form("InnerRig/cutOffRig[1]>>h%d_rig(1000,0.1,2.1)", i), Form("btstat_new==%d", i));
        for(auto h : {Form("h%d",i), Form("h%d_rig",i)}) {
            ((TH1D*)gDirectory->Get(h))->SetLineColor(colors[i-1]);
            ((TH1D*)gDirectory->Get(h))->SetLineWidth(3);
        }
    }

    // Three plots together
    vector<pair<string,string>> plots = {{"cutoffpi","hist"}, {"cutOffRig_hist","hist"}, {"cutOffRig_P","P"}};
    for(const auto& [name, opt] : plots) {
        for(int i = 1; i <= 3; i++) {
            string hname = (name.compare(0,8,"cutoffpi")==0) ? Form("h%d",i) : Form("h%d_rig",i);
            string drawOpt = (i==1) ? opt : opt + "same";
            ((TH1D*)gDirectory->Get(hname.c_str()))->Draw(drawOpt.c_str());
        }
        c->SaveAs((path + name + ".png").c_str());
    }

    // Two btstat plots
    for(int i = 1; i <= 2; i++) {
        ((TH1D*)gDirectory->Get(Form("h%d",i)))->Draw("hist");
        ((TH1D*)gDirectory->Get(Form("h%d_rig",i)))->Draw("histsame");
        c->SaveAs((path + Form("btstat%d.png",i)).c_str());
    }
}