void draw_yyh() {
    // 打开输入文件
    TFile *f = TFile::Open("/afs/cern.ch/user/y/yyou/public/BeBackground/Result.root");
    if (!f) return;

    // 创建PDF文件
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/Carbon_yyh_result.pdf["); // 打开PDF文件

    // 设置画布属性
    c->SetLogx(1);
    c->SetLogy(1);
    //gPad->SetGrid();

    // 检索器和同位素名称数组
    const char* dets[] = {"TOF", "NaF", "AGL"};
    const char* bes[] = {"Be7", "Be9", "Be10"};
    const char* types[] = {"NumEk", "genEk", "recEk", "genEk_weighted", "recEk_weighted"};

    // 循环遍历所有直方图
    for (const char* det : dets) {
        for (const char* be : bes) {
                TString histName = Form("Carbon_%s_%s", be, det);
                TH1F* h = (TH1F*)f->Get(histName);
                if (!h) continue;
                h->GetXaxis()->SetRangeUser(0.3,15);

                c->cd();
                //h->SetTitle(Form("%s;Energy (GeV/n);Events", histName.Data()));
                h->SetLineColor(kBlue);
                h->SetLineWidth(2);
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.8);
                h->Draw("hist");
                
                c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/Carbon_yyh_result.pdf");
        }
    }


    // 关闭PDF文件
    c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/Carbon_yyh_result.pdf]");

    delete c;
    f->Close();
}