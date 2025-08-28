void drawcpL2() {
    gStyle->SetOptStat(0);
    int z = 5;
    int ibin = 3;

    // 文件路径
    const char* fin1 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p2.root";
    const char* fin2 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p5.root";
    const char* fin3 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p8.root";
    const char* fout = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/pngresult/cpL2Temp_Boron_bin1.png";

    // 打开文件
    TFile *f1 = TFile::Open(fin1);
    TFile *f2 = TFile::Open(fin2);
    TFile *f3 = TFile::Open(fin3);
    if (!f1 || !f2) { printf("file open failed\n"); return; }

    // 取hist
    TH2* h2_L1 = (TH2*)f1->Get("L1Temp_Boron_TOF_Bkg");
    TH2* h2_1 = (TH2*)f1->Get("L2Temp_Boron_TOF_Bkg");
    TH2* h2_2 = (TH2*)f2->Get("L2Temp_Boron_TOF_Bkg");
    TH2* h2_3 = (TH2*)f3->Get("L2Temp_Boron_TOF_Bkg");
    //if (!h2_1 || !h2_2 || !h2_3 || !h2_L1) { printf("hist get failed\n"); return; }

    // ProjectionX，第1个y-bin
    TH1* h1_L1 = h2_L1->ProjectionX("h1_L1", ibin, ibin+1);
    TH1* h1_1 = h2_1->ProjectionX("h1_1", ibin, ibin+1);
    TH1* h1_2 = h2_2->ProjectionX("h1_2", ibin, ibin+1);
    TH1* h1_3 = h2_3->ProjectionX("h1_3", ibin, ibin+1);
    h1_L1->Rebin(2);
    h1_1->Rebin(2);
    h1_2->Rebin(2);
    h1_3->Rebin(2);
    cout<<h1_L1->GetMaximum()<<endl;
    cout<<h1_1->GetMaximum()<<endl;
    cout<<h1_2->GetMaximum()<<endl;
    cout<<h1_3->GetMaximum()<<endl;

    double int_L1 = h1_L1->GetMaximum();//Integral(1, 100*(z-2)+50);
    if (int_L1 > 0) h1_L1->Scale(1.0 / int_L1);
    double int1 = h1_1->GetMaximum();//Integral(1, 100*(z-2)+50);
    if (int1 > 0) h1_1->Scale(1.0 / int1);
    h1_2->Scale(1.0 / h1_2->GetMaximum());
    h1_3->Scale(1.0 / h1_3->GetMaximum());
    /*
    double int2 = h1_2->Integral(1, 100*(z-2)+50);
    if (int2 > 0) h1_2->Scale(1.0 / int2);
    double int3 = h1_3->Integral(1, 100*(z-2)+50);
    if (int3 > 0) h1_3->Scale(1.0 / int3);
    */

    // 设置样式
    h1_1->SetLineColor(kRed);     h1_1->SetLineWidth(2);
    h1_2->SetLineColor(kBlue);    h1_2->SetLineWidth(2);
    h1_3->SetLineColor(kGreen+2); h1_3->SetLineWidth(2);
    h1_L1->SetLineColor(kBlack);  h1_L1->SetLineWidth(2);

    // 画图
    TCanvas* c1 = new TCanvas("c1", "cpL2Temp", 800, 600);
    c1->SetLogy();
    h1_1->SetTitle("Boron TOF 0.61-0.86 GeV/n;Tracker Single Layer Q;Normalized Counts");
    h1_1->GetXaxis()->SetRangeUser(z-2,z+2);
    //h1_2->GetXaxis()->SetRangeUser(z-2,z+2);
    h1_1->Draw("hist");
    h1_2->Draw("hist same");
    h1_3->Draw("hist same");
    h1_L1->Draw("hist same");

    // legend
    auto leg = new TLegend(0.4, 0.2, 0.7, 0.5); // 右上但不贴边
    leg->SetFillStyle(0); // 透明
    leg->SetBorderSize(0); // 无边框
    leg->AddEntry(h1_1, "0.2*Qcut L2Temp", "l");
    leg->AddEntry(h1_2, "0.5*Qcut L2Temp", "l");
    leg->AddEntry(h1_3, "0.8*Qcut L2Temp", "l");
    leg->AddEntry(h1_L1, "0.2*Qcut L1Temp", "l");
    leg->Draw();

    // 保存
    c1->SaveAs(fout);

    // 关闭文件
    f1->Close(); f2->Close(); 
}