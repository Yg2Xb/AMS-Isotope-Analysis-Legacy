#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>
#include <iostream>
#include <TNetXNGFile.h>

void DrawPure() {
    // 设置ROOT样式
    gStyle->SetOptStat(0);
    
    // 打开ROOT文件
    TFile* file = TFile::Open("root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/ISS_FragmentBer_MassFitAnalysis_L1Inner_510cutoff_yzxQfit.root", "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    
    // 获取直方图
    TH1F* h_realfrac_sel_Be7_ratio = (TH1F*)file->Get("realfrac_sel_Be7_ratio");
    TH1F* h_realfrac_sel_Be9_ratio = (TH1F*)file->Get("realfrac_sel_Be9_ratio");
    TH1F* h_realfrac_sel_Be10_ratio = (TH1F*)file->Get("realfrac_sel_Be10_ratio");
    
    TH1F* h_fragment_Be_raw = (TH1F*)file->Get("h_fragment_Be_raw");
    TH1F* h_contamination_Be = (TH1F*)file->Get("h_contamination_Be");
    
    if (!h_realfrac_sel_Be7_ratio || !h_realfrac_sel_Be9_ratio || !h_realfrac_sel_Be10_ratio ||
        !h_fragment_Be_raw || !h_contamination_Be) {
        std::cerr << "Error: missing histograms in file" << std::endl;
        file->Close();
        return;
    }
    
    // 创建总Be纯度比例直方图
    TH1F* h_pure_be_ratio = (TH1F*)h_fragment_Be_raw->Clone("h_pure_be_ratio");
    h_pure_be_ratio->Add(h_contamination_Be, -1);
    
    // 使用Divide函数带B选项计算二项分布误差
    h_pure_be_ratio->Divide(h_pure_be_ratio, h_fragment_Be_raw, 1, 1, "B");
    
    // 绘制Be7比例
    TCanvas* c1 = new TCanvas("c1", "Be7 Ratio", 800, 600);
    h_realfrac_sel_Be7_ratio->SetTitle("^{7}Be after/before subtraction");
    h_realfrac_sel_Be7_ratio->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    h_realfrac_sel_Be7_ratio->GetYaxis()->SetTitle("ratio of ^{7}Be from L1 B");
    h_realfrac_sel_Be7_ratio->SetLineColor(kBlue);
    h_realfrac_sel_Be7_ratio->SetMarkerColor(kBlue);
    h_realfrac_sel_Be7_ratio->SetMarkerStyle(20);
    h_realfrac_sel_Be7_ratio->Draw("pz");
    c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be7_after_before_subtraction.png");
    
    // 绘制Be9比例
    TCanvas* c2 = new TCanvas("c2", "Be9 Ratio", 800, 600);
    h_realfrac_sel_Be9_ratio->SetTitle("^{9}Be after/before subtraction");
    h_realfrac_sel_Be9_ratio->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    h_realfrac_sel_Be9_ratio->GetYaxis()->SetTitle("ratio of ^{9}Be from L1 B");
    h_realfrac_sel_Be9_ratio->SetLineColor(kBlue);
    h_realfrac_sel_Be9_ratio->SetMarkerColor(kBlue);
    h_realfrac_sel_Be9_ratio->SetMarkerStyle(20);
    h_realfrac_sel_Be9_ratio->Draw("pz");
    c2->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be9_after_before_subtraction.png");
    
    // 绘制Be10比例
    TCanvas* c3 = new TCanvas("c3", "Be10 Ratio", 800, 600);
    h_realfrac_sel_Be10_ratio->SetTitle("^{10}Be after/before subtraction");
    h_realfrac_sel_Be10_ratio->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    h_realfrac_sel_Be10_ratio->GetYaxis()->SetRangeUser(0.5, 1.0);
    h_realfrac_sel_Be10_ratio->GetYaxis()->SetTitle("ratio of ^{10}Be from L1 B");
    h_realfrac_sel_Be10_ratio->GetYaxis()->SetNdivisions(508);
    h_realfrac_sel_Be10_ratio->SetLineColor(kRed);
    h_realfrac_sel_Be10_ratio->SetMarkerColor(kRed);
    h_realfrac_sel_Be10_ratio->SetMarkerStyle(20);
    h_realfrac_sel_Be10_ratio->SetMarkerSize(1.2);
    h_realfrac_sel_Be10_ratio->Draw("pz");
    c3->SetGridy();
    c3->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be10_after_before_subtraction.png");
    
    // 绘制总Be比例
    TCanvas* c4 = new TCanvas("c4", "Total Be Ratio", 800, 600);
    h_pure_be_ratio->SetTitle("Be after/before subtraction");
    h_pure_be_ratio->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    h_pure_be_ratio->GetYaxis()->SetTitle("ratio of Be from L1 B");
    h_pure_be_ratio->GetYaxis()->SetRangeUser(0.2, 0.8);
    h_pure_be_ratio->GetYaxis()->SetNdivisions(508);
    h_pure_be_ratio->SetLineColor(kBlack);
    h_pure_be_ratio->SetMarkerColor(kBlack);
    h_pure_be_ratio->SetMarkerStyle(20);
    h_pure_be_ratio->SetMarkerSize(1.2);
    h_pure_be_ratio->Draw("pz");
    c4->SetGridy();
    c4->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be_after_before_subtraction.png");
    
    // 关闭文件
    file->Close();
    
    std::cout << "Ratio plots saved as PNG files" << std::endl;
}