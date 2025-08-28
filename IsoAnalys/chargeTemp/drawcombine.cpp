#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <string>

void drawcombine() {
    // 设置ROOT样式
    gStyle->SetOptStat(0);
    
    // 打开输入文件
    TFile* egeFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_1010.root", "READ");
    TFile* egeLinearFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGELinear_1010.root", "READ");
    
    if (!egeFile || egeFile->IsZombie()) {
        std::cerr << "Error opening EGE file" << std::endl;
        return;
    }
    
    if (!egeLinearFile || egeLinearFile->IsZombie()) {
        std::cerr << "Error opening EGE Linear file" << std::endl;
        egeFile->Close();
        return;
    }
    
    // 获取输入直方图
    TH1D* h_be_ege = nullptr;
    TH1D* h_b_ege = nullptr;
    TH1D* h_c_ege = nullptr;
    
    TH1D* h_be_egelinear = nullptr;
    TH1D* h_b_egelinear = nullptr;
    TH1D* h_c_egelinear = nullptr;
    
    egeFile->GetObject("h_be_fraction", h_be_ege);
    egeFile->GetObject("h_b_fraction", h_b_ege);
    egeFile->GetObject("h_c_fraction", h_c_ege);
    
    egeLinearFile->GetObject("h_be_fraction", h_be_egelinear);
    egeLinearFile->GetObject("h_b_fraction", h_b_egelinear);
    egeLinearFile->GetObject("h_c_fraction", h_c_egelinear);
    
    if (!h_be_ege || !h_b_ege || !h_c_ege) {
        std::cerr << "Error: missing histograms in EGE file" << std::endl;
        egeFile->Close();
        egeLinearFile->Close();
        return;
    }
    
    if (!h_be_egelinear || !h_b_egelinear || !h_c_egelinear) {
        std::cerr << "Error: missing histograms in EGE Linear file" << std::endl;
        egeFile->Close();
        egeLinearFile->Close();
        return;
    }
    
    // 创建输出文件
    TFile* outputFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_combined.root", "RECREATE");
    
    // 创建合并的直方图
    TH1D* h_be_combined = (TH1D*)h_be_ege->Clone("h_be_fraction");
    TH1D* h_b_combined = (TH1D*)h_b_ege->Clone("h_b_fraction");
    TH1D* h_c_combined = (TH1D*)h_c_ege->Clone("h_c_fraction");
    
    h_be_combined->SetTitle("Be Fraction vs Energy (Combined)");
    h_b_combined->SetTitle("B Fraction vs Energy (Combined)");
    h_c_combined->SetTitle("C Fraction vs Energy (Combined)");
    
    // 合并模式: ege, egelinear, ege, ege, ege, egelinear, egelinear, 后面七个全是ege
    // 对于14个bin，这个模式如下：
    bool useEgeLinear[14] = {
        false,  // bin 1: ege
        true,   // bin 2: egelinear
        false,  // bin 3: ege
        false,  // bin 4: ege
        false,  // bin 5: ege
        true,   // bin 6: egelinear
        true,   // bin 7: egelinear
        false,  // bin 8: ege
        false,  // bin 9: ege
        false,  // bin 10: ege
        false,  // bin 11: ege
        false,  // bin 12: ege
        false,  // bin 13: ege
        false   // bin 14: ege
    };
    
    // 合并直方图
    int numBins = h_be_ege->GetNbinsX();
    if (numBins != 14) {
        std::cout << "Warning: Expected 14 bins, but found " << numBins << " bins" << std::endl;
    }
    
    for (int i = 1; i <= numBins; i++) {
        if (i <= 14 && useEgeLinear[i-1]) {
            // 使用EGE Linear的结果
            h_be_combined->SetBinContent(i, h_be_egelinear->GetBinContent(i));
            h_be_combined->SetBinError(i, h_be_egelinear->GetBinError(i));
            
            h_b_combined->SetBinContent(i, h_b_egelinear->GetBinContent(i));
            h_b_combined->SetBinError(i, h_b_egelinear->GetBinError(i));
            
            h_c_combined->SetBinContent(i, h_c_egelinear->GetBinContent(i));
            h_c_combined->SetBinError(i, h_c_egelinear->GetBinError(i));
        }
        // 其他bin使用EGE的结果（默认已经复制）
    }
    
    // 保存合并后的直方图
    outputFile->cd();
    h_be_combined->Write();
    h_b_combined->Write();
    h_c_combined->Write();
    
    // 同时复制其他必要的直方图或对象
    TH1D* h_signal_events = nullptr;
    egeFile->GetObject("h_signal_events", h_signal_events);
    if (h_signal_events) {
        h_signal_events->Clone("h_signal_events")->Write();
    }
    
    // 创建比较图像
    TCanvas* canvas_be = new TCanvas("canvas_be", "Be Fraction Comparison", 900, 600);
    h_be_ege->SetLineColor(kBlue);
    h_be_ege->SetMarkerColor(kBlue);
    h_be_ege->SetMarkerStyle(20);
    h_be_ege->GetYaxis()->SetRangeUser(0., 0.01);
    h_be_ege->Draw("P");
    
    h_be_egelinear->SetLineColor(kRed);
    h_be_egelinear->SetMarkerColor(kRed);
    h_be_egelinear->SetMarkerStyle(21);
    h_be_egelinear->Draw("P SAME");
    
    h_be_combined->SetLineColor(kGreen+2);
    h_be_combined->SetMarkerColor(kGreen+2);
    h_be_combined->SetMarkerStyle(22);
    h_be_combined->Draw("P SAME");
    
    TLegend* legend = new TLegend(0.65, 0.70, 0.88, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(h_be_ege, "EGE", "p");
    legend->AddEntry(h_be_egelinear, "EGE Linear", "p");
    legend->AddEntry(h_be_combined, "Combined", "p");
    legend->Draw();
    
    canvas_be->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/Comparison_Be_Fraction.png");
    
    // 创建Boron的比较图
    TCanvas* canvas_b = new TCanvas("canvas_b", "B Fraction Comparison", 900, 600);
    h_b_ege->SetLineColor(kBlue);
    h_b_ege->SetMarkerColor(kBlue);
    h_b_ege->SetMarkerStyle(20);
    h_b_ege->GetYaxis()->SetRangeUser(0.98, 1.01);
    h_b_ege->Draw("P");
    
    h_b_egelinear->SetLineColor(kRed);
    h_b_egelinear->SetMarkerColor(kRed);
    h_b_egelinear->SetMarkerStyle(21);
    h_b_egelinear->Draw("P SAME");
    
    h_b_combined->SetLineColor(kGreen+2);
    h_b_combined->SetMarkerColor(kGreen+2);
    h_b_combined->SetMarkerStyle(22);
    h_b_combined->Draw("P SAME");
    
    TLegend* legend_b = new TLegend(0.65, 0.70, 0.88, 0.88);
    legend_b->SetFillStyle(0);
    legend_b->SetBorderSize(0);
    legend_b->AddEntry(h_b_ege, "EGE", "p");
    legend_b->AddEntry(h_b_egelinear, "EGE Linear", "p");
    legend_b->AddEntry(h_b_combined, "Combined", "p");
    legend_b->Draw();
    
    canvas_b->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/Comparison_B_Fraction.png");
    
    // 关闭文件
    outputFile->Close();
    egeFile->Close();
    egeLinearFile->Close();
    
    std::cout << "Combined results saved to: /eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_combined.root" << std::endl;
    std::cout << "Comparison plots saved as PNG files" << std::endl;
}