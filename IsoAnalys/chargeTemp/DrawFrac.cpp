#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>

void DrawFrac() {
    // 设置ROOT样式
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);  // 不显示标题
    
    // 打开输入文件
    TFile* file = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_combined.root", "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    
    // 获取直方图
    TH1D* h_be_fraction = nullptr;
    TH1D* h_b_fraction = nullptr;
    TH1D* h_c_fraction = nullptr;
    
    file->GetObject("h_be_fraction", h_be_fraction);
    file->GetObject("h_b_fraction", h_b_fraction);
    file->GetObject("h_c_fraction", h_c_fraction);
    
    if (!h_be_fraction || !h_b_fraction || !h_c_fraction) {
        std::cerr << "Error: missing histograms in file" << std::endl;
        file->Close();
        return;
    }
    
    // 通用画布设置函数
    auto setupCanvas = [](TCanvas* canvas) {
        canvas->SetLeftMargin(0.1);
        canvas->SetRightMargin(0.03);
    };
    
    // 通用直方图设置函数
    auto setupHistogram = [](TH1D* hist) {
        // 去掉X和Y轴标题
        hist->GetXaxis()->SetTitle("");
        hist->GetYaxis()->SetTitle("");
        
        // 调大标签大小
        hist->GetXaxis()->SetLabelSize(0.1);
        hist->GetYaxis()->SetLabelSize(0.1);
        
        // 设置Y轴刻度分隔为505
        hist->GetYaxis()->SetNdivisions(505);
    };
    
    // 绘制Be分数
    TCanvas* canvas_be = new TCanvas("canvas_be", "Be Fraction", 900, 300);
    setupCanvas(canvas_be);
    
    h_be_fraction->SetLineColor(kBlue);
    h_be_fraction->SetMarkerColor(kBlue);
    h_be_fraction->SetMarkerStyle(20);
    h_be_fraction->SetMarkerSize(0.9);
    
    setupHistogram(h_be_fraction);
    
    // 设置Be的Y轴范围为0-0.006
    h_be_fraction->GetYaxis()->SetRangeUser(0, 0.006);
    
    h_be_fraction->Draw("pz");
    
    canvas_be->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/Be_Fraction.png");
    
    // 绘制B分数
    TCanvas* canvas_b = new TCanvas("canvas_b", "B Fraction", 900, 300);
    setupCanvas(canvas_b);
    
    h_b_fraction->SetLineColor(kGreen+2);
    h_b_fraction->SetMarkerColor(kGreen+2);
    h_b_fraction->SetMarkerStyle(20);
    h_b_fraction->SetMarkerSize(0.9);
    
    setupHistogram(h_b_fraction);
    
    // 设置B的Y轴范围
    h_b_fraction->GetYaxis()->SetRangeUser(0.985, 1.005);
    
    h_b_fraction->Draw("pz");
    
    canvas_b->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/B_Fraction.png");
    
    // 绘制C分数
    TCanvas* canvas_c = new TCanvas("canvas_c", "C Fraction", 900, 300);
    setupCanvas(canvas_c);
    
    h_c_fraction->SetLineColor(kMagenta);
    h_c_fraction->SetMarkerColor(kMagenta);
    h_c_fraction->SetMarkerStyle(20);
    h_c_fraction->SetMarkerSize(0.9);
    
    setupHistogram(h_c_fraction);
    
    // 设置C的Y轴范围
    h_c_fraction->GetYaxis()->SetRangeUser(0, 0.003);
    
    h_c_fraction->Draw("pz");
    
    canvas_c->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/C_Fraction.png");
    
    // 关闭文件
    file->Close();
    
    std::cout << "Fraction plots saved as PNG files" << std::endl;
}