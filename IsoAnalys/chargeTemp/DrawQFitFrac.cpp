#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>

void DrawQFitFrac() {
    // 设置ROOT样式
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    
    // 绘制Beryllium图
    {
        // 打开所有文件
        TFile* file_yyh = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/yyh_chargeFit.root", "READ");
        TFile* file_lg = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_lg1010.root", "READ");
        TFile* file_egelinear = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGELinear_1010.root", "READ");
        TFile* file_ege = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_1010.root", "READ");
        
        // 获取直方图
        TH1F* hist_yyh = nullptr;
        TH1D* hist_lg = nullptr;
        TH1D* hist_egelinear = nullptr;
        TH1D* hist_ege = nullptr;
        
        if (file_yyh) file_yyh->GetObject("hist_ratio_Brelium", hist_yyh);
        if (file_lg) file_lg->GetObject("h_be_fraction", hist_lg);
        if (file_egelinear) file_egelinear->GetObject("h_be_fraction", hist_egelinear);
        if (file_ege) file_ege->GetObject("h_be_fraction", hist_ege);
        
        // 检查是否成功获取所有直方图
        if (!hist_yyh) std::cerr << "Could not find histogram hist_ratio_Brelium in YYH file" << std::endl;
        if (!hist_lg) std::cerr << "Could not find histogram h_be_fraction in LG file" << std::endl;
        if (!hist_egelinear) std::cerr << "Could not find histogram h_be_fraction in EGE Linear file" << std::endl;
        if (!hist_ege) std::cerr << "Could not find histogram h_be_fraction in EGE file" << std::endl;
        
        // 创建画布
        TCanvas* canvas = new TCanvas("canvas_be", "Beryllium Fraction", 900, 700);
        canvas->SetTopMargin(0.07);
        canvas->SetLeftMargin(0.12);
        canvas->SetRightMargin(0.05);
        canvas->SetBottomMargin(0.12);
        
        // 创建图例
        TLegend* legend = new TLegend(0.65, 0.65, 0.88, 0.88);
        legend->SetFillStyle(0);  // 透明填充
        legend->SetBorderSize(0); // 无边框
        legend->SetTextSize(0.035);
        
        // 设置直方图样式并绘制
        double yMin = 0;
        double yMax = 0.01;
        
        // YYH直方图（LG Linear）
        if (hist_yyh) {
            // 创建副本以避免内存问题
            TH1F* hist_yyh_clone = (TH1F*)hist_yyh->Clone("hist_yyh_be");
            hist_yyh_clone->SetLineColor(kRed);
            hist_yyh_clone->SetMarkerColor(kRed);
            hist_yyh_clone->SetMarkerStyle(20);
            hist_yyh_clone->SetMarkerSize(1.0);
            hist_yyh_clone->SetLineWidth(2);
            hist_yyh_clone->SetTitle("InnerQ 4.5-5.5, L1Q 4.8-5.5");
            hist_yyh_clone->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
            hist_yyh_clone->GetYaxis()->SetTitle("Beryllium Fraction");
            hist_yyh_clone->GetYaxis()->SetRangeUser(yMin, yMax);
            hist_yyh_clone->GetXaxis()->SetTitleSize(0.045);
            hist_yyh_clone->GetYaxis()->SetTitleSize(0.045);
            hist_yyh_clone->GetXaxis()->SetTitleOffset(1.0);
            hist_yyh_clone->GetYaxis()->SetTitleOffset(1.2);
            hist_yyh_clone->GetXaxis()->SetLabelSize(0.04);
            hist_yyh_clone->GetYaxis()->SetLabelSize(0.04);
            hist_yyh_clone->Draw("pz");
            canvas->SetLogx(0);  // 确保使用线性X轴
            legend->AddEntry(hist_yyh_clone, "LG Linear", "pe");
        }
        
        // LG直方图
        if (hist_lg) {
            TH1D* hist_lg_clone = (TH1D*)hist_lg->Clone("hist_lg_be");
            hist_lg_clone->SetLineColor(kBlue);
            hist_lg_clone->SetMarkerColor(kBlue);
            hist_lg_clone->SetMarkerStyle(20);
            hist_lg_clone->SetMarkerSize(1.0);
            hist_lg_clone->SetLineWidth(2);
            hist_lg_clone->Draw("pz SAME");
            legend->AddEntry(hist_lg_clone, "LG", "pe");
        }
        
        // EGE Linear直方图
        if (hist_egelinear) {
            TH1D* hist_egelinear_clone = (TH1D*)hist_egelinear->Clone("hist_egelinear_be");
            hist_egelinear_clone->SetLineColor(kGreen+2);
            hist_egelinear_clone->SetMarkerColor(kGreen+2);
            hist_egelinear_clone->SetMarkerStyle(20);
            hist_egelinear_clone->SetMarkerSize(1.0);
            hist_egelinear_clone->SetLineWidth(2);
            hist_egelinear_clone->Draw("pz SAME");
            legend->AddEntry(hist_egelinear_clone, "EGE Linear", "pe");
        }
        
        // EGE直方图
        if (hist_ege) {
            TH1D* hist_ege_clone = (TH1D*)hist_ege->Clone("hist_ege_be");
            hist_ege_clone->SetLineColor(kMagenta);
            hist_ege_clone->SetMarkerColor(kMagenta);
            hist_ege_clone->SetMarkerStyle(20);
            hist_ege_clone->SetMarkerSize(1.0);
            hist_ege_clone->SetLineWidth(2);
            hist_ege_clone->Draw("pz SAME");
            legend->AddEntry(hist_ege_clone, "EGE", "pe");
        }
        
        // 绘制图例
        legend->Draw();
        
        // 保存画布
        canvas->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/QFitFrac_h_be_fraction.png");
        
        // 关闭文件
        if (file_yyh) file_yyh->Close();
        if (file_lg) file_lg->Close();
        if (file_egelinear) file_egelinear->Close();
        if (file_ege) file_ege->Close();
        
        std::cout << "Created plot for Beryllium Fraction" << std::endl;
    }
    
    // 绘制Boron图
    {
        // 打开所有文件
        TFile* file_yyh = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/yyh_chargeFit.root", "READ");
        TFile* file_lg = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_lg1010.root", "READ");
        TFile* file_egelinear = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGELinear_1010.root", "READ");
        TFile* file_ege = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_1010.root", "READ");
        
        // 获取直方图
        TH1F* hist_yyh = nullptr;
        TH1D* hist_lg = nullptr;
        TH1D* hist_egelinear = nullptr;
        TH1D* hist_ege = nullptr;
        
        if (file_yyh) file_yyh->GetObject("hist_ratio_Boron", hist_yyh);
        if (file_lg) file_lg->GetObject("h_b_fraction", hist_lg);
        if (file_egelinear) file_egelinear->GetObject("h_b_fraction", hist_egelinear);
        if (file_ege) file_ege->GetObject("h_b_fraction", hist_ege);
        
        // 检查是否成功获取所有直方图
        if (!hist_yyh) std::cerr << "Could not find histogram hist_ratio_Boron in YYH file" << std::endl;
        if (!hist_lg) std::cerr << "Could not find histogram h_b_fraction in LG file" << std::endl;
        if (!hist_egelinear) std::cerr << "Could not find histogram h_b_fraction in EGE Linear file" << std::endl;
        if (!hist_ege) std::cerr << "Could not find histogram h_b_fraction in EGE file" << std::endl;
        
        // 创建画布
        TCanvas* canvas = new TCanvas("canvas_b", "Boron Fraction", 900, 700);
        canvas->SetTopMargin(0.07);
        canvas->SetLeftMargin(0.12);
        canvas->SetRightMargin(0.05);
        canvas->SetBottomMargin(0.12);
        
        // 创建图例
        TLegend* legend = new TLegend(0.65, 0.65, 0.88, 0.88);
        legend->SetFillStyle(0);  // 透明填充
        legend->SetBorderSize(0); // 无边框
        legend->SetTextSize(0.035);
        
        // 设置直方图样式并绘制
        double yMin = 0.97;
        double yMax = 1.02;
        
        // YYH直方图（LG Linear）
        if (hist_yyh) {
            // 创建副本以避免内存问题
            TH1F* hist_yyh_clone = (TH1F*)hist_yyh->Clone("hist_yyh_b");
            hist_yyh_clone->SetLineColor(kRed);
            hist_yyh_clone->SetMarkerColor(kRed);
            hist_yyh_clone->SetMarkerStyle(20);
            hist_yyh_clone->SetMarkerSize(1.0);
            hist_yyh_clone->SetLineWidth(2);
            hist_yyh_clone->SetTitle("InnerQ 4.5-5.5, L1Q 4.8-5.5");
            hist_yyh_clone->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
            hist_yyh_clone->GetYaxis()->SetTitle("Boron Fraction");
            hist_yyh_clone->GetYaxis()->SetRangeUser(yMin, yMax);
            hist_yyh_clone->GetXaxis()->SetTitleSize(0.045);
            hist_yyh_clone->GetYaxis()->SetTitleSize(0.045);
            hist_yyh_clone->GetXaxis()->SetTitleOffset(1.0);
            hist_yyh_clone->GetYaxis()->SetTitleOffset(1.2);
            hist_yyh_clone->GetXaxis()->SetLabelSize(0.04);
            hist_yyh_clone->GetYaxis()->SetLabelSize(0.04);
            hist_yyh_clone->Draw("pz");
            canvas->SetLogx(0);  // 确保使用线性X轴
            legend->AddEntry(hist_yyh_clone, "LG Linear", "pe");
        }
        
        // LG直方图
        if (hist_lg) {
            TH1D* hist_lg_clone = (TH1D*)hist_lg->Clone("hist_lg_b");
            hist_lg_clone->SetLineColor(kBlue);
            hist_lg_clone->SetMarkerColor(kBlue);
            hist_lg_clone->SetMarkerStyle(20);
            hist_lg_clone->SetMarkerSize(1.0);
            hist_lg_clone->SetLineWidth(2);
            hist_lg_clone->Draw("pz SAME");
            legend->AddEntry(hist_lg_clone, "LG", "pe");
        }
        
        // EGE Linear直方图
        if (hist_egelinear) {
            TH1D* hist_egelinear_clone = (TH1D*)hist_egelinear->Clone("hist_egelinear_b");
            hist_egelinear_clone->SetLineColor(kGreen+2);
            hist_egelinear_clone->SetMarkerColor(kGreen+2);
            hist_egelinear_clone->SetMarkerStyle(20);
            hist_egelinear_clone->SetMarkerSize(1.0);
            hist_egelinear_clone->SetLineWidth(2);
            hist_egelinear_clone->Draw("pz SAME");
            legend->AddEntry(hist_egelinear_clone, "EGE Linear", "pe");
        }
        
        // EGE直方图
        if (hist_ege) {
            TH1D* hist_ege_clone = (TH1D*)hist_ege->Clone("hist_ege_b");
            hist_ege_clone->SetLineColor(kMagenta);
            hist_ege_clone->SetMarkerColor(kMagenta);
            hist_ege_clone->SetMarkerStyle(20);
            hist_ege_clone->SetMarkerSize(1.0);
            hist_ege_clone->SetLineWidth(2);
            hist_ege_clone->Draw("pz SAME");
            legend->AddEntry(hist_ege_clone, "EGE", "pe");
        }
        
        // 绘制图例
        legend->Draw();
        
        // 保存画布
        canvas->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/QFitFrac_h_b_fraction.png");
        
        // 关闭文件
        if (file_yyh) file_yyh->Close();
        if (file_lg) file_lg->Close();
        if (file_egelinear) file_egelinear->Close();
        if (file_ege) file_ege->Close();
        
        std::cout << "Created plot for Boron Fraction" << std::endl;
    }
    
    std::cout << "All plots created successfully!" << std::endl;
}