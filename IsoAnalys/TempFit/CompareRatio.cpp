#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <iostream>
#include <array>
#include <fstream>
#include <sstream>
#include "../Tool.h"
using namespace AMS_Iso;

// 定义探测器范围结构
struct DetectorRange {
    double min;
    double max;
    int color;
    const char* name;
};

// 定义绘图配置结构
struct DrawConfig {
    const char* xTitle;
    const char* yTitle;
    const char* drawOption;
    std::string savePath;
    double yMin;
    double yMax;
    double legX1;
    double legY1;
    double legX2;
    double legY2;
};

// 定义比率数据结构
struct RatioData {
    double ek_low;
    double ek_high;
    double ratio;
    double ratio_err;
};

// 探测器范围设置
namespace Setting {
    const std::array<DetectorRange, 3> DET_RANGES = {{
        {0.22, 1.5, kBlue, "ToF"},
        {0.8, 3.5, kOrange+1, "NaF"},
        {2.3, 15.0, kGreen+2, "AGL"}
    }};
}

// 读取Wei的比率数据
std::vector<RatioData> ReadRatioFile(const std::string& filename) {
    std::vector<RatioData> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Cannot open ratio file: " << filename << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // 跳过表头

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        RatioData rd;
        ss >> rd.ek_low >> rd.ek_high >> rd.ratio >> rd.ratio_err;
        data.push_back(rd);
    }
    return data;
}

// 从直方图创建比率直方图
TH1F* CreateRatioHist(TH1F* num, TH1F* den, const char* name) {
    TH1F* ratio = (TH1F*)num->Clone(name);
    ratio->Divide(num, den, 1.0, 1.0);
    
    // 计算误差传递
    for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
        double N = num->GetBinContent(i);
        double D = den->GetBinContent(i);
        double dN = num->GetBinError(i);
        double dD = den->GetBinError(i);
        
        if (D > 0) {
            double R = N/D;
            double dR = R * sqrt(pow(dN/N, 2) + pow(dD/D, 2));
            ratio->SetBinError(i, dR);
        } else {
            ratio->SetBinError(i, 0);
        }
    }
    return ratio;
}

// 从分数直方图创建Be10/Be9比率
TH1F* CreateBe10Be9RatioHist(TH1F* be7_hist, TH1F* be9_hist, const char* name) {
    TH1F* ratio = (TH1F*)be9_hist->Clone(name);
    
    for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
        double f7 = be7_hist->GetBinContent(i);
        double f9 = be9_hist->GetBinContent(i);
        double df7 = be7_hist->GetBinError(i);
        double df9 = be9_hist->GetBinError(i);
        
        if (f9 > 0) {
            double be10_frac = 1.0 - f7 - f9;
            double dbe10 = sqrt(df7*df7 + df9*df9);
            double R = be10_frac/f9;
            double dR = R * sqrt(pow(dbe10/be10_frac, 2) + pow(df9/f9, 2));
            
            ratio->SetBinContent(i, R);
            ratio->SetBinError(i, dR);
        } else {
            ratio->SetBinContent(i, 0);
            ratio->SetBinError(i, 0);
        }
    }
    return ratio;
}

void SetHistStyle(TH1F* hist, int color) {
    if (!hist) return;
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerSize(1.2);
}

// 绘制比较图
void DrawRatioComparison(const std::array<TH1F*, 3>& yan_hists, TGraphErrors* wei_graph,
                        const DrawConfig& config, int usemass) {
    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->SetLogx();
    gStyle->SetOptStat(0);

    bool first = true;
    for (int i = 0; i < 3; ++i) {
        if (!yan_hists[i]) continue;
        
        if (first) {
            yan_hists[i]->SetTitle(Form("Mass = %d Binning",usemass));
            yan_hists[i]->GetXaxis()->SetTitle(config.xTitle);
            yan_hists[i]->GetYaxis()->SetTitle(config.yTitle);
            yan_hists[i]->GetYaxis()->SetRangeUser(config.yMin, config.yMax);
            yan_hists[i]->GetXaxis()->SetRangeUser(0.4, 16.0);
            yan_hists[i]->Draw(config.drawOption);
            first = false;
        } else {
            yan_hists[i]->Draw(Form("%s SAME", config.drawOption));
        }
    }

    if (wei_graph) wei_graph->Draw("P SAME");

    TLegend* leg = new TLegend(config.legX1, config.legY1, config.legX2, config.legY2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.035);
    for (int i = 0; i < 3; ++i) {
        if (yan_hists[i]) {
            leg->AddEntry(yan_hists[i], Form("Yan %s", Setting::DET_RANGES[i].name), "p");
        }
    }
    if (wei_graph) leg->AddEntry(wei_graph, "Wei", "p");
    leg->Draw();

    c->SaveAs(config.savePath.c_str());
    delete c;
}

void CompareRatio() {
    const char* baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/";
    std::vector<int> masses = {7, 9, 10};
    
    // 读取Wei的比率数据
    std::vector<RatioData> wei_Be9Be7 = ReadRatioFile(Form("%sWei_ratio_Be9Be7.txt", baseDir));
    std::vector<RatioData> wei_Be10Be9 = ReadRatioFile(Form("%sWei_ratio_Be10Be9.txt", baseDir));
    
    // 创建Wei的图形
    TGraphErrors* wei_graph_9_7 = new TGraphErrors(wei_Be9Be7.size());
    TGraphErrors* wei_graph_10_9 = new TGraphErrors(wei_Be10Be9.size());
    
    for (size_t i = 0; i < wei_Be9Be7.size(); ++i) {
        double ek_center = (wei_Be9Be7[i].ek_low + wei_Be9Be7[i].ek_high) / 2;
        wei_graph_9_7->SetPoint(i, ek_center, wei_Be9Be7[i].ratio);
        wei_graph_9_7->SetPointError(i, 0, wei_Be9Be7[i].ratio_err);
    }
    
    for (size_t i = 0; i < wei_Be10Be9.size(); ++i) {
        double ek_center = (wei_Be10Be9[i].ek_low + wei_Be10Be9[i].ek_high) / 2;
        wei_graph_10_9->SetPoint(i, ek_center, wei_Be10Be9[i].ratio);
        wei_graph_10_9->SetPointError(i, 0, wei_Be10Be9[i].ratio_err);
    }
    
    wei_graph_9_7->SetMarkerStyle(21);
    wei_graph_9_7->SetMarkerColor(kRed);
    wei_graph_9_7->SetLineColor(kRed);
    wei_graph_10_9->SetMarkerStyle(21);
    wei_graph_10_9->SetMarkerColor(kRed);
    wei_graph_10_9->SetLineColor(kRed);

    // 对每个use mass进行处理
    for (int mass : masses) {
        TFile* yan_file = TFile::Open(Form("%sMassTF_Be_DiffBin_r1_NoA_Use%d.root", baseDir, mass));
        
        // 创建比率直方图数组
        std::array<TH1F*, 3> ratio_9_7_hists;
        std::array<TH1F*, 3> ratio_10_9_hists;
        
        for (int j = 0; j < 3; ++j) {
            TH1F* h7 = (TH1F*)yan_file->Get(Form("h_best_Be7_frac_%s", Setting::DET_RANGES[j].name));
            TH1F* h9 = (TH1F*)yan_file->Get(Form("h_best_Be9_frac_%s", Setting::DET_RANGES[j].name));
            
            if (h7 && h9) {
                ratio_9_7_hists[j] = CreateRatioHist(h9, h7, Form("ratio_9_7_%s", Setting::DET_RANGES[j].name));
                ratio_10_9_hists[j] = CreateBe10Be9RatioHist(h7, h9, Form("ratio_10_9_%s", Setting::DET_RANGES[j].name));
                
                SetHistStyle(ratio_9_7_hists[j], Setting::DET_RANGES[j].color);
                SetHistStyle(ratio_10_9_hists[j], Setting::DET_RANGES[j].color);
            } else {
                ratio_9_7_hists[j] = nullptr;
                ratio_10_9_hists[j] = nullptr;
            }
        }

        // 绘制Be9/Be7比率
        DrawConfig config_9_7 = {
            "Ek [GeV/n]",
            "Be9/Be7 Ratio",
            "P",
            Form("%sRatio_Be9Be7_Use%d.png", baseDir, mass),
            0.3, 0.9,
            0.2, 0.68, 0.45, 0.87
        };
        DrawRatioComparison(ratio_9_7_hists, wei_graph_9_7, config_9_7, mass);

        // 绘制Be10/Be9比率
        DrawConfig config_10_9 = {
            "Ek [GeV/n]",
            "Be10/Be9 Ratio",
            "P",
            Form("%sRatio_Be10Be9_Use%d.png", baseDir, mass),
            0.0, 0.6,
            0.2, 0.68, 0.45, 0.87
        };
        DrawRatioComparison(ratio_10_9_hists, wei_graph_10_9, config_10_9, mass);

        // 清理内存
        for (auto hist : ratio_9_7_hists) delete hist;
        for (auto hist : ratio_10_9_hists) delete hist;
        yan_file->Close();
        delete yan_file;
    }

    delete wei_graph_9_7;
    delete wei_graph_10_9;
}