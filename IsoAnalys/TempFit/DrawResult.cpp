#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <array>
#include <map>
#include "../Tool.h"
using namespace AMS_Iso;

// Isotope configuration struct
struct IsotopeConfig {
    std::string dataName;      
    std::vector<std::string> fractionNames;
    std::vector<std::pair<double, double>> fractionRanges;
    std::vector<int> masses;  // 存储所有同位素质量数
};

// Define detector range structure
struct DetectorRange {
    double min;
    double max;
    int color;
    const char* name;
};

// Global Setting and configurations
namespace Setting {
    // Detector ranges configuration
    const std::array<DetectorRange, 3> DET_RANGES = {{
        {0.45, 1.4, kBlue, "ToF"},
        {.9, 4, kOrange+1, "NaF"},
        {2.3, 16.5, kGreen+2, "AGL"}
    }};

    // Isotope configurations
    const std::map<int, IsotopeConfig> ISO_CONFIGS = {
        {3, {
            "Li",
            {"Li6", "Li7"},
            {{0.4, 0.65}, {0.35, 0.6}},
            {6, 7}
        }},
        {4, {
            "Be",
            {"Be7", "Be9", "Be10"},
            {{0.45, 0.7}, {0.15, 0.5}, {0.0, 0.15}},
            {7, 9, 10}
        }},
        {5, {
            "B",
            {"B10", "B11"},
            {{0.1, 0.5}, {0.5, 0.9}},
            {10, 11}
        }}
    };
}

struct DrawConfig {
    const char* xTitle;
    const char* yTitle;
    const char* drawOption;
    std::string savePath;
    bool useErrors;
    double yMin;
    double yMax;
    double legX1;
    double legY1;
    double legX2;
    double legY2;
};

void FillNewHist(TH1F* src, TH1F* dst, double xMin, double xMax, bool useErrors) {
    if (!src || !dst) {
        std::cout << "Null histogram pointer!" << std::endl;
        return;
    }
    
    for (int i = 1; i <= src->GetNbinsX(); ++i) {
        double binLowEdge = src->GetBinLowEdge(i);
        double binUpEdge = src->GetBinLowEdge(i) + src->GetBinWidth(i);
        double binCenter = src->GetBinCenter(i);
        double content = src->GetBinContent(i);
        
        if (binLowEdge < xMax && binUpEdge > xMin) {
            int dstBin = dst->FindBin(binCenter);
            if (dstBin > 0 && dstBin <= dst->GetNbinsX()) {
                dst->SetBinContent(dstBin, content);
                if (useErrors) {
                    dst->SetBinError(dstBin, src->GetBinError(i));
                }
            }
        }
    }
}

void SetHistStyle(TH1F* hist, int color) {
    if (!hist) return;
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerSize(1.2);
}
template<size_t N>
void DrawHist(const std::array<TH1F*, 3>& hists, const DrawConfig& config, 
              const std::array<double, N>& ek_bins) {
    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->SetLogx();
    gStyle->SetOptStat(0);

    double globalXMin = Setting::DET_RANGES[0].min;
    double globalXMax = 0;
    for (const auto& range : Setting::DET_RANGES) {
        if (range.max > globalXMax) globalXMax = range.max;
    }
    
    std::array<TH1F*, 3> newHists;
    for (int i = 0; i < 3; ++i) {
        if (!hists[i]) continue;
        
        newHists[i] = new TH1F(Form("h_%s_new", Setting::DET_RANGES[i].name),
                              "", ek_bins.size()-1, ek_bins.data());
        FillNewHist(hists[i], newHists[i], 
                   Setting::DET_RANGES[i].min, 
                   Setting::DET_RANGES[i].max, 
                   config.useErrors);
        SetHistStyle(newHists[i], Setting::DET_RANGES[i].color);
    }

    bool first = true;
    for (int i = 0; i < 3; ++i) {
        if (!newHists[i]) continue;
        
        if (first) {
            newHists[i]->GetYaxis()->SetNdivisions(505);
            newHists[i]->GetXaxis()->SetTitle(config.xTitle);
            newHists[i]->GetYaxis()->SetTitle(config.yTitle);
            if (config.yMin < config.yMax) {
                newHists[i]->GetYaxis()->SetRangeUser(config.yMin, config.yMax);
            }
            newHists[i]->GetXaxis()->SetRangeUser(globalXMin, globalXMax);
            newHists[i]->Draw(config.drawOption);
            first = false;
        } else {
            newHists[i]->Draw(Form("%s SAME", config.drawOption));
        }
    }

    TLegend* leg = new TLegend(config.legX1, config.legY1, config.legX2, config.legY2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.035);
    for (int i = 0; i < 3; ++i) {
        if (newHists[i]) {
            leg->AddEntry(newHists[i], Setting::DET_RANGES[i].name, "p");
        }
    }
    leg->Draw();

    c->SaveAs(config.savePath.c_str());

    for (auto hist : newHists) delete hist;
    delete c;
}

void DrawTFresult(int charge, int mass) {
    const auto& config = Setting::ISO_CONFIGS.at(charge);
    const auto& ek_bins = getKineticEnergyBins(charge, mass);
    
    const char* baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/";
    
    TFile* infile = TFile::Open(Form("%sMassTF_%s_DiffBin_NoA_NewMCISS_Use%d.root", baseDir, config.dataName.c_str(), mass));

    if (!infile) {
        std::cerr << "Error: Unable to open one or more files." << std::endl;
        return;
    }

    // 获取拟合的分数直方图
    std::vector<std::array<TH1F*, 3>> fractionHists;
    for (size_t i = 0; i < config.fractionNames.size() - 1; ++i) {
        std::array<TH1F*, 3> hists;
        for (int j = 0; j < 3; ++j) {
            hists[j] = (TH1F*)infile->Get(Form("h_best_%s_frac_%s", config.fractionNames[i].c_str(), Setting::DET_RANGES[j].name));
            cout<<hists[j]->GetMaximum()<<endl;
        }
        fractionHists.push_back(hists);
    }

    // 计算剩余同位素的比例
    std::array<TH1F*, 3> remainingHists;
    for (int i = 0; i < 3; ++i) {
        if (fractionHists[0][i]) {
            remainingHists[i] = (TH1F*)fractionHists[0][i]->Clone(
                Form("h_best_%s_frac_%d", config.fractionNames.back().c_str(), i));
            
            // 对每个bin计算剩余比例和误差
            for (int bin = 1; bin <= remainingHists[i]->GetNbinsX(); ++bin) {
                double total_frac = 0;
                double total_err2 = 0;
                
                // 累加已有比例
                for (size_t j = 0; j < fractionHists.size(); ++j) {
                    total_frac += fractionHists[j][i]->GetBinContent(bin);
                    double err = fractionHists[j][i]->GetBinError(bin);
                    total_err2 += err * err;
                }
                
                // 计算剩余比例和误差
                double remaining_frac = 1.0 - total_frac;
                double remaining_err = sqrt(total_err2);
                
                remainingHists[i]->SetBinContent(bin, remaining_frac);
                remainingHists[i]->SetBinError(bin, remaining_err);
            }
        } else {
            remainingHists[i] = nullptr;
        }
    }
    fractionHists.push_back(remainingHists);

    // 获取chi2直方图
    std::array<TH1F*, 3> chi2Hists;
    for (int i = 0; i < 3; ++i) {
        chi2Hists[i] = (TH1F*)infile->Get(Form("h_best_chi2_%s",Setting::DET_RANGES[i].name));
    }

    // 绘制所有分数
    for (size_t i = 0; i < config.fractionNames.size(); ++i) {
        DrawConfig fracConfig = {
            "Ek [GeV/n]", 
            Form("%s Fraction", config.fractionNames[i].c_str()),
            "E",
            Form("%sMassTF_%s_DiffBin_NoA_NewMCISS_Use%d_%s.png", baseDir, 
                config.dataName.c_str(), mass, config.fractionNames[i].c_str()),
            true,
            config.fractionRanges[i].first,
            config.fractionRanges[i].second,
            0.2, 0.65, 0.35, 0.85
        };
        DrawHist(fractionHists[i], fracConfig, ek_bins);
    }

    // 绘制chi2
    DrawConfig chi2Config = {
        "Ek [GeV/n]", "Chi2/NDF", "P",
        Form("%sMassTF_%s_DiffBin_NoA_NewMCISS_Use%d_Chi2.png", baseDir, config.dataName.c_str(), mass),
        false, 0, 3, 0.45, 0.65, 0.65, 0.85
    };
    DrawHist(chi2Hists, chi2Config, ek_bins);

        infile->Close();
        delete infile;
}

void DrawResult()
{
    DrawTFresult(4, 7);  // 示例：Be, mass=7
    //DrawTFresult(4, 9);  // 示例：Be, mass=7
    //DrawTFresult(4, 10);  // 示例：Be, mass=7
}