#ifndef HIST_TEMP_FIT_H
#define HIST_TEMP_FIT_H

#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include "RooAbsReal.h"
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <memory>
#include "../Tool.h"

struct IsotopeConfig {
    std::string name;           // Li, Be, B等
    double charge;
    std::string dataName;       // Lit, Ber, Bor
    std::vector<int> masses;    // 质量数
    std::vector<double> initFractions; // 初始分数
    std::vector<double> fitRangeLow;   // ToF, NaF, AGL的拟合下限
    std::vector<double> fitRangeUp;    // ToF, NaF, AGL的拟合上限
    int nBins;                  // bin数量
    double detRanges[2];        // TOF最大bin和AGL最小bin
    int nFitParams;             // 拟合参数数量
};

struct FitCacheEntry {
    double chi2_ndf;
    std::vector<double> fractions;
    std::vector<double> fraction_errors;
    double alpha;
    RooPlot* frame;
    std::vector<TH1D*> mc_hists;
    TH1D* data_hist;
    double significance;
    double signal_efficiency;
    double upperBound;
    double n_signal;
    double n_total;
    double n_entries;
    bool is_valid;
    
    FitCacheEntry() : 
        chi2_ndf(5000), alpha(0), frame(nullptr), data_hist(nullptr), 
        significance(0), signal_efficiency(0), upperBound(0), 
        n_signal(0), n_total(0), n_entries(0), is_valid(false) {}
    
    ~FitCacheEntry() {
        if (frame) delete frame;
        if (data_hist) delete data_hist;
        for (auto* hist : mc_hists) {
            if (hist) delete hist;
        }
    }
    
    FitCacheEntry(const FitCacheEntry&) = delete;
    FitCacheEntry& operator=(const FitCacheEntry&) = delete;
    // 添加移动构造函数和移动赋值运算符
    FitCacheEntry(FitCacheEntry&&) = default;
    FitCacheEntry& operator=(FitCacheEntry&&) = default;
};
class FitCache {
private:
    std::map<int, std::unique_ptr<FitCacheEntry>> best_fits;  // 键为能量bin

public:
    void UpdateCache(int ek_bin, double chi2_ndf, double alpha, RooFitResult* fit_res,
                    const std::vector<RooRealVar*>& fractions,
                    RooPlot* frame, const std::vector<TH1D*>& mc_hists, TH1D* data_hist, 
                    double significance = 0, double signal_efficiency = 0, double upperBound = 0,
                    double n_signal = 0, double n_total = 0, double n_entries = 0) {
        
        if (best_fits.find(ek_bin) == best_fits.end() || 
            chi2_ndf < best_fits[ek_bin]->chi2_ndf) {
            
            auto new_entry = std::make_unique<FitCacheEntry>();
            new_entry->chi2_ndf = chi2_ndf;
            new_entry->alpha = alpha;
            new_entry->is_valid = (fit_res->status() == 0);
            
            // 保存fraction结果
            new_entry->fractions.resize(fractions.size());
            new_entry->fraction_errors.resize(fractions.size());
            for (size_t i = 0; i < fractions.size(); ++i) {
                new_entry->fractions[i] = fractions[i]->getVal();
                new_entry->fraction_errors[i] = fractions[i]->getError();
            }
            
            // 深拷贝frame和直方图
            new_entry->frame = static_cast<RooPlot*>(frame->Clone());
            new_entry->data_hist = static_cast<TH1D*>(data_hist->Clone());
            for (const auto* hist : mc_hists) {
                new_entry->mc_hists.push_back(static_cast<TH1D*>(hist->Clone()));
            }
            
            // 保存其他结果
            new_entry->significance = significance;
            new_entry->signal_efficiency = signal_efficiency;
            new_entry->upperBound = upperBound;
            new_entry->n_signal = n_signal;
            new_entry->n_total = n_total;
            new_entry->n_entries = n_entries; 
            
            // 更新缓存
            best_fits[ek_bin] = std::move(new_entry);
        }
    }
    
    const FitCacheEntry* GetBestFit(int ek_bin) const {
        auto it = best_fits.find(ek_bin);
        return (it != best_fits.end()) ? it->second.get() : nullptr;
    }
    
    bool HasCache(int ek_bin) const {
        return best_fits.find(ek_bin) != best_fits.end();
    }
};

class IsoFitConstants {
public:
    static const std::map<std::string, IsotopeConfig> configs;
    static const char* x_titles[];
    static const char* x_titlesMC[];
};

void HistTempFit(const std::string& isotype = "Be", 
                 const std::string& suffix = "_",
                 const std::string& inputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData",
                 const std::string& outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit");

void HistTempFit_hist();

const char* IsoFitConstants::x_titles[] = {
    "TOFMass",
    "NaFMass",
    "AGLMass"
};
const char* IsoFitConstants::x_titlesMC[] = {
    "ToFMass",
    "NaFMass",
    "AGLMass"
};

const std::map<std::string, IsotopeConfig> IsoFitConstants::configs = {
    {"Li", {"Li", 3, "Lit", {6, 7}, {0.5, 0.5}, 
           {0.095, 0.095, 0.095}, {0.24, 0.242, 0.242},74, {18, 17}, 1}},

    {"Be", {"Be", 4, "Ber", {7, 9, 10}, {0.6, 0.33, 0.07}, 
           //{0.059, 0.063, 0.06}, {0.215, 0.217, 0.213}, 74, {18, 17}, 2}},//bkg
           //{0.06, 0.065, 0.062}, {0.205, 0.213, 0.21}, 74, {18, 17}, 2}},//bkg
           //{0.065, 0.068, 0.065}, {0.215, 0.216, 0.21}, 74, {18, 17}, 2}},
           {0.05, 0.05, 0.05}, {0.23, 0.23, 0.23}, 74, {18, 17}, 2}},//sdiat

    {"B", {"B", 5, "Bor", {10, 11}, {0.3, 0.7}, 
          {0.05, 0.05, 0.05}, {0.15, 0.15, 0.15},74, {18, 17}, 1}}
};

inline void setLegend(TLegend *legend) {
    legend->SetBorderSize(0);    // 移除边框
    legend->SetFillColor(0);     // 设置填充颜色为透明
    legend->SetFillStyle(0);     // 设置填充样式为透明
    legend->SetTextSize(0.035);  // 设置文字大小
}

inline void setPaveText(TPaveText *pt) {
    pt->SetBorderSize(0);    // 移除边框
    pt->SetFillColor(0);     // 设置填充颜色为透明
    pt->SetFillStyle(0);     // 设置填充样式为透明
    pt->SetTextSize(0.03);   // 设置文字大小
}

inline void setPullGraphStyle(TGraphErrors *pullGraph) {
    pullGraph->SetMarkerColor(kBlack);    // 黑色标记
    pullGraph->SetMarkerStyle(20);        // 实心圆点

    pullGraph->GetXaxis()->SetLabelSize(.12);     // X轴标签大小
    pullGraph->GetXaxis()->SetTitleSize(.15);     // X轴标题大小
    pullGraph->GetXaxis()->SetTitleOffset(.7);    // X轴标题偏移

    pullGraph->GetYaxis()->SetLabelSize(.09);     // Y轴标签大小
    pullGraph->GetYaxis()->SetTitleSize(.15);     // Y轴标题大小
    pullGraph->GetYaxis()->SetTitleOffset(.4);    // Y轴标题偏移
    pullGraph->GetYaxis()->SetNdivisions(505);    // Y轴刻度分隔
}

#endif // HIST_TEMP_FIT_H