#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <array>
#include <map>
#include <fstream>
#include <sstream>
#include "../Tool.h"
using namespace AMS_Iso;

// 结构体定义
struct IsotopeConfig {
    std::string dataName;      
    std::vector<std::string> fractionNames;
    std::vector<std::pair<double, double>> fractionRanges;
};

struct DetectorRange {
    double min;
    double max;
    int color;
    const char* name;
};

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

struct RatioData {
    double ek_low;
    double ek_high;
    double ratio;
    double ratio_err;
};

// Global constants
namespace Setting {
    const std::array<DetectorRange, 3> DET_RANGES = {{
        {0.22, 1.5, kBlue, "ToF"},
        {0.8, 3.5, kOrange+1, "NaF"},
        {2.3, 15.0, kGreen+2, "AGL"}
    }};
}

// Functions for reading data
TGraphErrors* CreateGraphFromFile(const std::string& filename) {
    std::vector<double> x, y, yerr;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Cannot open file: " << filename << std::endl;
        return nullptr;
    }

    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> values;

        while (std::getline(ss, value, ',')) {
            value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());
            values.push_back(std::stod(value));
        }

        if (values.size() >= 4) {
            double ek_center = (values[0] + values[1])/2;
            x.push_back(ek_center);
            y.push_back(values[2]);
            yerr.push_back(values[3]);
        }
    }

    TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), 
                                         nullptr, yerr.data());
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetMarkerSize(1.);

    return graph;
}

std::vector<RatioData> ReadRatioFile(const std::string& filename) {
    std::vector<RatioData> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Cannot open ratio file: " << filename << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        RatioData rd;
        ss >> rd.ek_low >> rd.ek_high >> rd.ratio >> rd.ratio_err;
        data.push_back(rd);
    }
    return data;
}

TGraphErrors* CreateFractionGraphFromRatios(const std::vector<RatioData>& Be9Be7_data, 
                                    const std::vector<RatioData>& Be10Be9_data,
                                    const std::string& fraction_type) {
    std::vector<double> x, y, yerr;
    
    for (size_t i = 0; i < Be9Be7_data.size(); ++i) {
        double ek_center = (Be9Be7_data[i].ek_low + Be9Be7_data[i].ek_high) / 2;
        
        double R97 = Be9Be7_data[i].ratio;
        double R109 = Be10Be9_data[i].ratio;
        double R97_err = Be9Be7_data[i].ratio_err;
        double R109_err = Be10Be9_data[i].ratio_err;
        
        double f_Be7 = 1.0/(1.0 + R97 + R97*R109);
        double f_Be9 = R97 * f_Be7;
        double f_Be10 = R109 * f_Be9;
        
        // Be7误差
        double df7_dR97 = -f_Be7*f_Be7*(1 + R109);
        double df7_dR109 = -f_Be7*f_Be7*R97;
        double sigma_Be7 = sqrt(pow(df7_dR97*R97_err, 2) + pow(df7_dR109*R109_err, 2));

        // Be9误差
        double df9_dR97 = f_Be7 + R97*df7_dR97;
        double df9_dR109 = R97*df7_dR109;
        double sigma_Be9 = sqrt(pow(df9_dR97*R97_err, 2) + pow(df9_dR109*R109_err, 2));
        
        // Be10误差
        double df10_dR97 = R109*df9_dR97;
        double df10_dR109 = f_Be9 + R109*df9_dR109;
        double sigma_Be10 = sqrt(pow(df10_dR97*R97_err, 2) + pow(df10_dR109*R109_err, 2));
        
        x.push_back(ek_center);
        if (fraction_type == "Be7") {
            y.push_back(f_Be7);
            yerr.push_back(sigma_Be7);
        } else if (fraction_type == "Be9") {
            y.push_back(f_Be9);
            yerr.push_back(sigma_Be9);
        } else if (fraction_type == "Be10") {
            y.push_back(f_Be10);
            yerr.push_back(sigma_Be10);
        }
    }
    
    TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), 
                                         nullptr, yerr.data());
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kRed);
    graph->SetMarkerSize(1.);
    
    return graph;
}

void FillNewHist(TH1F* src, TH1F* dst, double xMin, double xMax, bool useErrors) {
	if (!src || !dst) {
		std::cout << "Null histogram pointer!" << std::endl;
		return;
	}

	int filledBins = 0;
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
				filledBins++;
			}
		}
	}
}

TH1F* CalculateRemainingFraction(TH1F* be7_hists, TH1F* be9_hists, int j) {
    TH1F* be10_hist;
        if (be7_hists && be9_hists) {
            be10_hist = (TH1F*)be7_hists->Clone(Form("h_best_Be10_frac_%d",j));
            
            for (int bin = 1; bin <= be10_hist->GetNbinsX(); ++bin) {
                double be7_frac = be7_hists->GetBinContent(bin);
                double be9_frac = be9_hists->GetBinContent(bin);
                double be7_err = be7_hists->GetBinError(bin);
                double be9_err = be9_hists->GetBinError(bin);

                double be10_frac = 1.0 - be7_frac - be9_frac;
                double be10_err = sqrt(be7_err*be7_err + be9_err*be9_err);

                be10_hist->SetBinContent(bin, be10_frac);
                be10_hist->SetBinError(bin, be10_err);
            }
        } else {
            be10_hist = nullptr;
        }
    return be10_hist;
}

void SetHistStyle(TH1F* hist, int color) {
	if (!hist) return;
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(color);
	hist->SetLineColor(color);
	hist->SetMarkerSize(1.2);
}

template<size_t N>
void DrawComparisonHist(const std::array<TH1F*, 3>& yan_hists, TGraphErrors* wei_graph, TGraphErrors* geneva_graph,
                       const DrawConfig& config, const std::array<double, N>& ek_bins, int usemass) {
    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->SetLogx();
    gStyle->SetOptStat(0);

    // 获取合适的全局范围
    double globalXMin = Setting::DET_RANGES[0].min;
    double globalXMax = 0;
    for (const auto& range : Setting::DET_RANGES) {
        if (range.max > globalXMax) globalXMax = range.max;
    }

    // 处理Yan的直方图
    std::array<TH1F*, 3> newHists;
    for (int i = 0; i < 3; ++i) {
        if (!yan_hists[i]) continue;
        
        newHists[i] = new TH1F(Form("h_%s_new", Setting::DET_RANGES[i].name),
                              "", ek_bins.size()-1, ek_bins.data());
        FillNewHist(yan_hists[i], newHists[i], 
                   (i==0 && usemass == 7) ? 0.42 : Setting::DET_RANGES[i].min, 
                   Setting::DET_RANGES[i].max, 
                   config.useErrors);
        SetHistStyle(newHists[i], Setting::DET_RANGES[i].color);
    }

    // 绘图
    bool first = true;
    for (int i = 0; i < 3; ++i) {
        if (!newHists[i]) continue;
        
        if (first) {
            newHists[i]->GetXaxis()->SetTitle(config.xTitle);
            newHists[i]->GetYaxis()->SetTitle(config.yTitle);
            if (config.yMin < config.yMax) {
                newHists[i]->GetYaxis()->SetRangeUser(config.yMin, config.yMax);
            }
            newHists[i]->GetXaxis()->SetRangeUser(usemass == 7 ? 0.4 : globalXMin, 1.1*globalXMax);
            newHists[i]->Draw(config.drawOption);
            first = false;
        } else {
            newHists[i]->Draw(Form("%s SAME", config.drawOption));
        }
    }

    // 绘制Wei和Geneva的结果
    if (wei_graph) wei_graph->Draw("P SAME");
    if (geneva_graph) geneva_graph->Draw("P SAME");

    // 图例
    TLegend* leg = new TLegend(config.legX1, config.legY1, config.legX2, config.legY2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.035);
    for (int i = 0; i < 3; ++i) {
        if (newHists[i]) {
            leg->AddEntry(newHists[i], Form("Yan %s", Setting::DET_RANGES[i].name), "p");
        }
    }
    if (wei_graph) leg->AddEntry(wei_graph, "Wei", "p");
    if (geneva_graph) leg->AddEntry(geneva_graph, "Geneva", "p");
    leg->Draw();

    c->SaveAs(config.savePath.c_str());

    for (auto hist : newHists) delete hist;
    delete c;
}

void CompareBeResults() {
    const char* baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/";

    // Read Wei's ratios
    std::vector<RatioData> Be9Be7_data = ReadRatioFile(Form("%sWei_ratio_Be9Be7.txt", baseDir));
    std::vector<RatioData> Be10Be9_data = ReadRatioFile(Form("%sWei_ratio_Be10Be9.txt", baseDir));

    std::vector<std::string> isotopes = {"Be7", "Be9", "Be10"};
    std::vector<int> masses = {7, 9, 10};
	double yAxismin[] = {0.45, 0.15, 0.};
	double yAxismax[] = {0.7 , 0.5 , 0.15};
	
    for (size_t i = 0; i < isotopes.size(); ++i) {
		cout<<masses[i]<<endl;
        // Read Yan's results
        TFile* yan_file = TFile::Open(Form("%sMassTF_Be_DiffBin_r2_NoA_Use%d.root", 
                                         baseDir, masses[i]));
        
    	const auto& ek_bins = getKineticEnergyBins(4, masses[i]); // For Be
        std::array<TH1F*, 3> yan_hists;
        for (int j = 0; j < 3; ++j) {
			if(i < 2){
				yan_hists[j] = (TH1F*)yan_file->Get(Form("h_best_%s_frac_%s", 
							isotopes[i].c_str(), Setting::DET_RANGES[j].name));
			}
			else{
				TH1F *h7 = (TH1F*)yan_file->Get(Form("h_best_%s_frac_%s", "Be7", Setting::DET_RANGES[j].name));
				TH1F *h9 = (TH1F*)yan_file->Get(Form("h_best_%s_frac_%s", "Be9", Setting::DET_RANGES[j].name));
				yan_hists[j] = CalculateRemainingFraction(h7, h9, j);
			}
			cout<<yan_hists[j]->GetMaximum()<<endl;
        }

        // Read Wei's results
        TGraphErrors* wei_graph = CreateFractionGraphFromRatios(Be9Be7_data, Be10Be9_data, isotopes[i]);

        // Read Geneva's results
        TGraphErrors* geneva_graph = CreateGraphFromFile(Form("%sGeneva_%s_FluxFraction.txt", 
                                             baseDir, isotopes[i].c_str()));

        // Draw comparison
        DrawConfig config = {
            "Ek [GeV/n]", 
            Form("%s Fraction", isotopes[i].c_str()),
            "P",
            Form("%sComparison_r2_%s.png", baseDir, isotopes[i].c_str()),
            true,
            yAxismin[i], yAxismax[i],
            0.2, 0.68, 0.45, 0.87
        };

        DrawComparisonHist(yan_hists, wei_graph, geneva_graph, config, ek_bins, masses[i]);
		
		yan_file->Close();
		delete yan_file;
    }
}

void DrawCompare() {
    CompareBeResults();
}