#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>
#include <algorithm>
#include <iomanip>
#include "../Tool.h"
using namespace AMS_Iso;

// 配置
const std::vector<std::string> elements = {"Beryllium", "Boron"};//, "Carbon", "Nitrogen", "Oxygen"};
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::vector<std::string> methods = {"LG", "EGE", "EGELinear"};
const std::vector<int> cpcolors = {kRed, kGreen+2, kBlue};
const std::vector<int> markers = {21, 22, 23};

const std::string inputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/";
const std::string outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/";
const std::string l1FileName = "ChargeTemp_Hist_0p8.root";
const std::vector<std::string> l2FileNames = {
    "L2TempTuned_LG_detfuncT.root",
    "L2TempTuned_EGE_detfuncT.root",
    "L2TempTuned_EGELinear_NarrowBins.root"
};


// 各探测器能量范围
const std::map<std::string, std::pair<double, double>> detectorRange = {
    {"TOF", {0.51, 1.55}},
    {"NaF", {0.86, 4.00}},
    {"AGL", {2.7, 16.38}}
};

// 获取某探测器覆盖的bin index
std::vector<int> getDetectorBinIndices(const std::string& det) {
    std::vector<int> binIdx;
    auto range = detectorRange.at(det);
    for (size_t i = 0; i + 1 < Binning::NarrowBins.size(); ++i) {
        double low = Binning::NarrowBins[i], high = Binning::NarrowBins[i+1];
        if (high > range.first && low < range.second) binIdx.push_back(i);
    }
    return binIdx;
}

// 获取1D投影
TH1D* getEnergySlice(TH2D* h2d, double elow, double ehigh) {
    if (!h2d) return nullptr;
    TAxis* yAxis = h2d->GetYaxis();
    int binLow = yAxis->FindBin(elow + 1e-5);
    int binHigh = yAxis->FindBin(ehigh - 1e-5);
    if (binLow < 1) binLow = 1;
    if (binHigh > yAxis->GetNbins()) binHigh = yAxis->GetNbins();
    TString sliceName = TString::Format("%s_slice_%.2f_%.2f", h2d->GetName(), elow, ehigh);
    return h2d->ProjectionX(sliceName, binLow, binHigh);
}

// 归一化直方图
void normalizeHist(TH1D* hist, double low ,double up) {
    if (!hist) return;
    double integral = hist->Integral(hist->FindBin(low), hist->FindBin(up));
    if (integral > 0) hist->Scale(1.0 / hist->GetMaximum());//hist->Scale(1.0 / integral);
}

// 计算Pull，err按 pullValue * (err_l1 / val_l1)
TH1D* calculatePull(TH1D* h_l1, TH1D* h_l2, const std::string& name) {
    if (!h_l1 || !h_l2) return nullptr;
    TH1D* pull = (TH1D*)h_l2->Clone(name.c_str());
    pull->Reset();
    //pull->SetTitle("Pull: (L1 - L2) / #sigma_{L1}");
    for (int i = 1; i <= h_l2->GetNbinsX(); ++i) {
        int ibin = h_l1->FindBin(h_l2->GetBinCenter(i)); 
        if(ibin < 1 || ibin > h_l2->GetNbinsX()) continue;
        double val_l1 = h_l1->GetBinContent(ibin);
        double err_l1 = h_l1->GetBinError(ibin);
        double val_l2 = h_l2->GetBinContent(i);
        if (err_l1 > 0 && val_l1 != 0) {
            double pullValue = (val_l1 - val_l2) / err_l1;
            double pullErr = fabs(pullValue * (err_l1 / val_l1));
            pull->SetBinContent(i, pullValue);
            pull->SetBinError(i, pullErr);
            if(h_l1->GetBinCenter(ibin) > 5.1 && h_l1->GetBinCenter(ibin)<5.2){
                cout<<h_l2->GetName()<<endl;
                cout<<h_l1->GetBinCenter(ibin)<<endl;
                cout<<h_l2->GetBinCenter(i)<<endl;
                cout<<val_l1<<endl;
                cout<<val_l2<<endl;
                cout<<err_l1<<endl;
                cout<<pullValue<<endl;
            }
        } else {
            pull->SetBinContent(i, 0);
            pull->SetBinError(i, 0);
        }
    }
    return pull;
}

// 设置直方图样式
void setHistStyle(TH1D* hist, int color, int marker, double size = 0.8) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(marker);
    hist->SetMarkerSize(size);
    hist->SetLineWidth(2);
}

// 创建图例
TLegend* createLegend1() {
    TLegend* legend = new TLegend(0.65, 0.65, 0.9, 0.88);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    return legend;
}

void compareHistograms() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* l1File = new TFile((inputDir + l1FileName).c_str(), "READ");
    if (!l1File || l1File->IsZombie()) {
        std::cerr << "Error: Cannot open L1 file!" << std::endl;
        return;
    }
    std::vector<TFile*> l2Files;
    for (const auto& fileName : l2FileNames) {
        TFile* f = new TFile((inputDir + fileName).c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: Cannot open L2 file: " << fileName << std::endl;
            return;
        }
        l2Files.push_back(f);
    }

    std::string outputPDF = outputDir + "L1L2_Comparison_detfuncT.pdf";
    TCanvas* dummy = new TCanvas();
    dummy->Print((outputPDF + "[").c_str());
    delete dummy;
    
    int z = 4;
    for (const auto& element : elements) {
        for (const auto& detector : detectors) {
            std::vector<int> binIdxs = getDetectorBinIndices(detector);
            for (int idx : binIdxs) {
                double elow = Binning::NarrowBins[idx];
                double ehigh = Binning::NarrowBins[idx+1];
                // 获取L1直方图
                std::string l1HistName = "L1Temp_" + element + "_" + detector + "_Bkg";
                TH2D* h2d_l1 = (TH2D*)l1File->Get(l1HistName.c_str());
                if (!h2d_l1) {
                    std::cout << "Warning: Cannot find L1 histogram: " << l1HistName << std::endl;
                    continue;
                }
                TH1D* h1d_l1 = getEnergySlice(h2d_l1, elow, ehigh);
                if (!h1d_l1) continue;

                // 获取L2直方图
                std::vector<TH1D*> h1d_l2s; bool allL2Valid = true;
                for (int i = 0; i < methods.size(); ++i) {
                    std::string l2HistName = "L2TempTuned_" + methods[i] + "_" + element + "_" + detector;
                    TH2D* h2d_l2 = (TH2D*)l2Files[i]->Get(l2HistName.c_str());
                    if (!h2d_l2) {
                        std::cout << "Warning: Cannot find L2 histogram: " << l2HistName << std::endl;
                        allL2Valid = false; break;
                    }
                    TH1D* h1d_l2 = getEnergySlice(h2d_l2, elow, ehigh);
                    if (!h1d_l2) { allL2Valid = false; break; }
                    h1d_l2s.push_back(h1d_l2);
                }
                if (!allL2Valid) {
                    delete h1d_l1;
                    for (auto h : h1d_l2s) delete h;
                    continue;
                }
                // 归一化
                h1d_l1->Rebin(2);
                double low = z > 4 ? 2 : 1;
                double up  = z < 0 ? 2 : 0.45;
                normalizeHist(h1d_l1, z-low, z+up);
                for (auto h : h1d_l2s){
                    h->Rebin(2);
                    normalizeHist(h, z-low, z+up);
                }

                TCanvas* canvas = new TCanvas("canvas", "L1 L2 Comparison", 800, 600);
                canvas->Divide(1, 2);

                // 上方pad
                TPad* pad1 = (TPad*)canvas->cd(1);
                //pad1->SetLogy();
                pad1->SetPad(0, 0.27, 1, 1);
                pad1->SetBottomMargin(0.02);
                setHistStyle(h1d_l1, kBlack, 20, 0.8);
                for (int i = 0; i < methods.size(); ++i)
                    setHistStyle(h1d_l2s[i], cpcolors[i], markers[i], 0.8);

                std::ostringstream oss;
                oss << element << " " << detector
                    << std::fixed << std::setprecision(2)
                    << "  (E_{k}/n: " << elow << "-" << ehigh << " GeV/n)";
                h1d_l1->SetTitle(oss.str().c_str());
                h1d_l1->GetXaxis()->SetTitle("Charge");
                h1d_l1->GetXaxis()->SetRangeUser(z-1, z+1);
                h1d_l1->GetYaxis()->SetRangeUser(0.0002, 1.25*h1d_l1->GetMaximum());
                h1d_l1->GetYaxis()->SetTitle("Normalized Counts");
                h1d_l1->GetXaxis()->SetLabelSize(0);

                h1d_l1->Draw("hist");
                pad1->SetLogy(0);
                for (auto h : h1d_l2s) h->Draw("hist SAME");

                TLegend* legend = createLegend1();
                legend->AddEntry(h1d_l1, "L1Temp", "ep");
                legend->AddEntry(h1d_l2s[0], "LG Tuned L2Temp", "l");
                legend->AddEntry(h1d_l2s[1], "EGE Tuned L2Temp", "l");
                legend->AddEntry(h1d_l2s[2], "EGELinear L2Temp", "l");
                legend->Draw();

                // 下方pad - Pull分布
                TPad* pad2 = (TPad*)canvas->cd(2);
                pad2->SetPad(0, 0, 1, 0.25);
                pad2->SetBottomMargin(0.3);
                pad2->SetGridy();
                pad2->SetTopMargin(0.02);

                std::vector<TH1D*> pulls;
                for (int i = 0; i < h1d_l2s.size(); ++i) {
                    std::ostringstream pullName;
                    pullName << "pull_" << methods[i] << "_" << element << "_" << detector << "_" << idx;
                    TH1D* pull = calculatePull(h1d_l1, h1d_l2s[i], pullName.str());
                    if (pull) {
                        setHistStyle(pull, cpcolors[i], markers[i], 0.8);
                        pulls.push_back(pull);
                    }
                }
                if (!pulls.empty()) {
                    pulls[0]->GetXaxis()->SetTitle("Charge");
                    pulls[0]->GetYaxis()->SetTitle("Pull");
                    pulls[0]->GetXaxis()->SetTitleSize(0.15);
                    pulls[0]->GetYaxis()->SetTitleSize(0.15);
                    pulls[0]->GetXaxis()->SetLabelSize(0.15);
                    pulls[0]->GetYaxis()->SetLabelSize(0.15);
                    pulls[0]->GetXaxis()->SetTitleOffset(.9);
                    pulls[0]->GetYaxis()->SetTitleOffset(0.4);
                    pulls[0]->GetYaxis()->SetNdivisions(505);
                    pulls[0]->GetXaxis()->SetRangeUser(z-1, z+1);
                    pulls[0]->SetTitle("");

                    double maxPull = 0;
                    for (auto pull : pulls)
                        for (int i = 1; i <= pull->GetNbinsX(); ++i)
                            maxPull = std::max(maxPull, std::abs(pull->GetBinContent(i)));
                    pulls[0]->GetYaxis()->SetRangeUser(-10, 10);

                    pulls[0]->Draw("E");
                    for (size_t i = 1; i < pulls.size(); ++i)
                        pulls[i]->Draw("E SAME");

                    TLine* zeroLine = new TLine(-5, 0,
                                                5, 0);
                    zeroLine->SetLineStyle(2);
                    zeroLine->SetLineColor(kBlack);
                    zeroLine->Draw("SAME");
                    delete zeroLine;
                }
                canvas->Print(outputPDF.c_str());
                delete h1d_l1;
                for (auto h : h1d_l2s) delete h;
                for (auto h : pulls) delete h;
                delete canvas;

                std::cout << "Processed: " << element << "_" << detector
                          << " [" << elow << ", " << ehigh << "]" << std::endl;
            }
        }
        z=z+1;
    }

    TCanvas* dummy2 = new TCanvas();
    dummy2->Print((outputPDF + "]").c_str());
    delete dummy2;

    l1File->Close();
    for (auto f : l2Files) f->Close();
    delete l1File;
    for (auto f : l2Files) delete f;

    std::cout << "Comparison plots saved to: " << outputPDF << std::endl;
}

void cpL1L2() {
    compareHistograms();
}