#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TKey.h>
#include <TClass.h>
#include <TSystem.h>
#include <iostream>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooBifurGauss.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <vector>
#include <string>

using namespace std;

struct IsotopeConfig {
    string name;          // 元素名称
    vector<int> masses;   // 质量数
    double nafRange;      // NaF拟合范围
    double aglRange;      // AGL拟合范围
    string histname;   
};

double GetNGen(TH1D *hMCnum) {
    return (hMCnum) ? hMCnum->GetBinContent(1) : 1.0;
}

TH1* FilterHistogram(TH1* h, double minCounts) {
    TH1* hFiltered = (TH1*)h->Clone();
    hFiltered->Reset();

    int nBins = h->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double binContent = h->GetBinContent(i);
        double binError = h->GetBinError(i);
        if (binContent >= minCounts) {
            hFiltered->SetBinContent(i, binContent);
            hFiltered->SetBinError(i, binError);
        }
    }
    return hFiltered;
}

const char* DecName[2] = {"NaF", "AGL"};

void ProjectAndDraw(const char* inputFile, const char* outputName, 
                   double nafRange, double aglRange, const char* histname, bool isISS) {
    cout << "Processing file: " << inputFile << endl;
    cout << "Output name: " << outputName << endl;
    cout << "Detector type: " << (isISS ? "ISS" : "MC") << endl;
    
    double r[2] = {nafRange, aglRange};
    double Bmin[2] = {1 - nafRange, 1 - aglRange};
    double Bmax[2] = {1 + nafRange, 1 + aglRange};
    double CoreSmin[2] = {0.0005, 0.0003};
    double CoreSmax[2] = {0.002, 0.001};
    double frcmin[2] = {0.7, 0.7};
    double frcmax[2] = {.9, .95};
    double LRmin[2] = {1.1, 1.1};
    double LRmax[2] = {2.2, 2.4};
    double RRmin[2] = {1.1, 1.1};
    double RRmax[2] = {2.2, 2.4};

    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        cout << "Error opening file: " << inputFile << endl;
        return;
    }

    string outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/richBeta/";
    if (!gSystem->AccessPathName(outputDir.c_str())) {
        gSystem->mkdir(outputDir.c_str(), true);
    }
    
    TString shortName(outputName);
    TString pdfFilename = Form("%s/FitRich%s.pdf", outputDir.c_str(), outputName);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);

    gPad->SetLogy();
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.15);
    canvas->Print(pdfFilename + "[");

    for (int n = 0; n < 2; n++) {
        cout << "Processing " << DecName[n] << " detector" << endl;
        
        TH1D *h = (TH1D*)file->Get(Form("R_%sBeta_%s", DecName[n], histname));
        if (!h) {
            cout << "Error: histogram not found for " << DecName[n] << endl;
            continue;
        }
        
        h->Rebin(isISS ? 50 : 25);
        h->Sumw2();
        h->SetStats(0);
        h->SetTitle(Form("%s InnerRig>%dGV", isISS ? TString(outputName).Remove(3).Data() : outputName, 100 * (n + 1)));

        TH1 *hFiltered = FilterHistogram(h, std::max(0.,h->GetMaximum()/500.));
        
        double xmin = Bmin[n], xmax = Bmax[n];
        hFiltered->GetXaxis()->SetRangeUser(xmin, xmax);
        hFiltered->GetXaxis()->SetLabelSize(0.055);
        hFiltered->GetYaxis()->SetLabelSize(0.055);
        hFiltered->GetXaxis()->SetTitleSize(0.055);
        hFiltered->GetYaxis()->SetTitleSize(0.055);
        RooRealVar x("x", "Observable", xmin, xmax);
        RooDataHist data("data", "Dataset", x, hFiltered);

        RooRealVar mean("mean", "Mean", hFiltered->GetMean(), 
                       hFiltered->GetMean() - 0.001, hFiltered->GetMean() + 0.001);
        RooRealVar sigma_core("sigma_core", "Core sigma", 
                             hFiltered->GetRMS(), CoreSmin[n], CoreSmax[n]);
        RooGaussian core("core", "Core Gaussian", x, mean, sigma_core);

        RooRealVar LR("LR", "Left sigma ratio", 1.8, LRmin[n], LRmax[n]);
        RooRealVar RR("RR", "Right sigma ratio", 1.8, RRmin[n], RRmax[n]);
        RooFormulaVar sigma_L("sigma_L", "Left sigma", "@0*@1", RooArgList(LR, sigma_core));
        RooFormulaVar sigma_R("sigma_R", "Right sigma", "@0*@1", RooArgList(RR, sigma_core));
        RooBifurGauss asy_gaus("asy_gaus", "Asymmetric Gaussian", x, mean, sigma_L, sigma_R);

        RooRealVar frac("frac", "Fraction", 0.8);
        frac.setConstant(true);
        RooAddPdf model("model", "Gaussian + Asymmetric Gaussian", 
                       RooArgList(core, asy_gaus), frac);

        RooFitResult *fitResult = model.fitTo(data, RooFit::Save());
        RooPlot *frame = x.frame();
        data.plotOn(frame, RooFit::XErrorSize(0));
        model.plotOn(frame, RooFit::LineColor(kRed));
        model.plotOn(frame, RooFit::LineColor(kBlue), 
                    RooFit::LineStyle(kDashed), RooFit::Components(core));
        model.plotOn(frame, RooFit::LineColor(kGreen), 
                    RooFit::LineStyle(kDashed), RooFit::Components(asy_gaus));

        double chi2_value = frame->chiSquare();
        double numFitParameters = model.getParameters(data)->getSize();
        double numDataBins = data.numEntries();
        double degreesOfFreedom = numDataBins - numFitParameters;
        double chi2_ndf = chi2_value/degreesOfFreedom;

        hFiltered->Draw("e");
        hFiltered->GetXaxis()->SetNdivisions(505);
        frame->Draw("same");

        TPaveText *pt = new TPaveText(0.57, 0.66, .9, 0.88, "NDC");
        //TPaveText *pt = new TPaveText(0.22, 0.2, .73, 0.35, "NDC");
        pt->SetTextColor(kRed);
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);
        pt->SetFillColor(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.03);
      
        RooArgList params = fitResult->floatParsFinal();
        for (int i = 0; i < params.getSize(); ++i) {
            RooRealVar *param = (RooRealVar *)params.at(i);
            const char* name = param->GetName(); // Get parameter name
            double val = param->getVal();       // Get parameter value
            double err = param->getError();     // Get parameter error
        
            if (i < 2) { // 对于前两个参数 (index 0 and 1)
                pt->AddText(Form("%s = %.3f#pm%.3f", name, val, err));
            } else if (i == 2 || i == 3) { // 对于第三和第四个参数 (index 2 and 3)
                // 计算乘以 1e6 后的值和误差
                double err_scaled = err * 1e6;
                pt->AddText(Form("%s = %.6f#pm(%.0f)e-06",
                               name, val, err_scaled));
            } else {
                pt->AddText(Form("%s = %.6f#pm%.6f", name, val, err));
            }
        }
        
        // 这部分保持不变
        pt->AddText("fix core frac in 0.8");
        pt->AddText(Form("Chi2/NDF: %.2f", chi2_ndf));
        pt->Draw("same");
        
        canvas->Update();
        canvas->SaveAs(Form("%s/FitRich%s_%d.png", outputDir.c_str(), outputName, n));
        canvas->Print(pdfFilename);
    }
    
    canvas->Print(pdfFilename + "]");
    file->Close();
    delete canvas;
}

class RichBetaAnalyzer {
private:
    string inputDir_;
    string outputDir_;
    vector<IsotopeConfig> configs_;
    bool isISS_;

    string GetElementPrefix(const string& element, bool isISS) {
        if (isISS) {
            if (element == "Li") return "Lit";
            if (element == "Be") return "Ber";
            if (element == "B") return "Bor";
        } else {
            return element;
        }
        return "";
    }

public:
    RichBetaAnalyzer(const string& inputDir, const string& outputDir, bool isISS = false) 
        : inputDir_(inputDir)
        , outputDir_(outputDir)
        , isISS_(isISS) {
        configs_ = {
            {"Be", {7, 9, 10}, 0.006, 0.0025, "Beryllium"}
        };
    }

    void AnalyzeAll() {
        for (const auto& config : configs_) {
            const auto& masses = isISS_ ? vector<int>{config.masses.front()} : config.masses;
            for (int mass : masses) {
                string prefix = GetElementPrefix(config.name, isISS_);
                string inputFile = inputDir_ + prefix + to_string(mass) + ".root";
                string outputName = prefix + to_string(mass);
                
                cout << "\nProcessing mass " << mass << " for " << config.name << endl;
                ProjectAndDraw(inputFile.c_str(), outputName.c_str(), 
                             config.nafRange, config.aglRange, 
                             config.histname.c_str(), isISS_);
            }
        }
    }
};

void FitRichBeta() {
    string inputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/20250418NewData/";
    string outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/richBeta/";
    
    cout << "Starting analysis..." << endl;
    
    // 处理MC数据
    cout << "\nProcessing MC data..." << endl;
    RichBetaAnalyzer mcAnalyzer(inputDir, outputDir, false);
    mcAnalyzer.AnalyzeAll();
    
    // 处理ISS数据
    cout << "\nProcessing ISS data..." << endl;
    RichBetaAnalyzer issAnalyzer(inputDir, outputDir, true);
    issAnalyzer.AnalyzeAll();
    
    cout << "\nAnalysis completed." << endl;
}

#ifndef __CLING__
int main() {
    FitRichBeta();
    return 0;
}
#endif