#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

using namespace std;

struct FitConfig {
    std::string variable_name;
    double fit_min;
    double fit_max;
    double sig_min;
    double sig_max;
    double l_xmin;
    double l_xmax;
    int hist_suffix;
    std::string output_prefix;
};

double betaToEkn(double beta) {
    const double mass = 0.931;
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    return (gamma - 1.0) * mass;
}

double eknToBeta(double ekn) {
    const double mass = 0.931;
    double gamma = ekn/mass + 1;
    return sqrt(1 - 1/(gamma*gamma));
}

std::pair<int, int> getBetaBinRange(TH2D* hist, double ekn_low, double ekn_high) {
    double beta_low = eknToBeta(ekn_low);
    double beta_high = eknToBeta(ekn_high);
    cout<<" beta_low:"<<beta_low<<" beta_high:"<<beta_high<<endl;
    
    int bin_low = hist->GetXaxis()->FindBin(beta_low);
    int bin_high = hist->GetXaxis()->FindBin(beta_high);
    
    return std::make_pair(bin_low, bin_high-1);
}

bool isInDetectorRange(const std::string& detector, double ekn) {
    if (detector == "TOF") return (ekn >= 0.3 && ekn <= 2.5);
    if (detector == "NaF") return (ekn >= 0.4 && ekn <= 6.5);
    if (detector == "Agl") return (ekn >= 2. && ekn <= 15);
    return false;
}

const std::vector<double> Be_bins = {
    0.08, 0.13, 0.17, 0.21, 0.27, 0.33, 0.41, 0.49, 0.59, 0.70, 0.82, 0.96, 1.11, 
    1.28, 1.47, 1.68, 1.91, 2.16, 2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 
    5.64, 6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51, 11.45, 12.45, 13.50, 14.65
};

const std::vector<std::string> detectors = {"TOF", "NaF", "Agl"};
const std::map<std::string, int> detectorTypeMap = {
    {"TOF", 0},
    {"NaF", 1},
    {"Agl", 2}
};

void setupPullPlot(TGraphErrors* pullGraph, const std::string& variable_name, 
                  double fit_min, double fit_max) {
    pullGraph->SetTitle(Form(";%s;Pull", variable_name.c_str()));
    pullGraph->SetMarkerColor(kBlack);
    pullGraph->SetMarkerStyle(20);
    pullGraph->SetMarkerSize(0.8);

    pullGraph->GetXaxis()->SetLabelSize(0.12);
    pullGraph->GetYaxis()->SetLabelSize(0.09);
    pullGraph->GetYaxis()->SetNdivisions(010);
    pullGraph->GetXaxis()->SetTitleSize(0.15);
    pullGraph->GetYaxis()->SetTitleSize(0.15);
    pullGraph->GetYaxis()->SetTitleOffset(0.4);
    pullGraph->GetXaxis()->SetTitleOffset(0.7);
    pullGraph->GetXaxis()->SetRangeUser(fit_min, fit_max);
}

void calculatePull(RooPlot* frame, TGraphErrors* pullGraph, const char *DataName, const char* ModelName,  double fit_min, double fit_max) {
    RooHist* dataHist_frame = frame->getHist(Form("%s", DataName));
    RooCurve* modelCurve = frame->getCurve(Form("%s", ModelName));

    if (dataHist_frame && modelCurve) {
        double x_val, y_data;
        for (int i = 0; i < dataHist_frame->GetN(); ++i) {
            dataHist_frame->GetPoint(i, x_val, y_data);
            
            if (x_val >= fit_min && x_val <= fit_max) {
                double y_model = modelCurve->Eval(x_val);
                double data_err = dataHist_frame->GetErrorY(i);
                double pull_err = 0;
                
                if (y_model > 0 && y_data > 0) {
                    double pull_val = (y_data - y_model) / sqrt(y_data);
                    pull_err = pull_val * (data_err / y_data);
                    if (pull_val > -10 && pull_val < 10 && std::isfinite(pull_val) && std::isfinite(pull_err)) {
                        cout<<"1/mass: "<<x_val<<" pull:"<<pull_val<<endl;
                        pullGraph->SetPoint(i, x_val, pull_val);
                        pullGraph->SetPointError(i, 0, pull_err);  
                    }
                }
            }
        }
    }
}

TLegend* createLegend(RooPlot* frame) {
    TLegend* legend = new TLegend(0.2, 0.68, 0.45, 0.87);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->AddEntry(frame->findObject("Data"), "ISS Data", "p");
    legend->AddEntry(frame->findObject("Be7"), "Be7", "l");
    legend->AddEntry(frame->findObject("Be9"), "Be9", "l");
    legend->AddEntry(frame->findObject("Be10"), "Be10", "l");
    return legend;
}

TPaveText* createFitInfo(double l_xmin, double l_xmax, 
                        double fracBe7, double fracBe7Err,
                        double fracBe9, double fracBe9Err,
                        double chi2_ndf, int fitStatus,
                        double sig_min, double sig_max, double bestEfficiency,
                        double significance,
                        double N_sig, double N_tot) {
    TPaveText* pt = new TPaveText(l_xmin, 0.66, l_xmax, 0.88, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.035);
    if (fitStatus != 0) {
        pt->AddText(Form("Wrong Fit Status(!= 0)"));
        return pt;
    }
    pt->AddText(Form("Be7 Frac = %.3f #pm %.3f", fracBe7, fracBe7Err));
    pt->AddText(Form("Be9 Frac = %.3f #pm %.3f", fracBe9, fracBe9Err));
    pt->AddText(Form("Chi2/NDF = %.3f", chi2_ndf));
    pt->AddText(Form("%.2f Signal(Be10) Efficiency:", bestEfficiency));
    pt->AddText(Form("Significance(%.2f-%.2f): %.1f", sig_min, sig_max, significance));
    pt->AddText(Form("N_{sig}: %.1f", N_sig));
    pt->AddText(Form("N_{tot}: %.1f", N_tot));
    return pt;
}

double findUpperBound(RooRealVar& x, RooHistPdf* signalPdf, double lowerBound, double initialUpperBound, double targetFraction = 0.6) {
    if (!signalPdf) {
        std::cout << "Error: Null pointer passed to findUpperBound" << std::endl;
        return initialUpperBound;
    }
    
    double upperBound = initialUpperBound;
    double lowerSearch = lowerBound;
    double upperSearch = initialUpperBound;
    double tolerance = 1e-4;
    
    x.setRange("fullRange", lowerBound, initialUpperBound);
    RooAbsReal* full_integral = signalPdf->createIntegral(RooArgSet(x), 
                                                         RooFit::NormSet(RooArgSet(x)), 
                                                         RooFit::Range("fullRange"));
    double fullIntegral = full_integral->getVal();
    
    int iterations = 0;
    while (upperSearch - lowerSearch > tolerance && iterations < 100) {
        upperBound = (upperSearch + lowerSearch) / 2.0;
        
        x.setRange("testRange", lowerBound, upperBound);
        RooAbsReal* test_integral = signalPdf->createIntegral(RooArgSet(x), 
                                                             RooFit::NormSet(RooArgSet(x)), 
                                                             RooFit::Range("testRange"));
        double currentFraction = test_integral->getVal() / fullIntegral;
        
        if (currentFraction < targetFraction) {
            lowerSearch = upperBound;
        } else {
            upperSearch = upperBound;
        }
        
        delete test_integral;
        iterations++;
    }
    
    delete full_integral;
    return upperBound;
}

struct SignificanceResult {
    double efficiency;      // 最优信号效率
    double upperBound;      // 对应的积分上界
    double N_sig;          // 信号事例数
    double N_tot;          // 总事例数
    double significance;    // 显著度
};

SignificanceResult findOptimalSignificance(RooRealVar& x, RooHistPdf* signalPdf, RooAbsData& dataHist, 
                                         RooAddPdf& totalPdf, RooRealVar& fracBe7, RooRealVar& fracBe9,
                                         double fitLow, double initialUpper, 
                                         double minEfficiency = 0.5, double maxEfficiency = 0.99, 
                                         double efficiencyStep = 0.01) {
    SignificanceResult bestResult = {0, 0, 0, 0, 0};
    double maxSignificance = 0;
    
    // 遍历不同的信号效率
    for (double eff = minEfficiency; eff <= maxEfficiency; eff += efficiencyStep) {
        // 找到当前效率对应的上界
        double currentUpper = findUpperBound(x, signalPdf, fitLow, initialUpper, eff);
        
        // 计算总事例数
        x.setRange("sigRange", fitLow, currentUpper);
        std::unique_ptr<RooAbsReal> total_integral(
            totalPdf.createIntegral(RooArgSet(x), RooFit::NormSet(RooArgSet(x)), RooFit::Range("sigRange"))
        );
        double N_tot = total_integral->getVal() * dataHist.sumEntries();
        
        // 计算信号事例数
        std::unique_ptr<RooAbsReal> sig_integral(
            signalPdf->createIntegral(RooArgSet(x), RooFit::NormSet(RooArgSet(x)), RooFit::Range("sigRange"))
        );
        double fracBe10 = 1.0 - fracBe7.getVal() - fracBe9.getVal();
        double N_sig = sig_integral->getVal() * fracBe10 * dataHist.sumEntries();
        
        // 计算显著度
        double significance = N_sig / sqrt(N_tot);
        
        // 更新最优结果
        if (significance > maxSignificance) {
            maxSignificance = significance;
            bestResult.efficiency = eff;
            bestResult.upperBound = currentUpper;
            bestResult.N_sig = N_sig;
            bestResult.N_tot = N_tot;
            bestResult.significance = significance;
        }
    }
    
    return bestResult;
}