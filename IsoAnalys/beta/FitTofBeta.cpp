#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <algorithm>
#include "/afs/cern.ch/work/z/zuhao/public/yanzx/IsoAnalys/Tool.h"

using namespace std;
using namespace AMS_Iso;

const int NBin_ek = Constants::RIGIDITY_BINS;
const char* DecName[2] = {
    Detector::RichDetectorNames[0].c_str(),  // "NaF"
    Detector::RichDetectorNames[1].c_str()   // "AGL"
};

// 获取对应charge和mass的能量bins
const double* GetEnergyBins(int charge, int mass) {
    for (const auto& iso : IsotopeData) {
        if (iso.charge_ == charge) {
            for (size_t i = 0; i < Constants::MAX_ISOTOPES; ++i) {
                if (iso.mass_[i] == mass) {
                    return &Binning::KineticEnergyBins[charge-1][i][0];
                }
            }
        }
    }
    throw runtime_error("Invalid charge or mass");
}

// 查找bin index
int FindBinIndex(double value, const double* bins) {
    for (int i = 0; i < NBin_ek; i++) {
        if (value >= bins[i] && value < bins[i+1]) {
            return i;
        }
    }
    return -1;
}

// 绘制拟合结果并保存到画布
void DrawFitResults(TH1D *h, RooFitResult *fitResult, RooPlot *frame, TCanvas *canvas, 
    const char *name, RooAbsPdf &model, RooDataHist &data) {
    canvas->SetLogy();
    h->Draw("e");
    h->SetMarkerSize(1.);
    h->SetMarkerStyle(20);
    frame->Draw("same");

    TPaveText *pt = new TPaveText(0.65, 0.68, 0.76, .88, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.035);

    RooArgList params = fitResult->floatParsFinal();
    for (int i = 0; i < params.getSize(); ++i) {
        RooRealVar *param = (RooRealVar *)params.at(i);
        pt->AddText(Form("%s = %.5f +/- %.5f", param->GetName(), param->getVal(), param->getError()));
    }

    double chi2_value = frame->chiSquare();
    double numFitParameters = model.getParameters(data)->getSize();
    double numDataBins = data.numEntries();
    double degreesOfFreedom = numDataBins - numFitParameters;
    if(degreesOfFreedom>0)
        pt->AddText(Form("Chi2/NDF = %f", chi2_value / degreesOfFreedom));
    pt->Draw("same");
    canvas->Print(name);
}

// 填充直方图
void FillHistograms(int i, TH1D* hsig, TH1D* hsiga, TH1D* hMean, TH1D* hMeana, 
    double SigVal[], double SigErr[], double SigVala[], double SigErra[], 
    double MeanVal[], double MeanErr[], double MeanVala[], double MeanErra[]) {
    
    hsig->SetBinContent(i, SigVal[i]);
    hsig->SetBinError(i, SigErr[i]);
    hMean->SetBinContent(i, MeanVal[i]);
    hMean->SetBinError(i, MeanErr[i]);
    if (i > 15) {
        hsiga->SetBinContent(i, SigVala[i]);
        hsiga->SetBinError(i, SigErra[i]);
        hMeana->SetBinContent(i, MeanVala[i]);
        hMeana->SetBinError(i, MeanErra[i]);
    }
}

// 设置直方图样式
void SetHistogramStyle(TH1D* hsig, TH1D* hsiga, TH1D* hMean, TH1D* hMeana) {
    hsig->GetXaxis()->SetRangeUser(0.49, 14);
    hsig->SetMarkerColor(kRed);
    hsiga->SetMarkerColor(kBlue);
    hsig->SetMarkerStyle(20);
    hsiga->SetMarkerStyle(20);
    hsig->SetMarkerSize(1.5);
    hsiga->SetMarkerSize(1.5);
    hsig->SetLineColor(kRed);
    hsiga->SetLineColor(kBlue);
    hsig->GetYaxis()->SetLabelSize(0.04);
    hsig->GetYaxis()->SetTitleSize(0.05);
    hsig->GetYaxis()->SetTitleOffset(1.3);

    hMean->GetXaxis()->SetRangeUser(0.49, 14);
    hMean->SetMarkerColor(kRed);
    hMeana->SetMarkerColor(kBlue);
    hMean->SetMarkerStyle(20);
    hMeana->SetMarkerStyle(20);
    hMean->SetMarkerSize(1.5);
    hMeana->SetMarkerSize(1.5);
    hMean->SetLineColor(kRed);
    hMeana->SetLineColor(kBlue);
    hMean->GetYaxis()->SetLabelSize(0.04);
    hMean->GetYaxis()->SetTitleSize(0.05);
    hMean->GetYaxis()->SetTitleOffset(1.4);
}
void DrawHistograms(TCanvas* canvas, TH1D* hsig, TH1D* hsiga, TH1D* hMean, TH1D* hMeana, 
    TString pdfFilename) {
    
    canvas->cd();
    canvas->SetLogy(0);
    canvas->SetLogx();
    TLegend *legend = new TLegend(0.65, 0.15, 0.85, 0.45);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    hsig->Draw("e");
    hsiga->Draw("same");
    legend->AddEntry(hsig, "1/beta_{TOF}-1/beta_{NaF}", "ep");
    legend->AddEntry(hsiga, "1/beta_{TOF}-1/beta_{AGL}", "ep");
    legend->Draw("same");
    canvas->Print(pdfFilename);
    canvas->Clear();
    
    legend = new TLegend(0.65, 0.15, 0.85, 0.45);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    hMean->Draw("e");
    hMeana->Draw("same");
    legend->AddEntry(hMean, "1/beta_{TOF}-1/beta_{NaF}", "ep");
    legend->AddEntry(hMeana, "1/beta_{TOF}-1/beta_{AGL}", "ep");
    legend->Draw("same");
    canvas->Print(pdfFilename);
    canvas->SetLogx(0);
}

void FitBeta(const char *input, const char *output, const char *histname, int charge, int mass) {
    const double* energy_bins = GetEnergyBins(charge, mass);
    int startBin = FindBinIndex(0.3, energy_bins);
    int endBin = FindBinIndex(14.65, energy_bins);

    TFile *file = TFile::Open(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/%s.root", input));
    TString pdfFilename = Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/FitTOFinEkBin_%s.pdf", output);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1600, 1200);

    canvas->Print(pdfFilename + "[");

    double SigVal[NBin_ek], SigErr[NBin_ek], BinEdg[NBin_ek];
    double SigVala[NBin_ek], SigErra[NBin_ek], BinEdga[NBin_ek];
    double MeanVal[NBin_ek], MeanErr[NBin_ek];
    double MeanVala[NBin_ek], MeanErra[NBin_ek];

    for (int dec = 0; dec < 2; dec++) {
        cout << "Processing " << DecName[dec] << " detector" << endl;
        
        TH2D *h2 = (TH2D*)file->Get(Form("RBeta_Ek_TOF_%s_%s", DecName[dec], histname));
        if (!h2) {
            cout << "Error: histogram not found for " << DecName[dec] << endl;
            continue;
        }

        for (int i = startBin; i <= endBin; i++) {
            TH1D *h = h2->ProjectionX("h", i, i);
            h->Rebin(100);
            if (h->GetMaximum() < 10) { 
                delete h;
                continue; 
            }

            h->SetYTitle("events");
            h->SetTitle(Form("%s %.2f GeV/n~%.2f GeV/n", output, energy_bins[i-1], energy_bins[i]));

            gStyle->SetErrorX(0);
            h->SetStats(0);
            h->GetXaxis()->SetLabelSize(0.045);
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetLabelSize(0.05);
            h->GetYaxis()->SetTitleSize(0.06);
            h->GetXaxis()->SetLabelOffset(.006);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.3);

            TF1 *gaussFit = new TF1("gaussFit", "gaus", -0.05, 0.05);
            h->Fit(gaussFit, "R", "", -0.06, 0.06);
            h->GetXaxis()->SetRangeUser(-0.06, 0.06);

            double m = gaussFit->GetParameter(1);
            double sig = gaussFit->GetParameter(2);

            RooRealVar x("x", "Observable", m - 3 * sig, m + 3 * sig);
            RooDataHist data("data", "Dataset", x, h);
            RooRealVar mean("mean", "Mean", m, m - 2*sig, m + 2*sig);
            RooRealVar sigma_core("sigma_core", "Core sigma", sig, 0.5 * sig, 5 * sig);
            RooGaussian core("core", "Core Gaussian", x, mean, sigma_core);
            RooFitResult *fitResult = core.fitTo(data, RooFit::Save());
            RooPlot *frame = x.frame();
            data.plotOn(frame, RooFit::XErrorSize(0));
            core.plotOn(frame, RooFit::LineColor(kRed));

            DrawFitResults(h, fitResult, frame, canvas, pdfFilename, core, data);

            if (dec == 0) {
                SigVal[i] = sigma_core.getVal();
                SigErr[i] = sigma_core.getError();
                MeanVal[i] = mean.getVal();
                MeanErr[i] = mean.getError();
            } else {
                SigVala[i] = sigma_core.getVal();
                SigErra[i] = sigma_core.getError();
                MeanVala[i] = mean.getVal();
                MeanErra[i] = mean.getError();
            }

            delete h;
            delete gaussFit;
            delete fitResult;
            delete frame;
        }
    }

    TH1D *hsig = new TH1D(Form("TOF_NaF_Sig_%s", output), "Sigma;Ek/n[GeV/n];Sigma of 1/beta_{TOF}-1/beta_{RICH}", 
                          NBin_ek, energy_bins);
    TH1D *hsiga = new TH1D(Form("TOF_AGL_Sig_%s", output), "Sigma;Ek/n[GeV/n];Sigma of 1/beta_{TOF}-1/beta_{RICH}", 
                           NBin_ek, energy_bins);
    TH1D *hMean = new TH1D(Form("TOF_NaF_Mean_%s", output), "Mean;Ek/n[GeV/n];Mean of 1/beta_{TOF}-1/beta_{RICH}", 
                           NBin_ek, energy_bins);
    TH1D *hMeana = new TH1D(Form("TOF_AGL_Mean_%s", output), "Mean;Ek/n[GeV/n];Mean of 1/beta_{TOF}-1/beta_{RICH}", 
                            NBin_ek, energy_bins);

    for (int i = startBin; i <= endBin; i++) {
        FillHistograms(i, hsig, hsiga, hMean, hMeana, SigVal, SigErr, SigVala, SigErra, 
                      MeanVal, MeanErr, MeanVala, MeanErra);
    }

    SetHistogramStyle(hsig, hsiga, hMean, hMeana);
    DrawHistograms(canvas, hsig, hsiga, hMean, hMeana, pdfFilename);
    canvas->Print(pdfFilename + "]");

    TFile *fout = new TFile(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/FitTOFinEkBin_%s.root", 
                                output), "recreate");
    hsig->Write();    
    hsiga->Write();    
    hMean->Write();
    hMeana->Write();
    fout->Close();

    delete hsig;
    delete hsiga;
    delete hMean;
    delete hMeana;
    delete canvas;
    file->Close();
}

void FitTofBeta() {
    gStyle->SetErrorX(0);
    // Beryllium
    FitBeta("Be7", "Be7", "Berlium", 4, 7);
    FitBeta("Be9", "Be9", "Berlium", 4, 9);
    FitBeta("Be10", "Be10", "Berlium", 4, 10);
    // ISS Beryllium
    FitBeta("Ber7", "Ber7", "Berlium", 4, 7);
    FitBeta("Ber9", "Ber9", "Berlium", 4, 9);
    FitBeta("Ber10", "Ber10", "Berlium", 4, 10);
}

#ifndef __CLING__
int main() {
    FitTofBeta();
    return 0;
}
#endif