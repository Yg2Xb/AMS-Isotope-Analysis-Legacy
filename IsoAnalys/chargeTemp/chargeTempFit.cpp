#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <array>
#include <map>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TPaveText.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooFormulaVar.h> 
#include "../Tool.h" 

using namespace RooFit;
using namespace AMS_Iso;

// --- Configuration ---
const std::string l2FileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/L2TempTuned_EGE_detfuncT.root";
const std::string l1FileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p8.root";
const std::string outputPDF = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_detfuncT_unb_wide.pdf";
const std::string outputROOT = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_detfuncT_unb_wide.root";

const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::vector<std::string> elements = {"Beryllium", "Boron", "Carbon"};

const int N_TOYS = 1;

struct DetectorRange { int start_bin; int end_bin; double ek_min; double ek_max; };
const std::map<std::string, DetectorRange> detectorRanges = {
    {"TOF", {1, 6, 0.51, 1.17}},
    {"NaF", {5, 14, 1.17, 3.23}},
    {"AGL", {13, 30, 3.23, 16.38}}
};

struct FractionResult { double be_frac = 0, b_frac = 0, c_frac = 0; };

// --- Helper Function: Create Info Box ---
std::unique_ptr<TPaveText> createInfoBox(
    double chi2_ndf, int ndf, int fitStatus,
    const RooRealVar& frac_Be, const RooRealVar& frac_B,
    const FractionResult& narrowFractions,
    double Nsignal
) {
    auto info = std::make_unique<TPaveText>(0.2, 0.68, 0.45, 0.87, "NDC");
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextSize(0.035);
    info->SetTextAlign(12);
   
    if (fitStatus == 0) {
        info->AddText(Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2_ndf*ndf, ndf, chi2_ndf));
        info->AddText(Form("L1Q 3.5-6.3, Be: %.4f #pm%.4f, B: %.4f #pm%.4f", frac_Be.getVal(), frac_Be.getError(), frac_B.getVal(), frac_B.getError()));
        if (Nsignal > 0) {
             info->AddText(Form("L1Q 5.0-5.4, Be: %.4f #pm%.4f, B: %.4f #pm%.4f", 
                         narrowFractions.be_frac, sqrt(narrowFractions.be_frac*(1-narrowFractions.be_frac)/Nsignal),
                         narrowFractions.b_frac, sqrt(narrowFractions.b_frac*(1-narrowFractions.b_frac)/Nsignal)));
        }
    } else {
        info->AddText("Fit Failed");
    }
    return info;
}


// --- Fit Processor Class ---
class FitProcessor {
public:
    FitProcessor(int bin_idx, const std::string& det, TFile* l1File, TFile* l2File);
    ~FitProcessor() = default;

    bool initialize();
    RooPlot* runPrimaryFit();
    void runToySimulations(int n_toys, const std::map<std::string, TH2D*>& toy_hists);
    
    FractionResult getPrimaryNarrowFractions() const { return primaryNarrowFractions; }
    RooFitResult* getFitResult() const { return fitResult.get(); }
    double getChi2NDF() const { return chi2_ndf; }
    int getNDF() const { return ndf; }
    double getSignalEventsInRange() const { return Nsignal_range; }
    TH1D* getSignalHist() const { return h_signal_orig.get(); }
    RooRealVar& getChargeVar() { return charge; }
    RooRealVar& getFracBeVar() { return frac_Be; }
    RooRealVar& getFracBVar() { return frac_B; }

private:
    int binIndex;
    std::string detectorName;
    TFile *l1FilePtr, *l2FilePtr;

    std::unique_ptr<TH1D> h_signal_orig;
    std::vector<std::unique_ptr<TH1D>> templates_hist_orig;

    RooRealVar charge;
    RooRealVar frac_Be, frac_B;
    std::unique_ptr<RooFormulaVar> frac_C; 
    
    std::vector<std::unique_ptr<RooDataHist>> template_data_hists; 
    std::vector<std::unique_ptr<RooHistPdf>> template_pdfs;
    
    std::unique_ptr<RooAddPdf> total_pdf;
    std::unique_ptr<RooDataHist> data_hist;
    std::unique_ptr<RooFitResult> fitResult;

    double chi2_ndf = 0, Nsignal_range = 0;
    int ndf = 1;
    FractionResult primaryNarrowFractions;

    FractionResult calculateNarrowRangeFractions(RooAbsPdf& model, double frac_be_val, double frac_b_val);
};

// --- Main Analysis Function ---
void fitBoronCharge() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    auto l2File = std::unique_ptr<TFile>(TFile::Open(l2FileName.c_str()));
    auto l1File = std::unique_ptr<TFile>(TFile::Open(l1FileName.c_str()));
    if (!l2File || l2File->IsZombie() || !l1File || l1File->IsZombie()) {
        std::cerr << "Error: Cannot open input files" << std::endl;
        return;
    }
    
    int nbins = Binning::NarrowBins.size() - 1;
    auto h_be_fraction = std::make_unique<TH1D>("h_be_fraction", "Be Fraction (5.0-5.4);E_{k} [GeV/n];Fraction", nbins, Binning::NarrowBins.data());
    auto h_b_fraction = std::make_unique<TH1D>("h_b_fraction", "B Fraction (5.0-5.4);E_{k} [GeV/n];Fraction", nbins, Binning::NarrowBins.data());
    auto h_c_fraction = std::make_unique<TH1D>("h_c_fraction", "C Fraction (5.0-5.4);E_{k} [GeV/n];Fraction", nbins, Binning::NarrowBins.data());
    auto h_signal_events = std::make_unique<TH1D>("h_signal_events", "Signal Events in [5.0, 5.4];E_{k} [GeV/n];Events", nbins, Binning::NarrowBins.data());

    auto h_frac_be_val = std::make_unique<TH1D>("h_frac_be_val", "Fit Value of frac_Be;E_{k} [GeV/n];Value", nbins, Binning::NarrowBins.data());
    auto h_frac_be_err = std::make_unique<TH1D>("h_frac_be_err", "Fit Error of frac_Be;E_{k} [GeV/n];Error", nbins, Binning::NarrowBins.data());
    auto h_frac_b_val = std::make_unique<TH1D>("h_frac_b_val", "Fit Value of frac_B;E_{k} [GeV/n];Value", nbins, Binning::NarrowBins.data());
    auto h_frac_b_err = std::make_unique<TH1D>("h_frac_b_err", "Fit Error of frac_B;E_{k} [GeV/n];Error", nbins, Binning::NarrowBins.data());

    std::map<std::string, TH2D*> toy_hists;
    toy_hists["frac_be"] = new TH2D("h2_toy_frac_be", "Toy MC: Be Fraction;Be Fraction;E_{k} [GeV/n]", 10000, 0., 1.0, nbins, Binning::NarrowBins.data());
    toy_hists["frac_b"] = new TH2D("h2_toy_frac_b", "Toy MC: B Fraction;B Fraction;E_{k} [GeV/n]", 10000, 0., 1.0, nbins, Binning::NarrowBins.data());
    toy_hists["narrow_be"] = new TH2D("h2_toy_narrow_frac_be", "Toy MC: B Fraction (5.0-5.4);Be Fraction;E_{k} [GeV/n]", 2000, 0, 0.0015, nbins, Binning::NarrowBins.data());
    toy_hists["narrow_b"] = new TH2D("h2_toy_narrow_frac_b", "Toy MC: B Fraction (5.0-5.4);B Fraction;E_{k} [GeV/n]", 2000, 0.999, 1, nbins, Binning::NarrowBins.data());

    TCanvas* c0 = new TCanvas("c0", "pdf_canvas");
    c0->Print((outputPDF + "[").c_str());
    int processedPlots = 0;
    
    for (const auto& detector : detectors) {
        auto it = detectorRanges.find(detector);
        if (it == detectorRanges.end()) continue;
        const auto& range = it->second;
        
        for (int binIndex = range.start_bin+1; binIndex < range.end_bin; binIndex = binIndex + 2) {//here
            double eMin = Binning::NarrowBins[binIndex];
            double eMax = Binning::NarrowBins[binIndex + 2];//here
            if (eMax <= range.ek_min || eMin >= range.ek_max) continue;
            
            std::cout << "\nProcessing " << detector << " bin " << binIndex << ": " << eMin << "-" << eMax << " GeV/n" << std::endl;
            
            auto processor = std::make_unique<FitProcessor>(binIndex, detector, l1File.get(), l2File.get());
            if (!processor->initialize()) {
                std::cerr << "  -> Initialization failed, skipping bin." << std::endl;
                continue;
            }

            RooPlot* frame = processor->runPrimaryFit();
            if (!frame) {
                std::cerr << "  -> Primary fit failed, skipping bin." << std::endl;
                continue;
            }

            auto& frac_Be_var = processor->getFracBeVar();
            auto& frac_B_var = processor->getFracBVar();
            h_frac_be_val->SetBinContent(binIndex + 1, frac_Be_var.getVal());
            h_frac_be_val->SetBinError(binIndex + 1, frac_Be_var.getError());
            h_frac_b_val->SetBinContent(binIndex + 1, frac_B_var.getVal());
            h_frac_b_val->SetBinError(binIndex + 1, frac_B_var.getError());

            auto narrowFractions = processor->getPrimaryNarrowFractions();
            double Nsignal = processor->getSignalEventsInRange();
            
            h_b_fraction->SetBinContent(binIndex + 1, narrowFractions.b_frac);
            h_be_fraction->SetBinContent(binIndex + 1, narrowFractions.be_frac);
            h_c_fraction->SetBinContent(binIndex + 1, narrowFractions.c_frac);
            if (Nsignal > 0) {
                h_b_fraction->SetBinError(binIndex + 1, sqrt(narrowFractions.b_frac * (1 - narrowFractions.b_frac) / Nsignal));
                h_be_fraction->SetBinError(binIndex + 1, sqrt(narrowFractions.be_frac * (1 - narrowFractions.be_frac) / Nsignal));
                h_c_fraction->SetBinError(binIndex + 1, sqrt(narrowFractions.c_frac * (1 - narrowFractions.c_frac) / Nsignal));
            }
            h_signal_events->SetBinContent(binIndex + 1, Nsignal);

            TCanvas* canvas = new TCanvas(Form("c_%s_%d", detector.c_str(), binIndex), "Fit", 800, 600);
            canvas->Divide(1, 2);
            
            TPad* pad1 = (TPad*)canvas->cd(1);
            pad1->SetPad(0, 0.26, 1, 1);
            pad1->SetBottomMargin(0.02);
            pad1->SetLogy();
            
            double fit_min = processor->getChargeVar().getMin();
            double fit_max = processor->getChargeVar().getMax();
            frame->GetYaxis()->SetTitle("Events");
            frame->GetYaxis()->SetRangeUser(0.5, 2.5 * processor->getSignalHist()->GetMaximum());
            frame->GetXaxis()->SetRangeUser(fit_min, fit_max);
            frame->GetXaxis()->SetLabelSize(0);
            frame->Draw();
            
            auto legend = createLegend();
            legend->AddEntry("data_hist", "L1 Signal", "pze");
            //legend->AddEntry("Be_pdf", "Be Template", "l");
            //legend->AddEntry("B_pdf", "B Template", "l");
            //legend->AddEntry("C_pdf", "C Template", "l");
            legend->AddEntry("total_pdf", "Fit", "l");
            legend->Draw();
            
            auto info = createInfoBox(processor->getChi2NDF(), processor->getNDF(), processor->getFitResult()->status(),
                                      frac_Be_var, frac_B_var,
                                      narrowFractions, Nsignal);
            
            TPad* pad2 = (TPad*)canvas->cd(2);
            pad2->SetPad(0, 0, 1, 0.26);
            pad2->SetBottomMargin(0.3);
            pad2->SetGridy();
            pad2->SetTopMargin(0.02);
            
            auto pullGraph = std::make_unique<TGraphErrors>();
            calculatePull(frame, pullGraph.get(), fit_min, fit_max);
            setupPullPlot(pullGraph.get(), fit_min, fit_max);
            
            if (pullGraph->GetN() > 0) {
                pullGraph->Draw("AP");
                TLine zeroLine(fit_min, 0, fit_max, 0);
                zeroLine.SetLineStyle(2);
                zeroLine.SetLineColor(kGray + 1);
                zeroLine.Draw("SAME");
            }
            
            canvas->Print(outputPDF.c_str());
            canvas->Clear();
            info->Draw();
            canvas->Print(outputPDF.c_str());
            
            delete frame;
            delete canvas;
            
            std::cout << "  -> Running " << N_TOYS << " Toy MC simulations..." << std::endl;
            processor->runToySimulations(N_TOYS, toy_hists);
            std::cout << "  -> Toy MC finished." << std::endl;

            processedPlots++;
        }
    }
    
    c0->Print((outputPDF + "]").c_str());
    delete c0;

    auto fout = std::unique_ptr<TFile>(TFile::Open(outputROOT.c_str(), "RECREATE"));
    fout->cd();
    h_be_fraction->Write();
    h_b_fraction->Write();
    h_c_fraction->Write();
    h_signal_events->Write();
    h_frac_be_val->Write();
    h_frac_b_val->Write();

    for(auto const& [key, val] : toy_hists) {
        val->Write();
        delete val;
    }
    fout->Close();
    
    std::cout << "\nAnalysis complete. " << processedPlots << " plots saved to " << outputPDF << std::endl;
    std::cout << "Results and Toy MC data saved to " << outputROOT << std::endl;
}

// --- FitProcessor Method Implementations ---

FitProcessor::FitProcessor(int bin_idx, const std::string& det, TFile* l1File, TFile* l2File)
    : binIndex(bin_idx), detectorName(det), l1FilePtr(l1File), l2FilePtr(l2File),
      charge("charge", "Charge", 3.5, 6.3),
      frac_Be("frac_Be", "Be fraction", 0.2, 0.0, 0.9),
      frac_B("frac_B", "B fraction", 0.8, 0.0, 1.0) 
{
    // 【修正模型】使用 RooFormulaVar 定义 frac_C，确保总和为1
    frac_C = std::make_unique<RooFormulaVar>(
        "frac_C", "C fraction", "1.0 - frac_Be - frac_B", RooArgList(frac_Be, frac_B)
    );
}

bool FitProcessor::initialize() {
    TH2D* h2d_signal = (TH2D*)l1FilePtr->Get(("unbiasedL1Signal_Boron_" + detectorName + "_Bkg").c_str());
    if (!h2d_signal) { return false; }
    h_signal_orig = getEnergySlice(2, 2, h2d_signal, binIndex, "signal_orig");
    if (!h_signal_orig || h_signal_orig->GetEntries() < 50) { return false; }

    for (const auto& element : elements) {
        TH2D* h2d_template;
        if (element == "Carbon") {
            h2d_template = (TH2D*)l1FilePtr->Get(("L1Temp_" + element + "_" + detectorName + "_Bkg").c_str());
        } else {
            h2d_template = (TH2D*)l2FilePtr->Get(("L2TempTuned_EGE_" + element + "_" + detectorName + "_unb").c_str());
        }
        if (!h2d_template) { return false; }
        
        auto h_template = getEnergySlice(2, 2, h2d_template, binIndex, (element + "_template_orig").c_str());
        h_template->Smooth(1);
        if (!h_template || h_template->GetEntries() < 10) { return false; }
        templates_hist_orig.push_back(std::move(h_template));
    }

    for (size_t i = 0; i < templates_hist_orig.size(); ++i) {
        auto extended_hist = extendHistogram(templates_hist_orig[i].get(), charge.getMin(), charge.getMax());
        auto rooData = std::make_unique<RooDataHist>((elements[i] + "_data").c_str(), "", charge, extended_hist.get());
        template_pdfs.push_back(std::make_unique<RooHistPdf>((elements[i] + "_pdf").c_str(), "", charge, *rooData));
        template_data_hists.push_back(std::move(rooData));
    }

    total_pdf = std::make_unique<RooAddPdf>("total_pdf", "Total PDF",
                                            RooArgList(*template_pdfs[0], *template_pdfs[1], *template_pdfs[2]),
                                            RooArgList(frac_Be, frac_B, *frac_C),
                                            false); // <-- 必须为 false 或省略
    return true;
}

RooPlot* FitProcessor::runPrimaryFit() {
    data_hist = std::make_unique<RooDataHist>("data", "data", charge, h_signal_orig.get());
    
    //frac_Be.setVal(0.1);
    //frac_B.setVal(0.8);
    
    fitResult = std::unique_ptr<RooFitResult>(
        total_pdf->fitTo(*data_hist, Save(true), PrintLevel(-1), Strategy(2), Minimizer("Minuit2"))
    );

    if (fitResult->status() != 0) {
        return nullptr;
    }

    RooPlot* frame = charge.frame(Title(Form("%s, %.2f-%.2f GeV/n", detectorName.c_str(), Binning::NarrowBins[binIndex], Binning::NarrowBins[binIndex + 2])));
    data_hist->plotOn(frame, Name("data_hist"), MarkerStyle(20), MarkerSize(0.8), MarkerColor(kBlack), LineColor(kBlack), XErrorSize(0), DrawOption("PZ"));
    total_pdf->plotOn(frame, Name("total_pdf"), LineColor(kRed), LineWidth(2));
    total_pdf->plotOn(frame, Components(*template_pdfs[0]), Name("Be_pdf"), LineColor(kBlue), LineStyle(1), LineWidth(2));
    total_pdf->plotOn(frame, Components(*template_pdfs[1]), Name("B_pdf"), LineColor(kGreen+2), LineStyle(1), LineWidth(2));
    total_pdf->plotOn(frame, Components(*template_pdfs[2]), Name("C_pdf"), LineColor(kMagenta), LineStyle(1), LineWidth(2));
    
    double chi2 = calculateChi2(frame, "data_hist", "total_pdf", charge.getMin(), charge.getMax());
    int first_bin = h_signal_orig->FindBin(charge.getMin());
    int last_bin = h_signal_orig->FindBin(charge.getMax());
    ndf = (last_bin - first_bin + 1) - fitResult->floatParsFinal().getSize();
    chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
    
    primaryNarrowFractions = calculateNarrowRangeFractions(*total_pdf, frac_Be.getVal(), frac_B.getVal());
    
    int bin_low = h_signal_orig->FindBin(5.0);
    int bin_high = h_signal_orig->FindBin(5.4);
    Nsignal_range = h_signal_orig->Integral(bin_low, bin_high);

    return frame;
}

void FitProcessor::runToySimulations(int n_toys, const std::map<std::string, TH2D*>& toy_hists) {
    double y_val = toy_hists.at("frac_be")->GetYaxis()->GetBinCenter(binIndex + 1);

    for (int i = 0; i < n_toys; ++i) {
        auto toy_signal = std::unique_ptr<TH1D>((TH1D*)h_signal_orig->Clone(Form("toy_signal_%d", i)));
        toy_signal->Reset();
        toy_signal->FillRandom(h_signal_orig.get(), (int)h_signal_orig->GetEntries());

        std::vector<std::unique_ptr<RooDataHist>> toy_rdh_keepers;
        std::vector<std::unique_ptr<RooHistPdf>> toy_pdf_keepers;
        RooArgList toy_pdfs_for_model;

        for (size_t j = 0; j < templates_hist_orig.size(); ++j) {
            auto toy_h = std::unique_ptr<TH1D>((TH1D*)templates_hist_orig[j]->Clone(Form("toy_template_%d_%zu", i, j)));
            toy_h->Reset();
            toy_h->FillRandom(templates_hist_orig[j].get(), (int)templates_hist_orig[j]->GetEntries());
            
            auto extended_hist = extendHistogram(toy_h.get(), charge.getMin(), charge.getMax());
            
            std::string rdh_name = Form("toy_rdh_%d_%zu", i, j);
            std::string pdf_name = Form("toy_pdf_%d_%zu", i, j);
            
            auto toy_rdh = std::make_unique<RooDataHist>(rdh_name.c_str(), rdh_name.c_str(), charge, extended_hist.get());
            auto toy_pdf = std::make_unique<RooHistPdf>(pdf_name.c_str(), pdf_name.c_str(), charge, *toy_rdh);
            
            toy_pdfs_for_model.add(*toy_pdf);
            toy_rdh_keepers.push_back(std::move(toy_rdh));
            toy_pdf_keepers.push_back(std::move(toy_pdf));
        }
        
        RooAddPdf toy_model(Form("toy_model_%d", i), "Toy Model PDF", toy_pdfs_for_model, RooArgList(frac_Be, frac_B, *frac_C), false);
        
        RooDataHist toy_data_for_fit(Form("toy_data_for_fit_%d", i), "Toy Data For Fit", charge, toy_signal.get());
        auto toy_fit_result = std::unique_ptr<RooFitResult>(
            toy_model.fitTo(toy_data_for_fit, PrintLevel(-1), Strategy(2), Save(true))
        );
        
        if (toy_fit_result && toy_fit_result->status() == 0) {
            FractionResult toy_narrow = calculateNarrowRangeFractions(toy_model, frac_Be.getVal(), frac_B.getVal());
            toy_hists.at("frac_be")->Fill(frac_Be.getVal(), y_val);
            toy_hists.at("frac_b")->Fill(frac_B.getVal(), y_val);
            toy_hists.at("narrow_be")->Fill(toy_narrow.be_frac, y_val);
            toy_hists.at("narrow_b")->Fill(toy_narrow.b_frac, y_val);
        }
    }
}

FractionResult FitProcessor::calculateNarrowRangeFractions(RooAbsPdf& model, double frac_be_val, double frac_b_val) {
    FractionResult result;
    charge.setRange("narrow", 5.0, 5.4);
    
    auto& pdf_list = static_cast<RooAddPdf&>(model).pdfList();

    auto be_integral_obj = std::unique_ptr<RooAbsReal>(static_cast<RooAbsPdf*>(pdf_list.at(0))->createIntegral(charge, NormSet(charge), Range("narrow")));
    auto b_integral_obj  = std::unique_ptr<RooAbsReal>(static_cast<RooAbsPdf*>(pdf_list.at(1))->createIntegral(charge, NormSet(charge), Range("narrow")));
    auto c_integral_obj  = std::unique_ptr<RooAbsReal>(static_cast<RooAbsPdf*>(pdf_list.at(2))->createIntegral(charge, NormSet(charge), Range("narrow")));

    double be_integral = be_integral_obj->getVal();
    double b_integral  = b_integral_obj->getVal();
    double c_integral  = c_integral_obj->getVal();
    
    double frac_c_val = 1.0 - frac_be_val - frac_b_val;
    if (frac_c_val < 0) frac_c_val = 0;

    double total_integral = frac_be_val * be_integral + frac_b_val * b_integral + frac_c_val * c_integral;
    
    if (total_integral > 1e-9) { 
        result.be_frac = frac_be_val * be_integral / total_integral;
        result.b_frac = frac_b_val * b_integral / total_integral;
        result.c_frac = frac_c_val * c_integral / total_integral;
    }
    
    charge.setRange("full", 3.5, 6.3); 
    return result;
}

void chargeTempFit() {
    fitBoronCharge();
}