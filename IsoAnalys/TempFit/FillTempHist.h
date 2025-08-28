/***********************************************************
 *  Mass Template Generation for Isotope Analysis header file
 *  
 *  Author: Z.Yan
 *  Date: 2024.10.31
 *  Update Date: 2025.5.13 Happy Birth Day My Sweety
 *  
 *  Purpose1: Generate mass templates for isotope separation
 *  Purpose2: select pure be10, validate iss and mc consistancy
 *  Update: Backgroud Study
 ***********************************************************/

#ifndef FILL_TEMP_HIST_H
#define FILL_TEMP_HIST_H

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPDF.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>
#include <iostream>
#include <array>
#include <memory>
#include <vector>
#include <map>
#include "./richTune.cc"
#include "../Tool.h"

using namespace AMS_Iso;

namespace {
    // Cut status masks
    constexpr unsigned int TRACKER_TOF_MASK = 0x03;     // 0b000011 (bit 0和1: TOFCut)
    constexpr unsigned int TRACKER_RICH_MASK = 0x05;     // 0b000101 (bit 0和2: RICHCut)
    // TOF Bkg Study Mask: require bit 7, 8, 9, 10, 12, 21, 22, 23, 26, 27 == 1
    //constexpr unsigned int BKG_TOF_MASK = 0xCE01780;
    // RICH Bkg Study Mask: require bit 7, 8, 9, 10, 12, 21, 22, 23, 28, 29 == 1
    //constexpr unsigned int BKG_RICH_MASK = 0x30E01780;

    // basic mask: bits 7,8,9,10,12,21,22,23
    constexpr unsigned int BKG_BASIC_MASK = 0x0E01780;
    // basic mask: bits 7,8,9,10,12,16,17,18,23, plus L1 Norm cut, no l1unbias cut
    //constexpr unsigned int BKG_BASIC_MASK = 0x871780;
    //bits 7,8,9,10,11,12,13,14,15,16,17,18,23
    constexpr unsigned int L1INNER_MASK = 0x87FF80;
    //bits 7,8,9,10,11,12,13,19,20,21,22,23
    constexpr unsigned int UNBL1INNER_MASK = 0xF83F80;

    // TOF Cut mask: bits 26,27
    constexpr unsigned int BKG_TOF_MASK = 0xC000000; 
    // RICH Cut mask: bits 28,29,30
    constexpr unsigned int BKG_RICH_MASK = 0x30000000;  
    constexpr unsigned int NORM_RICH_MASK = 0x70000000;  

    // Nuclei Config
    struct IsotopeConfig {
        std::string name;           // Li, B, Be...
        std::string mcPrefix;       // MC rootfile name prefix
        std::string dataName;       // data rootfile name prefix
        int charge;                 // 
        std::vector<int> masses;    // 
        double alphaMin;            // alpha tuning by Dr.Yi Jia
        double alphaMax;
        int nAlpha;
        std::vector<int> ParticleID; 
    };

    const std::map<std::string, IsotopeConfig> isotope_configs = {
        {"Li", {
            "Li", "Li", "Lit", 3, {6, 7}, 
            0.98, 1.02, 21, {61, 62}
        }},
        {"Be", {
            "Be", "Be", "Ber", 4, {7, 9, 10}, 
            0.98, 1.02, 21, {63, 64, 114}
        }},
        {"B", {
            "B", "B", "Bor", 5, {10, 11}, 
            0.98, 1.02, 21, {65, 66}
        }}
    };
}

int getIsotopeID(const IsotopeConfig& config, int mass) {
    auto it = std::find(config.masses.begin(), config.masses.end(), mass);
    if (it != config.masses.end()) {
        size_t index = std::distance(config.masses.begin(), it);
        if (index < config.ParticleID.size()) {
            return config.ParticleID[index];
        }
    }
    return -1;
}

std::array<double, 13> kSourceFluxFitPoints{
    2., 2.6, 3.2, 4, 6, 13.0, 30, 50.0, 100.0, 211.0, 330.0, 643.0, 2000.0
};

// beta detector
enum DetectorType {TOF, NaF, AGL, DETECTOR_COUNT};
const std::string DetectorNames[DETECTOR_COUNT] = {"TOF", "NaF", "AGL"};

//using same wide ek bin now
struct HistogramSet {
    // basic: mass distribution
    std::vector<std::unique_ptr<TH2D>> massHist[DETECTOR_COUNT];
    
    // only mc: pure L2 isotope mass dis.
    std::vector<std::unique_ptr<TH2D>> origMassHist[DETECTOR_COUNT];
    
    // only iss now(25.5.13): isotope estimator
    std::vector<std::unique_ptr<TH2D>> estHist[DETECTOR_COUNT];
    
    // Pure FragNuc(eg. Be10) Mass Distri.
    std::vector<std::unique_ptr<TH2D>> pureFragNucMassHist[DETECTOR_COUNT];

    // create mc hist for certain alpha
    void createMCHistograms(int processMass, const std::vector<int>& masses, int alphaIndex, const std::vector<std::array<double, 31>>& all_ek_bins) {
        // each det.
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            //each masses
            for (size_t k = 0; k < masses.size(); ++k) {
                int mass = masses[k];
                const auto& ek_bins = all_ek_bins[k];
                //det, processMCMass, bin for mass, alpha
                massHist[det].push_back(std::make_unique<TH2D>(
                    Form("h_inv_%sMass_mass%d_bin%d_alpha%d_FragToBe", DetectorNames[det].c_str(), processMass, mass, alphaIndex),//change!!
                    Form("inv_%sMass", DetectorNames[det].c_str()), 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));

                origMassHist[det].push_back(std::make_unique<TH2D>(
                    Form("hL2original_inv_%sMass_mass%d_bin%d_alpha%d", DetectorNames[det].c_str(), processMass, mass, alphaIndex),
                    Form("inv_%sMass", DetectorNames[det].c_str()), 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));
            }
            //frag Be7 9 10
            int bemass[3] = {7,9,10};
            for (size_t k = 0; k < 3; ++k){
                pureFragNucMassHist[det].push_back(std::make_unique<TH2D>(
                    Form("h_inv_%sMass_mass%d_alpha%d_FragMass%d", DetectorNames[det].c_str(), processMass, alphaIndex, bemass[k]),//change!!
                    Form("inv_%sMass", DetectorNames[det].c_str()), 200, 0.0, 0.5, Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
                ));
            }
        }
    }
    
    // create iss hist 
    void createISSHistograms(const std::string& dataName, const std::vector<int>& masses, const std::vector<std::array<double, 31>>& all_ek_bins) {
        for (size_t j = 0; j < masses.size(); ++j) {
            const auto& ek_bins = all_ek_bins[j];
            
            // each det.
            for (int det = 0; det < DETECTOR_COUNT; det++) {
                // mass
                massHist[det].push_back(std::make_unique<TH2D>(
                    Form("h_inv_%sMass_%sUseMass%d", DetectorNames[det].c_str(), dataName.c_str(), masses[j]),
                    Form("inv_%sMass", DetectorNames[det].c_str()), 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));
                
                // est
                estHist[det].push_back(std::make_unique<TH2D>(
                    Form("h_inv_%sEst_%sUseMass%d", DetectorNames[det].c_str(), dataName.c_str(), masses[j]),
                    Form("inv_%sEst", DetectorNames[det].c_str()), 500, 0.0, 1.0, ek_bins.size()-1, ek_bins.data()
                ));
            }
        }
    }
    
    // 
    void writeToFile(TFile* outputFile, bool isMC) {
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            for (size_t i = 0; i < massHist[det].size(); ++i) {
                massHist[det][i]->Write();
                
                if (isMC && i < origMassHist[det].size()) {
                    origMassHist[det][i]->Write();
                } else if (!isMC && i < estHist[det].size()) {
                    estHist[det][i]->Write();
                }
            }
            for (size_t i = 0; i < 3; ++i){
                if(isMC) pureFragNucMassHist[det][i]->Write();
            }
        }
    }
};

// MC temp validation 
struct ValidHistograms {
    // 
    std::unique_ptr<TH2F> CutoffHist[DETECTOR_COUNT];
    std::unique_ptr<TH2F> pureCutoffHist[DETECTOR_COUNT];
    
    ValidHistograms() {
        // 
        CutoffHist[TOF] = std::make_unique<TH2F>("CutoffHist_TOF", "Pass TOF Selection;cutoff rigidity[GV];TOF #beta", 250, 0, 8, 250, 0.6, 0.9);
        CutoffHist[NaF] = std::make_unique<TH2F>("CutoffHist_NaF", "Pass NaF Selection;cutoff rigidity[GV];NaF #beta", 250, 0, 20, 250, 0.75, 1.001);
        CutoffHist[AGL] = std::make_unique<TH2F>("CutoffHist_AGL", "Pass AGL Selection;cutoff rigidity[GV];AGL #beta", 250, 4, 30, 250, 0.954, 1.001);
        
        // pure Isotope
        pureCutoffHist[TOF] = std::make_unique<TH2F>("CutoffHist_TOF_pure", "Pass TOF and Pure Iso Selection;cutoff rigidity[GV];TOF #beta", 250, 0, 8, 250, 0.6, 0.9);
        pureCutoffHist[NaF] = std::make_unique<TH2F>("CutoffHist_NaF_pure", "Pass NaF and Pure Iso Selection;cutoff rigidity[GV];NaF #beta", 250, 0, 20, 250, 0.75, 1.001);
        pureCutoffHist[AGL] = std::make_unique<TH2F>("CutoffHist_AGL_pure", "Pass AGL and Pure Iso Selection;cutoff rigidity[GV];AGL #beta", 250, 4, 30, 250, 0.954, 1.001);
    }
    
    //
    void writeToFile(TFile* outputFile) {
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            CutoffHist[det]->Write();
            pureCutoffHist[det]->Write();
        }
    }
    
    // draw
    void drawAndSave() {
        // diff beta range for diff det
        const std::array<double, DETECTOR_COUNT> betaMin = {0.6, 0.75, 0.95};
        const std::array<double, DETECTOR_COUNT> betaMax = {0.9, 1.005, 1.005};
        
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            auto canvas = std::make_unique<TCanvas>(
                Form("c_%s", DetectorNames[det].c_str()),
                Form("%s Beta vs Cutoff", DetectorNames[det].c_str()),
                1200, 900
            );
            canvas->SetRightMargin(0.18);
            canvas->SetLogz();
            
            CutoffHist[det]->Draw("colz");
            DrawTheoryLines(canvas.get(), betaMin[det], betaMax[det]);
            //canvas->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/%s_beta_cutoff.png", DetectorNames[det].c_str()));
            
            pureCutoffHist[det]->GetZaxis()->SetRangeUser(1, CutoffHist[det]->GetMaximum());
            pureCutoffHist[det]->Draw("colz");
            DrawTheoryLines(canvas.get(), betaMin[det], betaMax[det]);
            //canvas->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/%s_beta_cutoff_pure.png", DetectorNames[det].c_str()));
        }
    }
};

// Background study for isotope 
struct BackgroundHistograms { 
    static constexpr int NUM_ISOTOPES = 3; 
    const std::array<std::string, NUM_ISOTOPES> isotopeNames = {"Be7", "Be9", "Be10"};
    const std::array<std::string, 2> SourceTypeName = {"Truth", "Cut"};
    
    // Total source counts per detector
    std::array<std::array<std::unique_ptr<TH1F>, 2>, DETECTOR_COUNT> SourceCounts;
    // Total cut inner frag events counts per det., NO Q CUT FOR MASS FIT(eg. rich q) 
    std::array<std::unique_ptr<TH1F>, DETECTOR_COUNT> FragNucCounts;
    // Counts of source fragmenting to specific isotopes, per detector and isotope
    std::array<std::array<std::unique_ptr<TH1F>, NUM_ISOTOPES>, DETECTOR_COUNT> SourceToIsoCounts;
    // Ratio of fragments to source, per detector and isotope
    std::array<std::array<std::unique_ptr<TH1F>, NUM_ISOTOPES>, DETECTOR_COUNT> SourceToIsoRatio;
    // Fragmented Particle ID in each Ek bin
    std::array<std::unique_ptr<TH2F>, DETECTOR_COUNT> FragNucID_Ek;
    
    BackgroundHistograms() {
        std::cout << "Creating BackgroundHistograms..." << std::endl;
        // Initialize histograms for each detector type
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            // Total source counts in Layer 1
            for(int t= 0; t < 2; t++){
                SourceCounts[det][t] = std::make_unique<TH1F>(
                    Form("%sL1SourceCounts_%s", SourceTypeName[t].c_str(), DetectorNames[det].c_str()),
                    Form("%s Source Particle in Layer1 selected by %s; events; E_{k}/n", DetectorNames[det].c_str(), SourceTypeName[t].c_str()),
                    Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
                );
                SourceCounts[det][t]->Sumw2();
                //std::cout << "Initialized SourceCounts[" << det << "][" << t << "] -> " << SourceCounts[det][t].get() << std::endl;
            }
            //
            FragNucCounts[det] = std::make_unique<TH1F>(
                Form("FragNucCounts_%s", DetectorNames[det].c_str()),
                Form("%s Fragmented Nuclei Counts selected by cuts; events; E_{k}/n", DetectorNames[det].c_str()),
                Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
            );
            FragNucCounts[det]->Sumw2();
            
            //id ek
            FragNucID_Ek[det] = std::make_unique<TH2F>(
                Form("FragNucID_Ek_%s", DetectorNames[det].c_str()),
                Form("%s Fragmented Nuclei ID in each Ek/n Bin; ParticleID; E_{k}/n", DetectorNames[det].c_str()),
                2500, -0.5, 2499.5, Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
            );

            // Create histograms for each isotope
            for (int iso = 0; iso < NUM_ISOTOPES; iso++) {
                // Source fragments to specific isotope
                SourceToIsoCounts[det][iso] = std::make_unique<TH1F>(
                    Form("SourceToIsoCounts_%s_%s", DetectorNames[det].c_str(), isotopeNames[iso].c_str()),
                    Form("%s Source Fragment to %s; events; E_{k}/n", 
                         DetectorNames[det].c_str(), isotopeNames[iso].c_str()),
                    Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
                );
                SourceToIsoCounts[det][iso]->Sumw2();
                
                // Ratio of fragments to source events
                SourceToIsoRatio[det][iso] = std::make_unique<TH1F>(
                    Form("SourceToIsoRatio_%s_%s", DetectorNames[det].c_str(), isotopeNames[iso].c_str()),
                    Form("%s Frag%s/Source; events; E_{k}/n", 
                         DetectorNames[det].c_str(), isotopeNames[iso].c_str()),
                    Binning::NarrowBins.size()-1, Binning::NarrowBins.data()
                );
            }
        }
    }
       
    // Calculate ratio of fragments to source
    void calculateRatios() {
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            for (int iso = 0; iso < NUM_ISOTOPES; iso++) {
                // Divide fragment counts by source counts to get ratio
                SourceToIsoRatio[det][iso]->Divide(
                    SourceToIsoCounts[det][iso].get(), 
                    SourceCounts[det][0].get()
                );
            }
        }
    }

    // Write all histograms to output file
    void writeToFile(TFile* outputFile) {
        for (int det = 0; det < DETECTOR_COUNT; det++) {
            //std::cout << "Writing SourceCounts[" << det << "] histograms to file..." << std::endl;
            SourceCounts[det][0]->Write(); // L14.5-5.5, inner3.5-5.5 sel. counts, extra truth cut: L1 id=Boron for mc
            SourceCounts[det][1]->Write(); // iss: L14.5-5.5, inner3.5-4.5 sel. counts; mc: no extra truth cut 
            FragNucCounts[det]->Write();
            FragNucID_Ek[det]->Write();
            
            for (int iso = 0; iso < NUM_ISOTOPES; iso++) {
                SourceToIsoCounts[det][iso]->Write();
                SourceToIsoRatio[det][iso]->Write();
            }
        }
    }
};

// get beta and ek/n for each det.
struct DetectorMeasure {
    double beta;
    double ek;
    bool isValid;
    
    static DetectorMeasure get(DetectorType type, 
                           const std::array<double, DETECTOR_COUNT>& betas, 
                           const std::array<double, DETECTOR_COUNT>& energies) {
        DetectorMeasure data;
        data.beta = betas[type];
        data.ek = energies[type];
        data.isValid = isValidBeta(data.beta);
        return data;
    }
};

// declare
void ProcessTreeWithMass(TFile* file, const char* treeName, const char* outputFileName, 
                        bool isMC, const IsotopeConfig& config, int processMass = 0);
void BuildTempHist(const std::string& isotype);
void FillTempHist(const std::string& NucName);

#endif // FILL_TEMP_HIST_H