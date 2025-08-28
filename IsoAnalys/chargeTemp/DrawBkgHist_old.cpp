#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include "../Tool.h"
using namespace AMS_Iso;

void DrawBkgHist_old() {
    std::cout << "Creating background histograms from saveTree..." << std::endl;

    // ====== Constants and Settings ======
    struct ElementInfo {
        int charge;    // Z
        int mass;      // A
        std::string name;
    };
    const std::vector<ElementInfo> elements = {
        {4,  7,  "Beryllium"},
        {5, 10,  "Boron"},
        {6, 12,  "Carbon"},
        {7, 14,  "Nitrogen"},
        {8, 16,  "Oxygen"}
    };
    const std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};
    constexpr int NDET = 3;
    constexpr double SAFETY_FACTORS[NDET] = {1.06, 1.005, 1.0005}; // TOF, NaF, AGL

    // Cut parameters
    constexpr double Trk_InnerQCut = 0.45;
    constexpr double Trk_L1QHighCut = 0.65;
    constexpr double Trk_L1QLowCut  = 0.46;
    constexpr double Tof_QLowCut    = 0.6;
    constexpr double Tof_QHighCut   = 1.2;

    // Bit masks for basic and detector cuts
    constexpr unsigned int BKG_BASIC_MASK     = 0x870780;
    constexpr unsigned int BKG_UNB_BASIC_MASK = 0xE00780; // For unbiased L1
    constexpr unsigned int BKG_TOF_MASK   = 0xC000000;    // Bits 26,27
    constexpr unsigned int BKG_RICH_MASK  = 0x30000000;   // Bits 28,29

    // Histogram binning
    const int nQbins = 600;
    const double qMin = 3;
    const double qMax = 9;
    const auto& wideBins = Binning::WideBins;

    // ====== Beta Bins Preparation ======
    std::vector<double> Rbins_beta;
    for (auto ek : wideBins)
        Rbins_beta.push_back(kineticEnergyToBeta(ek));

    // Histogram types and their info (name, title, xLabel, yLabel, useUnbL1Q, useUnbMask)
    struct HistInfo {
        std::string name;
        std::string title;
        std::string xLabel;
        std::string yLabel;
        bool useUnbL1Q;   // If true, fill with trk_ql1_unbias or trk_ql2_unbias
        bool useUnbMask;  // If true, use unbiased mask; else use basic mask
    };
    const std::vector<HistInfo> histTypes = {
        //           name                title                        xLabel      yLabel         useUnbL1Q useUnbMask
        {"L1Signal",     "L1 Q Signal",   "L1 Q",     "E_{k}/n [GeV/n]",    false,  false},
        {"unbiasedL1Signal", "L1 Q Signal (Unbiased)", "L1 Q (unb)", "E_{k}/n [GeV/n]",    true,   true},
        {"L1Temp",       "L1 Q Temp",     "L1 Q",     "E_{k}/n [GeV/n]",    false,  false},
        {"L1Temp_unb",   "L1 Q Temp (Unbiased)", "L1 Q (unb)", "E_{k}/n [GeV/n]",    true,   true},
        {"L2Temp",       "L2 Temp Q",     "L2 Q",     "E_{k}/n [GeV/n]",    false,  false},
        {"L2Temp_unb",   "L2 Temp Q (UnbiasedL1Cut)", "L2 Q", "E_{k}/n [GeV/n]", false,  true}
    };

    // ====== Utility Lambdas ======
    auto in_cut = [](double val, double z, double coe, double low, double high) {
        return val > z - coe*low && val < z + coe*high;
    };

    // ====== Input/Output File Handling ======
    TFile* inFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root");
    if (!inFile || inFile->IsZombie()) {
        std::cout << "Error opening input file" << std::endl;
        return;
    }
    TTree* saveTree = (TTree*)inFile->Get("saveTree");
    if (!saveTree) {
        std::cout << "Error: Cannot find saveTree in the file" << std::endl;
        inFile->Close();
        delete inFile;
        return;
    }
    TFile* outFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist.root", "RECREATE");

    // ====== Histogram Creation: Use a single loop to avoid redundancy ======
    std::map<std::string, TH2D*> histMap;
    for (const auto& elem : elements) {
        for (int idet = 0; idet < NDET; ++idet) {
            std::string elemName = elem.name;
            std::string det = detNames[idet];
            for (const auto& hinfo : histTypes) {
                std::string hname = Form("%s_%s_%s_Bkg", hinfo.name.c_str(), elemName.c_str(), det.c_str());
                std::string htitle = Form("%s for %s @%s;%s;%s",
                        hinfo.title.c_str(), elemName.c_str(), det.c_str(),
                        hinfo.xLabel.c_str(), hinfo.yLabel.c_str());
                histMap[hname] = new TH2D(hname.c_str(), htitle.c_str(),
                    nQbins, qMin, qMax, wideBins.size()-1, wideBins.data());
            }
        }
    }

    // ====== Branch Address Setup ======
    Float_t trk_ql1, trk_ql1_unbias, trk_ql2;
    Float_t QTrkInner[3], rich_q[2], tof_ql[4], tk_qrmn[2][3];
    bool ChargeCutsSelectJudge[5][3];
    bool betacutoffcut[5];
    double usedEk, TOFEk, NaFEk, AGLEk;
    double TOFBeta, NaFBeta, AGLBeta, cutOffRig;
    unsigned int cutStatus = 0;
    bool isMC = false;

    saveTree->SetBranchAddress("trk_ql1", &trk_ql1);
    saveTree->SetBranchAddress("trk_ql1_unbias", &trk_ql1_unbias);
    saveTree->SetBranchAddress("trk_ql2", &trk_ql2);
    saveTree->SetBranchAddress("QTrkInner", QTrkInner);
    saveTree->SetBranchAddress("rich_q", rich_q);
    saveTree->SetBranchAddress("tof_ql", tof_ql);
    saveTree->SetBranchAddress("tk_qrmn", tk_qrmn);
    saveTree->SetBranchAddress("ChargeCutsSelectJudge", ChargeCutsSelectJudge);
    saveTree->SetBranchAddress("betacutoffcut", betacutoffcut);
    saveTree->SetBranchAddress("usedEk", &usedEk);
    saveTree->SetBranchAddress("TOFEk", &TOFEk);
    saveTree->SetBranchAddress("NaFEk", &NaFEk);
    saveTree->SetBranchAddress("AGLEk", &AGLEk);
    saveTree->SetBranchAddress("TOFBeta", &TOFBeta);
    saveTree->SetBranchAddress("NaFBeta", &NaFBeta);
    saveTree->SetBranchAddress("AGLBeta", &AGLBeta);
    saveTree->SetBranchAddress("cutOffRig", &cutOffRig);
    saveTree->SetBranchAddress("cutStatus", &cutStatus);

    // ====== Entry Loop ======
    Long64_t nEntries = saveTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        if (i % 1000000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << std::endl;
        }

        if(i>1000001) break;

        saveTree->GetEntry(i);

        // Standard and unbiased masks
        bool passbasic = ((cutStatus & BKG_BASIC_MASK) == BKG_BASIC_MASK);
        bool passunb   = ((cutStatus & BKG_UNB_BASIC_MASK) == BKG_UNB_BASIC_MASK);

        if(!passbasic && !passunb) continue;

        bool passtof  = (cutStatus & BKG_TOF_MASK) == BKG_TOF_MASK;
        bool passrich = (cutStatus & BKG_RICH_MASK) == BKG_RICH_MASK;

        double beta_det[NDET] = {TOFBeta, NaFBeta, AGLBeta};
        double ek_det[NDET]   = {TOFEk, NaFEk, AGLEk};
        bool pass_det[NDET] = {
            isValidBeta(TOFBeta) && passtof,
            isValidBeta(NaFBeta) && passrich,
            isValidBeta(AGLBeta) && passrich,
        };
        
        double tof_qup  = 0.5 * (tof_ql[0] + tof_ql[1]); if(tof_ql[0]*tof_ql[1]==0) tof_qup = tof_qup*2; 
        double tof_qlow = 0.5 * (tof_ql[2] + tof_ql[3]); if(tof_ql[2]*tof_ql[3]==0) tof_qlow = tof_qlow*2; 

        for (int idx = 0; idx < elements.size(); idx++) {
            const auto& elem = elements[idx];
            int z = elem.charge;
            int mass = elem.mass;
            std::string elemName = elem.name;

            
            for (int idet = 0; idet < NDET; ++idet) {
                std::string det = detNames[idet];
                if (!pass_det[idet]) continue;

                int betaBin = findBin(Rbins_beta, beta_det[idet]);
                bool cutoffcut = (betaBin >= 0) ?
                    isBeyondCutoff(Rbins_beta[betaBin], cutOffRig, SAFETY_FACTORS[idet], z, mass, isMC) : false;
                if (!cutoffcut) continue;

                bool richqcut = idet == 0 ? true : sqrt(rich_q[0]) > z - 1 && sqrt(rich_q[0]) < z + 2 && tof_qlow > z - 0.6;

                // Loop over histogram types, decide which to fill
                for (const auto& hinfo : histTypes) {
                    // Mask selection
                    bool passmask = hinfo.useUnbMask ? passunb : passbasic;
                    if (!passmask) continue;

                    std::string hname = Form("%s_%s_%s_Bkg", hinfo.name.c_str(), elemName.c_str(), det.c_str());

                    // The cuts and fill value for each histogram type
                    if (hinfo.name == "L1Signal") {
                        if (ChargeCutsSelectJudge[idx][0] &&
                            QTrkInner[0] > 3.5 && QTrkInner[0] < z+0.5)
                            histMap[hname]->Fill(trk_ql1, ek_det[idet]);
                    }
                    else if (hinfo.name == "unbiasedL1Signal") {
                        if (ChargeCutsSelectJudge[idx][0] &&
                            QTrkInner[0] > 3.45 && QTrkInner[0] < z+0.45)
                            histMap[hname]->Fill(trk_ql1_unbias, ek_det[idet]);
                    }
                    else if (hinfo.name == "L1Temp") {
                        if (ChargeCutsSelectJudge[idx][1] &&
                            QTrkInner[1] > z-0.3 && QTrkInner[1] < z+0.3 &&
                            in_cut(tof_qup, z, 0.4, Tof_QLowCut, Tof_QHighCut))
                            histMap[hname]->Fill(trk_ql1, ek_det[idet]);
                    }
                    else if (hinfo.name == "L1Temp_unb") {
                        if (ChargeCutsSelectJudge[idx][1] &&
                            QTrkInner[1] > z-0.3 && QTrkInner[1] < z+0.3 &&
                            in_cut(tof_qup, z, 0.4, Tof_QLowCut, Tof_QHighCut))
                            histMap[hname]->Fill(trk_ql1_unbias, ek_det[idet]);
                    }
                    else if (hinfo.name == "L2Temp") {
                        if (ChargeCutsSelectJudge[idx][2] &&
                            QTrkInner[2] > z-0.3 && QTrkInner[2] < z+0.3 &&
                            in_cut(trk_ql1, z, 0.4, Trk_L1QLowCut, Trk_L1QHighCut) &&
                            in_cut(tof_qup, z, 0.4, Tof_QLowCut, Tof_QHighCut))
                            histMap[hname]->Fill(trk_ql2, ek_det[idet]);
                    }
                    else if (hinfo.name == "L2Temp_unb") {
                        if (ChargeCutsSelectJudge[idx][2] &&
                            QTrkInner[2] > z-0.3 && QTrkInner[2] < z+0.3 &&
                            in_cut(trk_ql1_unbias, z, 0.4, Trk_L1QLowCut, Trk_L1QHighCut) &&
                            in_cut(tof_qup, z, 0.4, Tof_QLowCut, Tof_QHighCut) &&
                            richqcut)
                            histMap[hname]->Fill(trk_ql2, ek_det[idet]);
                    }
                }
            }
        }
    }

    // ====== Write and Clean Up ======
    outFile->cd();
    for (auto& kv : histMap) kv.second->Write();
    outFile->Close();
    inFile->Close();
    delete inFile;
    delete outFile;
    std::cout << "Background histograms creation completed." << std::endl;
}