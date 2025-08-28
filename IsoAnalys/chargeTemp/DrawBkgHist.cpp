#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <vector>
#include <string>
#include "EventProcessor.hh"

using namespace AMS_Iso;

enum HistTypeIdx {
    L1Signal_Idx = 0,
    unbiasedL1Signal_Idx,
    L1Temp_Idx,
    L1Temp_unb_Idx,
    L2Temp_Idx,
    L2Temp_unb_Idx
};

const int NUM_HIST_TYPES = 6;

struct HistInfo {
    std::string name;
    std::string title;
    std::string xLabel;
    std::string yLabel;
    bool useUnbMask;
    HistTypeIdx idx;
};

const std::vector<HistInfo> histTypes = {
    {"L1Signal",         "L1 Q Signal",           "L1 Q",         "E_{k}/n {GeV/n}", false, L1Signal_Idx},
    {"unbiasedL1Signal", "L1 Q Signal (Unbiased)", "L1 Q (unb)",   "E_{k}/n {GeV/n}", true,  unbiasedL1Signal_Idx},
    {"L1Temp",           "L1 Q Temp",             "L1 Q",         "E_{k}/n {GeV/n}", false, L1Temp_Idx},
    {"L1Temp_unb",       "L1 Q Temp (Unbiased)",  "L1 Q (unb)",   "E_{k}/n {GeV/n}", true,  L1Temp_unb_Idx},
    {"L2Temp",           "L2 Temp Q",             "L2 Q",         "E_{k}/n {GeV/n}", false, L2Temp_Idx},
    {"L2Temp_unb",       "L2 Temp Q (Unbiased)",  "L2 Q",         "E_{k}/n {GeV/n}", true,  L2Temp_unb_Idx}
};

void DrawBkgHist() {
    double coe = 0.8;
    std::cout << "Creating background histograms using EventProcessor... coe = " << coe << std::endl;

    const auto& NarrowBins = Binning::NarrowBins;
    std::vector<double> Rbins_beta;
    for (auto ek : NarrowBins) Rbins_beta.push_back(kineticEnergyToBeta(ek));

    TH2D* hists[NUM_ELEMENTS][NUM_DETECTORS][NUM_HIST_TYPES] = {};
    for (int elem_idx = 0; elem_idx < NUM_ELEMENTS; ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int nQbins = 400;
        double qMin = elem.charge - 2;
        double qMax = elem.charge + 2;
        for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
            for (int h = 0; h < NUM_HIST_TYPES; ++h) {
                if(h < 2) {qMin = 3; qMax = 9; nQbins = 600;}
                else {qMin = elem.charge - 2; qMax = elem.charge + 2; nQbins = 400;}
                const auto& hinfo = histTypes[h];
                TString hname = Form("%s_%s_%s_Bkg", hinfo.name.c_str(), elem.name.c_str(), detNames[idet].c_str());
                TString htitle = Form("%s for %s @%s;%s;%s",
                    hinfo.title.c_str(), elem.name.c_str(), detNames[idet].c_str(),
                    hinfo.xLabel.c_str(), hinfo.yLabel.c_str());
                hists[elem_idx][idet][h] = new TH2D(hname, htitle, nQbins, qMin, qMax, NarrowBins.size()-1, NarrowBins.data());
            }
        }
    }
    std::clock_t start_clock = std::clock();
    time_t start_time = time(nullptr);
    
    TFile* inFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root");
    TTree* tree = (TTree*)inFile->Get("saveTree");
    EventProcessor ep;
    ep.setBranchAddresses(tree);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 1000000 == 0) {
            std::clock_t now_clock = std::clock();
            time_t now_time = time(nullptr);
            double elapsed_sec = double(now_clock - start_clock) / CLOCKS_PER_SEC;
            double wall_sec = difftime(now_time, start_time);
            double mem_mb = getCurrentRSS_MB();
            printf("Entry %lld / %lld | CPU(s): %.1f | Wall(s): %.1f | Mem: %.1f MB\n",
                i, nEntries, elapsed_sec, wall_sec, mem_mb);
        }
        
        //if(i > 10000) break;

        ep.getEntry(i);
        if (!ep.passBasicCut() && !ep.passUnbiasedCut()) continue;

        for (int elem_idx = 0; elem_idx < NUM_ELEMENTS; ++elem_idx) {
            const auto& elem = elements[elem_idx];
            int Z = elem.charge, A = elem.mass;

            for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
                if (!ep.passDetectorCut(idet)) continue;
                if (!ep.passBeyondCutoff(idet, Z, A, Rbins_beta)) continue;
                bool richQCut = (idet == 0) || ep.passRichQCut(Z);

                for (int htype = 0; htype < NUM_HIST_TYPES; ++htype) {
                    const auto& hinfo = histTypes[htype];
                    if (hinfo.useUnbMask && !ep.passUnbiasedCut()) continue;
                    if (!hinfo.useUnbMask && !ep.passBasicCut()) continue;

                    TH2D* hist = hists[elem_idx][idet][htype];
                    double Ek = (idet == 0) ? ep.TOFEk : ((idet == 1) ? ep.NaFEk : ep.AGLEk);

                    switch (hinfo.idx) {
                        case L1Signal_Idx:
                            if (ep.getNormalL1XY() && ep.passInnerRMSCut() && ep.passQTrkInner0L1Cut(Z))
                                hist->Fill(ep.trk_ql1, Ek);
                            break;
                        case unbiasedL1Signal_Idx:
                            if (ep.passInnerRMSCut() && ep.passQTrkInner0UnbiasedL1Cut(Z))
                                hist->Fill(ep.trk_ql1_unbias, Ek);
                            break;
                        case L1Temp_Idx:
                            if (ep.getNormalL1XY() && ep.passInnerRMSCut() && ep.passQTrkInnerCut(1, Z, coe) && ep.passTofQUpCut(Z, coe) && richQCut)
                                hist->Fill(ep.trk_ql1, Ek);
                            break;
                        case L1Temp_unb_Idx:
                            if (ep.passInnerRMSCut() && ep.passQTrkInnerCut(1, Z, coe) && ep.passTofQUpCut(Z, coe) && richQCut)
                                hist->Fill(ep.trk_ql1_unbias, Ek);
                            break;
                        case L2Temp_Idx:
                            if (ep.getNormalL1XY() && ep.getL2XY() && ep.getQl2StatusCut() && ep.passQTrkInnerCut(2, Z, coe) && ep.passTrkQL1Cut_L2Template(Z, coe) && ep.passTofQUpCut(Z, coe) && richQCut)
                                hist->Fill(ep.trk_ql2, Ek);
                            break;
                        case L2Temp_unb_Idx:
                            if (ep.getL2XY() && ep.getQl2StatusCut() && ep.passQTrkInnerCut(2, Z, coe) && ep.passTrkQL1UnbiasCut_L2Template(Z, coe) && ep.passTofQUpCut(Z, coe) && richQCut)
                                hist->Fill(ep.trk_ql2, Ek);
                            break;
                    }
                }
            }
        }
    }
    
    std::clock_t end_clock = std::clock();
    time_t end_time = time(nullptr);
    double elapsed_sec = double(end_clock - start_clock) / CLOCKS_PER_SEC;
    double wall_sec = difftime(end_time, start_time);
    double mem_mb = getCurrentRSS_MB();
    printf("Total CPU time: %.1f s, Wall time: %.1f s, Max memory: %.1f MB\n",
        elapsed_sec, wall_sec, mem_mb);

    TFile* outFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p8.root", "RECREATE");
    for (int i = 0; i < NUM_ELEMENTS; ++i)
        for (int j = 0; j < NUM_DETECTORS; ++j)
            for (int k = 0; k < NUM_HIST_TYPES; ++k)
                hists[i][j][k]->Write();

    outFile->Close();
    inFile->Close();
    std::cout << "Histograms written successfully." << std::endl;
}
