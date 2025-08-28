#ifndef EVENT_PROCESSOR_HH
#define EVENT_PROCESSOR_HH

#include <iostream>
#include <vector>
#include <string>
#include <bitset>
#include <TTree.h>
#include <cmath>
#include "./Tool.h"

namespace AMS_Iso {

struct ElementInfo {
    int charge;
    int mass;
    std::string name;
};

const std::vector<ElementInfo> elements = {
    {4,  7,  "Beryllium"},
    {5, 10,  "Boron"},
    {6, 12,  "Carbon"}
};

const std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};
constexpr int NUM_ELEMENTS = 3;
constexpr int NUM_DETECTORS = 3;
constexpr double SAFETY_FACTORS[NUM_DETECTORS] = {1.06, 1.005, 1.0005};
constexpr double Trk_L1QLowCut  = 0.46;
constexpr double Trk_L1QHighCut = 0.65;
constexpr double Tof_QLowCut    = 0.6;
constexpr double Tof_QHighCut   = 1.2;
constexpr unsigned int BKG_BASIC_MASK     = 0x870780;
constexpr unsigned int BKG_UNB_BASIC_MASK = 0xE00780;
constexpr unsigned int BKG_TOF_MASK   = 0xC000000;
constexpr unsigned int BKG_RICH_MASK  = 0x30000000;

class EventProcessor {
public:
    Float_t trk_ql1, trk_ql1_unbias, trk_ql2;
    Float_t QTrkInner[3], rich_q[2], tof_ql[4], tk_qrmn[2][3], tk_qln[2][9][3];
    int tk_qls[9], tk_hitb[2];
    bool ChargeCutsSelectJudge[5][3];
    bool betacutoffcut[5];
    double usedEk, TOFEk, NaFEk, AGLEk, L1InnerRig;
    double TOFBeta, NaFBeta, AGLBeta, cutOffRig;
    unsigned int cutStatus = 0;
    bool isMC = false;

private:
    TTree* fTree = nullptr;

    bool in_cut(double val, double z_val, double coe, double low, double high) const {
        return val > z_val - coe * low && val < z_val + coe * high;
    }

public:
    void setBranchAddresses(TTree* tree) {
        fTree = tree;
        if (!fTree) return;

        fTree->SetBranchAddress("trk_ql1", &trk_ql1);
        fTree->SetBranchAddress("trk_ql1_unbias", &trk_ql1_unbias);
        fTree->SetBranchAddress("trk_ql2", &trk_ql2);
        fTree->SetBranchAddress("L1InnerRig", &L1InnerRig);
        fTree->SetBranchAddress("usedEk", &usedEk);
        fTree->SetBranchAddress("TOFEk", &TOFEk);
        fTree->SetBranchAddress("NaFEk", &NaFEk);
        fTree->SetBranchAddress("AGLEk", &AGLEk);
        fTree->SetBranchAddress("TOFBeta", &TOFBeta);
        fTree->SetBranchAddress("NaFBeta", &NaFBeta);
        fTree->SetBranchAddress("AGLBeta", &AGLBeta);
        fTree->SetBranchAddress("cutOffRig", &cutOffRig);
        fTree->SetBranchAddress("cutStatus", &cutStatus);

        fTree->SetBranchAddress("QTrkInner", QTrkInner);
        fTree->SetBranchAddress("rich_q", rich_q);
        fTree->SetBranchAddress("tof_ql", tof_ql);
        fTree->SetBranchAddress("tk_qls", tk_qls);
        fTree->SetBranchAddress("tk_hitb", tk_hitb);
        fTree->SetBranchAddress("tk_qrmn", tk_qrmn);
        //fTree->SetBranchAddress("tk_qln", tk_qln);
        fTree->SetBranchAddress("ChargeCutsSelectJudge", ChargeCutsSelectJudge);
        fTree->SetBranchAddress("betacutoffcut", betacutoffcut);
    }

    void getEntry(Long64_t i) {
        if (fTree) fTree->GetEntry(i);
    }

    bool passBasicCut() const { return (cutStatus & BKG_BASIC_MASK) == BKG_BASIC_MASK; }
    bool passUnbiasedCut() const { return (cutStatus & BKG_UNB_BASIC_MASK) == BKG_UNB_BASIC_MASK; }
    bool passTofCut() const { return (cutStatus & BKG_TOF_MASK) == BKG_TOF_MASK; }
    bool passRichCut() const { return (cutStatus & BKG_RICH_MASK) == BKG_RICH_MASK; }

    bool getNormalL1XY() const { return std::bitset<32>(tk_hitb[0]).test(0); }
    bool getL2XY() const { return std::bitset<32>(tk_hitb[0]).test(1); }
    bool getQl2StatusCut() const { return ((tk_qls[1] & 0x10013D) == 0); }

    double getTofQUp() const {
        double val = 0.5 * (tof_ql[0] + tof_ql[1]);
        return (tof_ql[0] * tof_ql[1] == 0) ? val * 2 : val;
    }

    double getTofQLow() const {
        double val = 0.5 * (tof_ql[2] + tof_ql[3]);
        return (tof_ql[2] * tof_ql[3] == 0) ? val * 2 : val;
    }

    double getRichQ() const { return sqrt(rich_q[0]); }

    bool passDetectorCut(int idet) const {
        if (idet == 0) return isValidBeta(TOFBeta) && passTofCut();
        if (idet == 1) return isValidBeta(NaFBeta) && passRichCut();
        if (idet == 2) return isValidBeta(AGLBeta) && passRichCut();
        return false;
    }

    bool passBeyondCutoff(int idet, int z, int mass, const std::vector<double>& rbins_beta) const {
        double beta = (idet == 0) ? TOFBeta : ((idet == 1) ? NaFBeta : AGLBeta);
        int betaBin = findBin(rbins_beta, beta);
        return (betaBin >= 0) ? isBeyondCutoff(rbins_beta[betaBin], cutOffRig, SAFETY_FACTORS[idet], z, mass, isMC) : false;
    }

    bool passInnerRMSCut() const {
        return tk_qrmn[0][2] < 0.55;
    }

    bool passQTrkInner0UnbiasedL1Cut(int z) const {
        return QTrkInner[0] > 3.45 && QTrkInner[0] < z + 0.45;
    }

    bool passQTrkInner0L1Cut(int z) const {
        return QTrkInner[0] > 3.5 && QTrkInner[0] < z + 0.5;
    }

    bool passQTrkInnerCut(int index, int z, double coe = 0.5, double low = 0.45, double high = 0.45) const {
        return in_cut(QTrkInner[index], z, coe, low, high);
    }

    bool passTofQUpCut(int z, double coe = 0.5, double low = Tof_QLowCut, double high = Tof_QHighCut) const {
        return in_cut(getTofQUp(), z, coe, low, high);
    }

    bool passTrkQL1Cut_L2Template(int z, double coe = 0.5, double low = Trk_L1QLowCut, double high = Trk_L1QHighCut) const {
        return in_cut(trk_ql1, z, coe, low, high);
    }

    bool passTrkQL1UnbiasCut_L2Template(int z, double coe = 0.5, double low = Trk_L1QLowCut, double high = Trk_L1QHighCut) const {
        return in_cut(trk_ql1_unbias, z, coe, low, high);
    }

    bool passRichQCut(int z) const {
        return getRichQ() > z - 1 && getRichQ() < z + 2 && getTofQLow() > z - 0.6;
    }
};

} // namespace AMS_Iso

#endif // EVENT_PROCESSOR_HH
