#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>
#include <TMath.h>
#include "../Tool.h"
using namespace AMS_Iso;

// PDF parameter mapping
const std::vector<std::pair<std::string, int>> LG_paramList = {
    {"Width", 0}, {"MPV", 1}, {"Area", 2}, {"Sigma", 3}
};
const std::vector<std::pair<std::string, int>> EGE_paramList = {
    {"Peak", 0}, {"SigmaL", 1}, {"AlphaL", 2},
    {"SigmaR", 3}, {"AlphaR", 4}, {"Norm", 5}, {"xmin", 6}, {"xmax", 7}
};

const std::map<int, std::string> elemMap = {{4, "Beryllium"}};//, {5, "Boron"}, {6, "Carbon"}, {7, "Nitrogen"}, {8, "Oxygen"}};
const std::vector<int> charges = {4};//{4, 5, 6, 7, 8};
const std::vector<int> masses  = {7};//{7, 10, 12, 14, 16};
const std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};
const int NDET = 3;
constexpr double SAFETY_FACTORS[NDET] = {1.06, 1.005, 1.0005};

const std::string paramFileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_1010.root";
const std::string dataFileName  = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/BeToOxy_Bkg.root";
const std::string outFileName_LG  = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/L2TempTuned_LG_Be.root";
const std::string outFileName_EGE = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/L2TempTuned_EGE_Be.root";
const std::string outFileName_EGELinear = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/L2TempTuned_EGELinear_Be.root";

// Landau-Gauss convolution
Double_t langaufun(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.3989422804014;
    Double_t mpshift = -0.22278298;
    Double_t np = 600.0;
    Double_t sc = 5.0;
    Double_t xx, mpc, fland, sum = 0.0, xlow, xupp, step, i;
    mpc = par[1] - mpshift * par[0];
    xlow = x[0]*par[4] - sc*par[3];
    xupp = x[0]*par[4] + sc*par[3];
    step = (xupp - xlow) / np;
    for (i = 1.0; i <= np/2; i++) {
        xx = xlow + (i-.5)*step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0]*par[4], xx, par[3]);
        xx = xupp - (i-.5)*step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0]*par[4], xx, par[3]);
    }
    return (par[2] * step * sum * invsq2pi / par[3]);
}

// ExpGausExp
double funcExpGausExp(double *x, double *p) {
    constexpr double sqrtPiOver2 = 1.2533141373;
    constexpr double sqrt2 = 1.4142135624;
    const double &x0 = p[0], &sigmaL = p[1], &alphaL = p[2], &sigmaR = p[3], &alphaR = p[4], &norm = p[5], &xmin = p[6], &xmax = p[7];
    double t = (x[0] - x0) / (x[0] < x0 ? sigmaL : sigmaR);
    double v = 0;
    if (t < -alphaL) {
        double a = 0.5 * alphaL * alphaL, b = alphaL * t;
        v = std::exp(a + b);
    } else if (t <= alphaR) {
        v = std::exp(-0.5 * t * t);
    } else {
        double a = 0.5 * alphaR * alphaR, b = alphaR * (-t);
        v = std::exp(a + b);
    }
    double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
    double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
    double sum = 0;
    if (tmin < -alphaL) {
        double a = 0.5 * alphaL * alphaL;
        double lv = tmin, uv = std::min(tmax, -alphaL);
        sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
    }
    if (tmax > alphaR) {
        double a = 0.5 * alphaR * alphaR;
        double lv = std::max(tmin, alphaR), uv = tmax;
        sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
    }
    if (tmin < alphaR && tmax > -alphaL) {
        double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
        double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
        sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
    }
    return norm * (v / sum);
}

// Parameter loading helpers
void getLGParams(TFile* fin, const std::string& elem, const std::string& det, const std::string& temp, int bin, double* par) {
    for (const auto& p : LG_paramList) {
        std::string hname = elem + "_" + det + "_" + temp + "_LG_" + p.first;
        TH1* h = (TH1*)fin->Get(hname.c_str());
        par[p.second] = h ? h->GetBinContent(bin) : 0;
    }
    par[4] = 1.0;
}

void getEGEParams(TFile* fin, const std::string& elem, const std::string& det, const std::string& temp, int bin, double* par) {
    for (const auto& p : EGE_paramList) {
        std::string hname = elem + "_" + det + "_" + temp + "_EGE_" + p.first;
        TH1* h = (TH1*)fin->Get(hname.c_str());
        par[p.second] = h ? h->GetBinContent(bin) : 0;
    }
}
                       
// Numerical inverse CDF
double InvCDF(TF1& f, double target_cdf, double xlow, double xhigh, double tol=1e-4) {
    double left = xlow, right = xhigh;
    double total = f.Integral(xlow, xhigh, 1e-7);
    if (total <= 0) return xlow;
    while (right - left > tol) {
        double mid = 0.5 * (left + right);
        double cdf = f.Integral(xlow, mid, 1e-7) / total;
        if (cdf < target_cdf) left = mid;
        else right = mid;
    }
    return 0.5 * (left + right);
}

void L2TempTuning() {
    TFile finParam(paramFileName.c_str());
    TFile inFile(dataFileName.c_str());
    if (finParam.IsZombie() || inFile.IsZombie()) {
        std::cerr << "Error opening input files!" << std::endl; 
        return;
    }

    // Prepare output histograms
    std::map<std::string, std::unique_ptr<TH2D>> h2_LG, h2_EGE, h2_EGELinear;
    for (int z : charges) {
        std::string elemName = elemMap.at(z);
        for (const std::string& det : detNames) {
            std::string key = elemName + "_" + det;
            h2_LG[key] = std::make_unique<TH2D>(
                Form("L2TempTuned_LG_%s", key.c_str()),
                Form("L2 Q LG Tuned %s;Q;E_{k}/n", key.c_str()), 
                200, z-1, z+1, Binning::WideBins.size()-1, Binning::WideBins.data());
            h2_EGE[key] = std::make_unique<TH2D>(
                Form("L2TempTuned_EGE_%s", key.c_str()),
                Form("L2 Q EGE Tuned %s;Q;E_{k}/n", key.c_str()), 
                200, z-1, z+1, Binning::WideBins.size()-1, Binning::WideBins.data());
            h2_EGELinear[key] = std::make_unique<TH2D>(
                Form("L2TempTuned_EGELinear_%s", key.c_str()),
                Form("L2 Q EGE Linear Tuned %s;Q;E_{k}/n", key.c_str()), 
                200, z-1, z+1, Binning::WideBins.size()-1, Binning::WideBins.data());
        }
    }

    // Precompute PDFs
    struct PdfSet { TF1 lg_l1, lg_l2, ege_l1, ege_l2; };
    std::map<std::string, std::vector<PdfSet>> pdf_map;
    
    for (int z : charges) {
        std::string elemName = elemMap.at(z);
        for (const std::string& det : detNames) {
            std::vector<PdfSet> pdfvec;
            int nEk = Binning::WideBins.size()-1;
            for (int iy = 1; iy <= nEk; ++iy) {
                double par_lg_l1[5]={0}, par_lg_l2[5]={0}, par_ege_l1[8]={0}, par_ege_l2[8]={0};
                getLGParams(&finParam, elemName, det, "L1Temp", iy, par_lg_l1);
                getLGParams(&finParam, elemName, det, "L2Temp", iy, par_lg_l2);
                getEGEParams(&finParam, elemName, det, "L1Temp", iy, par_ege_l1);
                getEGEParams(&finParam, elemName, det, "L2Temp", iy, par_ege_l2);

                PdfSet pset = {
                    TF1(Form("lg_l1_%d",iy), langaufun, z-1, z+1, 5),
                    TF1(Form("lg_l2_%d",iy), langaufun, z-1, z+1, 5),
                    TF1(Form("ege_l1_%d",iy), funcExpGausExp, z-1, z+1, 8),
                    TF1(Form("ege_l2_%d",iy), funcExpGausExp, z-1, z+1, 8)
                };
                pset.lg_l1.SetParameters(par_lg_l1); 
                pset.lg_l2.SetParameters(par_lg_l2);
                pset.ege_l1.SetParameters(par_ege_l1); 
                pset.ege_l2.SetParameters(par_ege_l2);
                pdfvec.push_back(pset);
            }
            pdf_map[elemName+"_"+det] = std::move(pdfvec);
        }
    }

    // Load data tree
    TTree* saveTree = (TTree*)inFile.Get("saveTree");
    if (!saveTree) { 
        std::cerr << "Error: Cannot find saveTree" << std::endl; 
        return; 
    }
    
    Float_t trk_ql1, trk_ql2, QTrkInner[3];
    bool ChargeCutsSelectJudge[5][3];
    double TOFEk, NaFEk, AGLEk, TOFBeta, NaFBeta, AGLBeta, cutOffRig;
    unsigned int cutStatus = 0;
    bool isMC = false;
    
    saveTree->SetBranchAddress("trk_ql1", &trk_ql1);
    saveTree->SetBranchAddress("trk_ql2", &trk_ql2);
    saveTree->SetBranchAddress("QTrkInner", QTrkInner);
    saveTree->SetBranchAddress("ChargeCutsSelectJudge", ChargeCutsSelectJudge);
    saveTree->SetBranchAddress("TOFEk", &TOFEk);
    saveTree->SetBranchAddress("NaFEk", &NaFEk);
    saveTree->SetBranchAddress("AGLEk", &AGLEk);
    saveTree->SetBranchAddress("TOFBeta", &TOFBeta);
    saveTree->SetBranchAddress("NaFBeta", &NaFBeta);
    saveTree->SetBranchAddress("AGLBeta", &AGLBeta);
    saveTree->SetBranchAddress("cutOffRig", &cutOffRig);
    saveTree->SetBranchAddress("cutStatus", &cutStatus);

    // Build beta bins
    std::vector<double> Rbins_beta;
    for (size_t ip = 0; ip < Binning::WideBins.size(); ++ip)
        Rbins_beta.push_back(kineticEnergyToBeta(Binning::WideBins[ip]));

    std::vector<double> WideBinsVec(Binning::WideBins.begin(), Binning::WideBins.end());
    
    unsigned int BKG_BASIC_MASK = 0x870780;
    unsigned int BKG_TOF_MASK   = 0xC000000;
    unsigned int BKG_RICH_MASK  = 0x30000000;

    Long64_t nEntries = saveTree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        if (i % 1000000 == 0) std::cout << "Entry " << i << "/" << nEntries << std::endl;
        //if (i > 300000) break;  // Debug limit
       
        saveTree->GetEntry(i);
        
        bool passbasic = (cutStatus & BKG_BASIC_MASK) == BKG_BASIC_MASK;
        if (!passbasic) continue;
        
        bool passtof  = (cutStatus & BKG_TOF_MASK) == BKG_TOF_MASK;
        bool passrich = (cutStatus & BKG_RICH_MASK) == BKG_RICH_MASK;
        
        double beta_det[NDET] = {TOFBeta, NaFBeta, AGLBeta};
        double ek_det[NDET]   = {TOFEk, NaFEk, AGLEk};
        bool pass_det[NDET] = {
            isValidBeta(TOFBeta) && passtof && TOFEk >= 0.61 && TOFEk <= 1.7,
            isValidBeta(NaFBeta) && passrich && NaFEk >= 0.61 && NaFEk <= 4.0,
            isValidBeta(AGLBeta) && passrich && AGLEk > 2.57 && AGLEk <= 16.3,
        };

        for (int idx = 0; idx < 1; idx++) {
            int z = charges[idx], mass = masses[idx];
            if (trk_ql2 > z+1 || trk_ql2 < z-1) continue;
            
            std::string elemName = elemMap.at(z);
            for (int idet = 0; idet < NDET; ++idet) {
                if (!pass_det[idet]) continue;
                
                int betaBin = findBin(Rbins_beta, beta_det[idet]);
                bool cutoffcut = (betaBin >= 0) ?
                    isBeyondCutoff(Rbins_beta[betaBin], cutOffRig, SAFETY_FACTORS[idet], z, mass, isMC) : false;
                if (!cutoffcut) continue;

                // L2Temp cut
                if (!(ChargeCutsSelectJudge[idx][2] && QTrkInner[2] > z-0.4 && QTrkInner[2] < z+0.4))
                    continue;

                // Find energy bin
                int ybin = findBin2(WideBinsVec, ek_det[idet]);
                if (ybin < 0 || ybin >= int(Binning::WideBins.size())-1) continue;
                
                std::string key = elemName + "_" + detNames[idet];
                PdfSet& pset = pdf_map[key][ybin];
                
                // LG mapping
                double total_LG = pset.lg_l2.Integral(pset.lg_l2.GetXmin(), pset.lg_l2.GetXmax(), 1e-7);
                double cdf_LG = total_LG > 0 ? pset.lg_l2.Integral(pset.lg_l2.GetXmin(), trk_ql2, 1e-7) / total_LG : 0;
                double qquantile_LG = pset.lg_l1.GetXmin();
                if (cdf_LG > 0 && cdf_LG < 1)
                    qquantile_LG = InvCDF(pset.lg_l1, cdf_LG, pset.lg_l1.GetXmin(), pset.lg_l1.GetXmax());
                h2_LG[key]->Fill(qquantile_LG, ek_det[idet]);
                
                // EGE mapping
                double total_EGE = pset.ege_l2.Integral(pset.ege_l2.GetXmin(), pset.ege_l2.GetXmax(), 1e-7);
                double cdf_EGE = total_EGE > 0 ? pset.ege_l2.Integral(pset.ege_l2.GetXmin(), trk_ql2, 1e-7) / total_EGE : 0;
                double qquantile_EGE = pset.ege_l1.GetXmin();
                if (cdf_EGE > 0 && cdf_EGE < 1)
                    qquantile_EGE = InvCDF(pset.ege_l1, cdf_EGE, pset.ege_l1.GetXmin(), pset.ege_l1.GetXmax());
                h2_EGE[key]->Fill(qquantile_EGE, ek_det[idet]);
                
                //EGE Linear Mapping
                double peak_L1 = pset.ege_l1.GetParameter(0);
                double peak_L2 = pset.ege_l2.GetParameter(0);
                double sigmaL_L1 = pset.ege_l1.GetParameter(1);
                double sigmaR_L1 = pset.ege_l1.GetParameter(3);
                double sigmaL_L2 = pset.ege_l2.GetParameter(1);
                double sigmaR_L2 = pset.ege_l2.GetParameter(3);

                double sigma_L1 = (trk_ql2 < peak_L2) ? sigmaL_L1 : sigmaR_L1;
                double sigma_L2 = (trk_ql2 < peak_L2) ? sigmaL_L2 : sigmaR_L2;

                double qLinear = peak_L1 + (trk_ql2 - peak_L2) * (sigma_L1 / sigma_L2);

                h2_EGELinear[key]->Fill(qLinear, ek_det[idet]);
            }
        }
    }

    // Save outputs
    TFile fout_LG(outFileName_LG.c_str(), "RECREATE");
    for (auto& kv : h2_LG) kv.second->Write();
    fout_LG.Close();
    
    TFile fout_EGE(outFileName_EGE.c_str(), "RECREATE");
    for (auto& kv : h2_EGE) kv.second->Write();
    fout_EGE.Close();
    
    TFile fout_EGELinear(outFileName_EGELinear.c_str(), "RECREATE");
    for (auto& kv : h2_EGELinear) kv.second->Write();
    fout_EGELinear.Close();

    std::cout << "All tuned L2 Q histograms saved." << std::endl;
}
