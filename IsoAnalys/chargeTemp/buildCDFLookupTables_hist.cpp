#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TString.h>
#include <TMath.h>
#include <TVectorD.h>
#include "../Tool.h" // Assumes Binning::NarrowBins is defined here

using namespace AMS_Iso;

// --- Global Constants ---

// PDF parameter name to array index mapping.
const std::vector<std::pair<std::string, int>> LG_paramList = {
    {"Width", 0}, {"MPV", 1}, {"Area", 2}, {"Sigma", 3}
};
const std::vector<std::pair<std::string, int>> EGE_paramList = {
    {"Peak", 0}, {"SigmaL", 1}, {"AlphaL", 2},
    {"SigmaR", 3}, {"AlphaR", 4}, {"Norm", 5}, {"xmin", 6}, {"xmax", 7}
};

// Map from atomic number (Z) to element name.
const std::map<int, std::string> elemMap = {{4, "Beryllium"}, {5, "Boron"}};
const std::vector<int> charges = {4, 5};
const std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};

// Input parameter file and output lookup table file.
const std::string histParamFileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_0p8.root";
const std::string outFileName_hist = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/CDFLookupTable_fromHist_0p8.root";


// --- Data Structures and Functions ---

// Stores the generated CDF lookup table.
struct LookupTable {
    std::vector<double> q_values, cdf_l1, cdf_l2;
    
    // Builds the CDF lookup table from two PDF functions.
    void build(TF1& f_l1, TF1& f_l2, double qmin, double qmax, int npts = 2000) {
        if (npts <= 0) {
             std::cerr << "Error: npts must be positive." << std::endl;
             return;
        }
        double dq = (qmax - qmin) / npts;
        q_values.assign(npts + 1, 0.0);
        cdf_l1.assign(npts + 1, 0.0);
        cdf_l2.assign(npts + 1, 0.0);
        
        f_l1.SetRange(qmin, qmax);
        f_l2.SetRange(qmin, qmax);
        
        double total_l1 = f_l1.Integral(qmin, qmax, 1e-9);
        double total_l2 = f_l2.Integral(qmin, qmax, 1e-9);
        
        if (total_l1 <= 0 || total_l2 <= 0) {
            std::cerr << "Warning: Integral of a PDF is zero or negative. total_l1=" << total_l1 << ", total_l2=" << total_l2 << std::endl;
        }
        
        for (int i = 0; i <= npts; ++i) {
            double q = qmin + i * dq;
            q_values[i] = q;
            double int1 = (total_l1 > 0) ? f_l1.Integral(qmin, q, 1e-9) : 0;
            double int2 = (total_l2 > 0) ? f_l2.Integral(qmin, q, 1e-9) : 0;
            cdf_l1[i] = (total_l1 > 0) ? int1 / total_l1 : 0;
            cdf_l2[i] = (total_l2 > 0) ? int2 / total_l2 : 0;
        }
    }
};

// Helper function to get a parameter from a specific histogram bin.
double getHistParam(TFile* fin, const std::string& histName, int binIndex) {
    TH1D* hist = (TH1D*)fin->Get(histName.c_str());
    if (!hist) {
        std::cerr << "Warning: Could not find histogram: " << histName << std::endl;
        return 0;
    }
    // ROOT bins are 1-indexed, so we add 1 to the 0-indexed binIndex.
    return hist->GetBinContent(binIndex + 1);
}

// Gets Landau-Gauss parameters for a given energy bin from histograms.
void getLGParamsFromHist(TFile* fin, const std::string& elem, const std::string& det, 
                         const std::string& temp, int binIndex, double* par) {
    for (const auto& p : LG_paramList) {
        std::string histName = elem + "_" + det + "_" + temp + "_LG_" + p.first;
        par[p.second] = getHistParam(fin, histName, binIndex);
        cout<<par[p.second]<<" ";
    }
    cout<<""<<endl;
    // Per user request, fix the x-scaling parameter to 1.0.
    par[4] = 1.0;
}

// Gets EGE parameters for a given energy bin from histograms.
void getEGEParamsFromHist(TFile* fin, const std::string& elem, const std::string& det, 
                          const std::string& temp, int binIndex, double* par) {
    for (const auto& p : EGE_paramList) {
        std::string histName = elem + "_" + det + "_" + temp + "_EGE_" + p.first;
        par[p.second] = getHistParam(fin, histName, binIndex);
        cout<<par[p.second]<<" ";
    }
    cout<<""<<endl;
}

// Finds the starting bin index for Ek >= 0.51 GeV.
int findStartBin() {
    for (size_t i = 0; i < Binning::NarrowBins.size(); ++i) {
        if (Binning::NarrowBins[i] >= 0.51) {
            return i;
        }
    }
    return 0;
}

// --- Main Function to Build the Tables ---

void buildCDFLookupTables_hist() {
    std::unique_ptr<TFile> finParam(TFile::Open(histParamFileName.c_str()));
    if (!finParam || finParam->IsZombie()) {
        std::cerr << "Error opening parameter file: " << histParamFileName << std::endl;
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(outFileName_hist.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error creating output file: " << outFileName_hist << std::endl;
        return;
    }
    
    int startBin = findStartBin();
    int nEk = Binning::NarrowBins.size() - 1;

    std::cout << "Building CDF lookup tables from histograms..." << std::endl;
    std::cout << "Reading from: " << histParamFileName << std::endl;
    std::cout << "Writing to: " << outFileName_hist << std::endl;
    std::cout << "Processing Ek bins from index " << startBin << " to " << nEk - 1 << std::endl;
    
    for (int z : charges) {
        std::string elemName = elemMap.at(z);
        // Set global integration range per user request.
        double qmin_global = z - 2.0;
        double qmax_global = z + 2.0;
        
        for (const std::string& det : detNames) {
            std::cout << "Procossing " << elemName << " for detector " << det << "..." << std::endl;
            
            for (int iy = startBin; iy < nEk; ++iy) {

                double ekLow = Binning::NarrowBins[iy];
                if ((det=="TOF" && ekLow>1.55) || (det=="NaF" && (ekLow<0.86 || ekLow>4.91)) || (det=="AGL" && ekLow<2.88)) continue;
                cout<<det<<" "<<iy<<endl;
                
                double par_lg_l1[5] = {0}, par_lg_l2[5] = {0};
                double par_ege_l1[8] = {0}, par_ege_l2[8] = {0};
                
                getLGParamsFromHist(finParam.get(), elemName, det, "L1Temp", iy, par_lg_l1);
                getLGParamsFromHist(finParam.get(), elemName, det, "L2Temp", iy, par_lg_l2);
                getEGEParamsFromHist(finParam.get(), elemName, det, "L1Temp", iy, par_ege_l1);
                getEGEParamsFromHist(finParam.get(), elemName, det, "L2Temp", iy, par_ege_l2);

                // Create TF1 objects with the global range.
                TF1 lg_l1("lg_l1", langaufun, qmin_global, qmax_global, 5);
                TF1 lg_l2("lg_l2", langaufun, qmin_global, qmax_global, 5);
                TF1 ege_l1("ege_l1", funcExpGausExp, qmin_global, qmax_global, 8);
                TF1 ege_l2("ege_l2", funcExpGausExp, qmin_global, qmax_global, 8);
                
                lg_l1.SetParameters(par_lg_l1); 
                lg_l2.SetParameters(par_lg_l2);
                ege_l1.SetParameters(par_ege_l1); 
                ege_l2.SetParameters(par_ege_l2);
                
                LookupTable lg_lookup, ege_lookup;
                lg_lookup.build(lg_l1, lg_l2, qmin_global, qmax_global);
                ege_lookup.build(ege_l1, ege_l2, qmin_global, qmax_global);
                
                // --- Save results to output file ---
                fout->cd();
                std::string baseName = Form("%s_%s_bin%d", elemName.c_str(), det.c_str(), iy);
                
                // Save LG lookup table
                TVectorD(lg_lookup.q_values.size(), lg_lookup.q_values.data()).Write((baseName + "_LG_q").c_str());
                TVectorD(lg_lookup.cdf_l1.size(), lg_lookup.cdf_l1.data()).Write((baseName + "_LG_cdf_l1").c_str());
                TVectorD(lg_lookup.cdf_l2.size(), lg_lookup.cdf_l2.data()).Write((baseName + "_LG_cdf_l2").c_str());
                
                // Save EGE lookup table
                TVectorD(ege_lookup.q_values.size(), ege_lookup.q_values.data()).Write((baseName + "_EGE_q").c_str());
                TVectorD(ege_lookup.cdf_l1.size(), ege_lookup.cdf_l1.data()).Write((baseName + "_EGE_cdf_l1").c_str());
                TVectorD(ege_lookup.cdf_l2.size(), ege_lookup.cdf_l2.data()).Write((baseName + "_EGE_cdf_l2").c_str());
                
                // Save LG parameters for potential linear mapping later
                TVectorD(5, par_lg_l1).Write((baseName + "_LG_params_l1").c_str());
                TVectorD(5, par_lg_l2).Write((baseName + "_LG_params_l2").c_str());
                // Save EGE parameters for potential linear mapping later
                TVectorD(8, par_ege_l1).Write((baseName + "_EGE_params_l1").c_str());
                TVectorD(8, par_ege_l2).Write((baseName + "_EGE_params_l2").c_str());
            }
        }
    }
    
    fout->Close();
    std::cout << "CDF lookup tables saved to: " << outFileName_hist << std::endl;
}
