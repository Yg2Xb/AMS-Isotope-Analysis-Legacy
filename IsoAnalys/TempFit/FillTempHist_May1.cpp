/***********************************************************
 *  Mass Template Generation for Isotope Analysis cpp file
 *  
 *  Author: Z.Yan
 *  Date: 2024.10.31
 *  
 *  Purpose: select pure be10, validate iss and mc consistancy
 ***********************************************************/

#include "FillTempHist.h"
using namespace AMS_Iso;


void ProcessTreeWithMass(TFile* file, const char* treeName, const char* outputFileName, 
                        bool isMC, const IsotopeConfig& config, int processMass = 0) {
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Invalid input file!" << std::endl;
        return;
    }

    auto* tree = static_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        std::cerr << "Error: Cannot find tree " << treeName << std::endl;
        return;
    }

    // Setup branches
    double InnerRig = -1, L1InnerRig = -1, TOFBeta = -1, NaFBeta = -1, AGLBeta = -1, TOFEk = -1, NaFEk = -1, AGLEk = -1;
    bool rich_NaF = false;
    unsigned int cutStatus = 0;
    int ntrack = 0;
    
        double generatedRig = -1;
        int mtrpar[9];
    if (isMC) 
    {
        tree->SetBranchAddress("generatedRig", &generatedRig);
        tree->SetBranchAddress("mtrpar", &mtrpar);
    }
    float cutoffpi[2] = {-999, -999}, iso_ll[3][3], mcutoffi[4][2], tk_qin[2][3], tk_exqln[2][2][3];
    double cutOffRig = -1;
    int btstat_new = -1;
    if(!isMC)
    {
        tree->SetBranchAddress("cutoffpi", &cutoffpi);
        tree->SetBranchAddress("iso_ll", &iso_ll);
        tree->SetBranchAddress("cutOffRig", &cutOffRig);
        tree->SetBranchAddress("btstat_new", &btstat_new);
        tree->SetBranchAddress("tk_qin", &tk_qin);
        tree->SetBranchAddress("tk_exqln", &tk_exqln);
    }
        
    tree->SetBranchAddress("InnerRig", &InnerRig);
    tree->SetBranchAddress("NaFBeta", &NaFBeta);
    tree->SetBranchAddress("AGLBeta", &AGLBeta);
    tree->SetBranchAddress("TOFBeta", &TOFBeta);
    tree->SetBranchAddress("TOFEk", &TOFEk);
    tree->SetBranchAddress("AGLEk", &AGLEk);
    tree->SetBranchAddress("NaFEk", &NaFEk);
    tree->SetBranchAddress("rich_NaF", &rich_NaF);
    tree->SetBranchAddress("cutStatus", &cutStatus);
    tree->SetBranchAddress("ntrack", &ntrack);

    // Create histograms
    std::vector<std::unique_ptr<TH2D>> h_inv_ToFMass;
    std::vector<std::unique_ptr<TH2D>> h_inv_NaFMass;
    std::vector<std::unique_ptr<TH2D>> h_inv_AGLMass;
    std::vector<std::unique_ptr<TH2D>> h_inv_ToFEst;
    std::vector<std::unique_ptr<TH2D>> h_inv_NaFEst;
    std::vector<std::unique_ptr<TH2D>> h_inv_AGLEst;

    // 获取所有mass对应的动能bins
    std::vector<std::array<double, AMS_Iso::Constants::RIGIDITY_BINS + 1>> all_ek_bins;
    //std::vector<std::array<double, 15>> all_ek_bins;
    std::vector<std::vector<double>> beta_bins;
        for (int mass : config.masses) {
            cout<<mass<<endl;
            const auto& ek_bins = getKineticEnergyBins(config.charge, mass);//Binning::WideBins;//
            cout<<ek_bins[10]<<endl;
            std::vector<double> Rbins_beta;
            for (double ek : ek_bins) {
                double beta = kineticEnergyToBeta(ek);
                if (isValidBeta(beta)) {
                    Rbins_beta.push_back(beta);
                }
            }
            beta_bins.push_back(Rbins_beta);
            all_ek_bins.push_back(ek_bins);
        }

    // Safety factors
    constexpr double SAFETY_TOF = 1.06;
    constexpr double SAFETY_NAF = 1.005;
    constexpr double SAFETY_AGL = 1.0005;
    double alphaStep = (config.alphaMax - config.alphaMin) / (config.nAlpha - 1);
    if (isMC) {
        // MC模式：为每个alpha和每个mass创建直方图
        for (int j = 0; j < config.nAlpha ; ++j) {
            for (size_t k = 0; k < config.masses.size(); ++k) {
                int mass = config.masses[k];
                const auto& ek_bins = all_ek_bins[k];
                if(j==0) cout<<ek_bins[10]<<endl;
                
                // ToF直方图
                h_inv_ToFMass.push_back(std::make_unique<TH2D>(
                    Form("h_inv_ToFMass_mass%d_bin%d_alpha%d", processMass, mass, j),
                    "inv_ToFMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));
                
                // NaF直方图
                h_inv_NaFMass.push_back(std::make_unique<TH2D>(
                    Form("h_inv_NaFMass_mass%d_bin%d_alpha%d", processMass, mass, j),
                    "inv_NaFMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));
                
                // AGL直方图
                h_inv_AGLMass.push_back(std::make_unique<TH2D>(
                    Form("h_inv_AGLMass_mass%d_bin%d_alpha%d", processMass, mass, j),
                    "inv_AGLMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
                ));
            }
        }
    }  else {
        // ISS模式：为每个mass创建一个直方图
        for (size_t j = 0; j < config.masses.size(); ++j) {
            cout<<"begin iss hist"<<endl;
            const auto& ek_bins = all_ek_bins[j];
            cout<<ek_bins[10]<<endl;
            
            h_inv_ToFMass.push_back(std::make_unique<TH2D>(
                Form("h_inv_ToFMass_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_ToFMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
            ));
            
            h_inv_NaFMass.push_back(std::make_unique<TH2D>(
                Form("h_inv_NaFMass_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_NaFMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
            ));
            
            h_inv_AGLMass.push_back(std::make_unique<TH2D>(
                Form("h_inv_AGLMass_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_AGLMass", 200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data()
            ));

            h_inv_ToFEst.push_back(std::make_unique<TH2D>(
                Form("h_inv_ToFEst_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_ToFMass", 500, 0.0, 1.0, ek_bins.size()-1, ek_bins.data()
            ));
            
            h_inv_NaFEst.push_back(std::make_unique<TH2D>(
                Form("h_inv_NaFEst_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_NaFMass", 500, 0.0, 1.0, ek_bins.size()-1, ek_bins.data()
            ));
            
            h_inv_AGLEst.push_back(std::make_unique<TH2D>(
                Form("h_inv_AGLEst_%sUseMass%d", config.dataName.c_str(), config.masses[j]),
                "inv_AGLMass", 500, 0.0, 1.0, ek_bins.size()-1, ek_bins.data()
            ));
        }
    }
    //hist for pure selection study
    TH2F *hToF = new TH2F("hToF","Pass ToF Selection;cutoff rigidity[GV];ToF #beta", 250, 0, 8, 250, 0.6, 0.9);
    TH2F *hNaF = new TH2F("hNaF","Pass NaF Selection;cutoff rigidity[GV];NaF #beta", 250, 0, 20, 250, 0.75, 1.001);
    TH2F *hAGL = new TH2F("hAGL","Pass AGL Selection;cutoff rigidity[GV];AGL #beta", 250, 4, 30, 250, 0.954, 1.001);
    
    TH2F *hToF_pure = new TH2F("hToF_pure","Pass ToF and Pure Be10 Selection;cutoff rigidity[GV];ToF #beta", 250, 0, 8, 250, 0.6, 0.9);
    TH2F *hNaF_pure = new TH2F("hNaF_pure","Pass NaF and Pure Be10 Selection;cutoff rigidity[GV];NaF #beta", 250, 0, 20, 250, 0.75, 1.001);
    TH2F *hAGL_pure = new TH2F("hAGL_pure","Pass AGL and Pure Be10 Selection;cutoff rigidity[GV];AGL #beta", 250, 4, 30, 250, 0.954, 1.001);


    std::unique_ptr<TFile> fluxFile(TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BerylliumFlux.root"));
    if(!fluxFile || !fluxFile->Get("graph1")) throw std::runtime_error("Cannot access flux file/graph");

    TF1 *sourceFlux_ = SplineFit(static_cast<TGraph*>(fluxFile->Get("graph1")), 
            kSourceFluxFitPoints.data(), 
            kSourceFluxFitPoints.size(), 
            0x38, "b1e1", "f_SourceFlux", 1., 2000.0);

    int ia = processMass == 7 ? 0 : (processMass == 9 ? 1 : 2);
    // Process entries
    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessed = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        //if(i > 10000) break;
        
        if (i % 1000000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries 
                     << " (processed: " << nProcessed << ")" << std::endl;
        }

        // Calculate weight for MC
        double weight = 1.0;
        if (isMC) {
            if (generatedRig <= 0) continue;
            double flux_weight = sourceFlux_->Eval(generatedRig) / sourceFlux_->Integral(1., 2000.);
            double mc_weight = getMCFunction().Eval(generatedRig) / getMCNorm();
            weight = flux_weight / mc_weight;
            if (weight <= 0.0) continue;
        }
        // Process conditions using named constants
        bool processTOF = (cutStatus & TRACKER_TOF_MASK) == TRACKER_TOF_MASK && 
                         isValidBeta(TOFBeta);
        bool processNaF = (cutStatus & TRACKER_RICH_MASK) == TRACKER_RICH_MASK && 
                         isValidBeta(NaFBeta);
        bool processAGL = (cutStatus & TRACKER_RICH_MASK) == TRACKER_RICH_MASK && 
                         isValidBeta(AGLBeta);

        // Fill histograms
        int isotopeCode = (processMass == 7) ? 63 : (processMass == 9) ? 64 : (processMass == 10) ? 114 : -1;  // Be7=63, Be9=64, Be10=114
        if (isMC) {
            //select below L2 pure Isotope
            bool pureMC = mtrpar[1] == isotopeCode;
            if(!pureMC) continue;

            if (processTOF) { 
                for (int j = 0; j < config.nAlpha; ++j) {
                    double alpha = config.alphaMin + j * alphaStep;
                    auto result = calculateMass(TOFBeta, alpha, InnerRig, config.charge);
                    if (result.isValid()) {
                        for (size_t k = 0; k < config.masses.size(); ++k) {
                            int hist_idx = j * config.masses.size() + k;
                            h_inv_ToFMass[hist_idx]->Fill(result.invMass, TOFEk, weight);
                        }
                    }
                }
            }
            if (processNaF) {
                double smearBeta = GetSmearRichBeta(config.charge, ia, NaFBeta, true);
                if (isValidBeta(smearBeta)) {
                    for (int j = 0; j < config.nAlpha; ++j) {
                        double alpha = config.alphaMin + j * alphaStep;
                        auto result = calculateMass(smearBeta, alpha, InnerRig, config.charge);
                        if (result.isValid()) {
                            for (size_t k = 0; k < config.masses.size(); ++k) {
                                int hist_idx = j * config.masses.size() + k;
                                h_inv_NaFMass[hist_idx]->Fill(result.invMass, NaFEk, weight);
                            }
                        }
                    }
                }
            }
            if (processAGL) {
                double smearBeta = GetSmearRichBeta(config.charge, ia, AGLBeta, false);
                if (isValidBeta(smearBeta)) {
                    for (int j = 0; j < config.nAlpha; ++j) {
                        double alpha = config.alphaMin + j * alphaStep;
                        auto result = calculateMass(smearBeta, alpha, InnerRig, config.charge);
                        if (result.isValid()) {
                            for (size_t k = 0; k < config.masses.size(); ++k) {
                                int hist_idx = j * config.masses.size() + k;
                                h_inv_AGLMass[hist_idx]->Fill(result.invMass, AGLEk, weight);
                            }
                        }
                    }
                }
            }
        } else {
            // ISS数据处理
            
            //if(InnerRig <= cutOffRig) continue;
            
            //double q_l1_unbiased = tk_exqln[0][0][2];
            //bool Qcut = (tk_qin[0][2] > 3.5 && tk_qin[0][2] < 4.5) && (q_l1_unbiased > 4.5 && q_l1_unbiased < 5.5);
            //if(!Qcut) continue;
            
            for (size_t j = 0; j < config.masses.size(); ++j) {
                const auto& mass_beta_bins = beta_bins[j];
                int TOFBin = findBin(mass_beta_bins, TOFBeta);
                int NaFBin = findBin(mass_beta_bins, NaFBeta);
                int AGLBin = findBin(mass_beta_bins, AGLBeta);
                
                if (processTOF) {
                    //if(j==0){hToF->Fill(cutoffpi[1], TOFBeta);}
                    //if (selectPure(TOFBeta, cutoffpi[1], SAFETY_TOF, config.charge, InnerRig, false)) {
                    //if(Qcut){
                    if(isBeyondCutoff(mass_beta_bins[TOFBin], cutOffRig, SAFETY_TOF, config.charge, config.masses[j], !isMC)){
                        //cout<<InnerRig<<" "<<tk_qin[0][2]<<" "<<q_l1_unbiased<<" "<<TOFEk<<endl;
                        auto result = calculateMass(TOFBeta, 1.0, InnerRig, config.charge);
                        //if(j==0){hToF_pure->Fill(cutoffpi[1], TOFBeta);}
                        if (result.isValid()) {
                            h_inv_ToFMass[j]->Fill(result.invMass, TOFEk);
                        }
                        if(iso_ll[0][0] >= 0 && iso_ll[1][0] >= 0 && iso_ll[2][0] >= 0)
                        {
                            h_inv_ToFEst[j]->Fill(iso_ll[2][0]/(iso_ll[0][0] + iso_ll[1][0] + iso_ll[2][0]), TOFEk);
                        }
                    }
                }
                
                if (processNaF) {
                    //if(j==0){hNaF->Fill(cutoffpi[1], NaFBeta);}
                    //if (selectPure(NaFBeta, cutoffpi[1], SAFETY_NAF, config.charge, InnerRig, false)) {
                    //if(Qcut){
                    if(isBeyondCutoff(mass_beta_bins[NaFBin], cutOffRig, SAFETY_NAF, config.charge, config.masses[j], !isMC)){
                        auto result = calculateMass(NaFBeta, 1.0, InnerRig, config.charge);
                        //if(j==0){hNaF_pure->Fill(cutoffpi[1], NaFBeta);}
                        if (result.isValid()) {
                            h_inv_NaFMass[j]->Fill(result.invMass, NaFEk);
                        }
                        if(iso_ll[0][1] >= 0 && iso_ll[1][1] >= 0 && iso_ll[2][1] >= 0)
                        {
                            h_inv_NaFEst[j]->Fill(iso_ll[2][1]/(iso_ll[0][1] + iso_ll[1][1] + iso_ll[2][1]), NaFEk);
                        }
                    }
                }
                
                if (processAGL) {
                    //if(j==0){hAGL->Fill(cutoffpi[1], AGLBeta);}
                    //if (selectPure(AGLBeta, cutoffpi[1], SAFETY_AGL, config.charge InnerRig, false)) {
                    //if(Qcut){
                    if(isBeyondCutoff(mass_beta_bins[AGLBin], cutOffRig, SAFETY_AGL, config.charge, config.masses[j], !isMC)){
                        auto result = calculateMass(AGLBeta, 1.0, InnerRig, config.charge);
                        //if(j==0){hAGL_pure->Fill(cutoffpi[1], AGLBeta);}
                        if (result.isValid()) {
                            h_inv_AGLMass[j]->Fill(result.invMass, AGLEk);
                        }
                        if(iso_ll[0][2] >= 0 && iso_ll[1][2] >= 0 && iso_ll[2][2] >= 0)
                        {
                            h_inv_AGLEst[j]->Fill(iso_ll[2][2]/(iso_ll[0][2] + iso_ll[1][2] + iso_ll[2][2]), AGLEk);
                        }
                    }
                }
            }
        }
        
        nProcessed++;
    }



    std::cout << "Processed " << nProcessed << "/" << nEntries << " entries" << std::endl;

    // Save histograms
    // Save histograms
    auto outputFile = std::make_unique<TFile>(outputFileName, "RECREATE");
    for (size_t i = 0; i < h_inv_ToFMass.size(); ++i) {
        h_inv_ToFMass[i]->Write();
        h_inv_NaFMass[i]->Write();
        h_inv_AGLMass[i]->Write();
        if(!isMC){
        h_inv_ToFEst[i]->Write();
        h_inv_NaFEst[i]->Write();
        h_inv_AGLEst[i]->Write();
        }
    }
    hToF->Write();
    hNaF->Write();
    hAGL->Write();
    hToF_pure->Write();
    hNaF_pure->Write();
    hAGL_pure->Write();
    outputFile->Close();

    // ToF画布
    TCanvas *c1 = new TCanvas("c1","ToF Beta vs Cutoff", 1200,900);
    c1->SetRightMargin(0.18);
    c1->SetLogz();
    hToF->Draw("colz");
    DrawTheoryLines(c1, 0.6, 0.9);
    //c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/ToF_beta_cutoff.png");
    hToF_pure->GetZaxis()->SetRangeUser(1,hToF->GetMaximum());
    hToF_pure->Draw("colz");
    DrawTheoryLines(c1, 0.6, 0.9);
    //c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/ToF_beta_cutoff_pure.png");

    // NaF画布
    TCanvas *c2 = new TCanvas("c2","NaF Beta vs Cutoff", 1200,900);
    c2->SetRightMargin(0.18);
    c2->SetLogz();
    hNaF->Draw("colz");
    DrawTheoryLines(c2, 0.75, 1.005);
    //c2->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/NaF_beta_cutoff.png");
    hNaF_pure->GetZaxis()->SetRangeUser(1,hNaF->GetMaximum());
    hNaF_pure->Draw("colz");
    DrawTheoryLines(c2, 0.75, 1.005);
    //c2->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/NaF_beta_cutoff_pure.png");

    // AGL画布
    TCanvas *c3 = new TCanvas("c3","AGL Beta vs Cutoff", 1200,900);
    c3->SetRightMargin(0.18);
    c3->SetLogz();
    hAGL->Draw("colz");
    DrawTheoryLines(c3, 0.95, 1.005);
    //c3->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/AGL_beta_cutoff.png");
    hAGL_pure->GetZaxis()->SetRangeUser(1,hAGL->GetMaximum());
    hAGL_pure->Draw("colz");
    DrawTheoryLines(c3, 0.95, 1.005);
    //c3->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/AGL_beta_cutoff_pure.png");

    // 清理内存
    delete c1;
    delete c2;
    delete c3;

}

void BuildTempHist(const std::string& isotype) {
    const std::string baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/";
    
    try {
        if (isotope_configs.find(isotype) == isotope_configs.end()) {
            std::cerr << "Unknown isotope type: " << isotype << std::endl;
            return;
        }

        const auto& config = isotope_configs.at(isotype);

        // Process ISS data (using only first mass in filename)
        std::string issFile = baseDir + config.dataName + std::to_string(config.masses[0]) + ".root";
        //std::string issFile = baseDir + config.dataName + "andBforBkg" + ".root";
        //std::string outFile = baseDir + config.dataName + "_temp_bkgBe.root";
        std::string outFile = baseDir + config.dataName + "_temp.root";
        
        auto file = std::make_unique<TFile>(issFile.c_str());
        if (file && !file->IsZombie()) {
            cout<<"begin iss"<<endl;
            ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), false, config);
        } else {
            std::cerr << "Failed to open ISS file: " << issFile << std::endl;
        }

        // Process MC for each mass
        for (int mass : config.masses) {
            std::string mcFile = baseDir + config.mcPrefix + std::to_string(mass) + ".root";
            std::string outFile = baseDir + config.mcPrefix + std::to_string(mass) + "_temp_bkg.root";
            
            auto file = std::make_unique<TFile>(mcFile.c_str());
            if (file && !file->IsZombie()) {
                ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), true, config, mass);
            } else {
                std::cerr << "Failed to open MC file: " << mcFile << std::endl;
            }
        }
        

        
    } catch (const std::exception& e) {
        std::cerr << "Error in BuildTempHist: " << e.what() << std::endl;
    }
}

void FillTempHist() {
    BuildTempHist("Be");  // 处理铍
}