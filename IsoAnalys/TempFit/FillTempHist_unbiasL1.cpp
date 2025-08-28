/***********************************************************
 *  Mass Template Generation for Isotope Analysis cpp file
 *  
 *  Author: Z.Yan
 *  Date: 2024.10.31
 *  Update Date: 2025.5.13 Happy Birth Day My Sweety
 *  
 *  Purpose1: Generate mass templates for isotope separation
 *  Purpose2: select pure be10, validate iss and mc consistancy
 *  Update: Backgroud Study
 ***********************************************************/

#include "FillTempHist.h"
using namespace AMS_Iso;

void ProcessTreeWithMass(TFile* file, const char* treeName, const char* outputFileName, 
                        bool isMC, const IsotopeConfig& config, int processMass) {
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Invalid input file!" << std::endl;
        return;
    }

    auto* tree = static_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        std::cerr << "Error: Cannot find tree " << treeName << std::endl;
        return;
    }
    
    //Setup used charge
    double charge = 4.0; //change!! config.charge;
    //double charge = config.charge;
    int GeneID = getIsotopeID(config, processMass);
    cout<<GeneID<<endl;

    // Setup branches
    double InnerRig = -1, TOFBeta = -1, NaFBeta = -1, AGLBeta = -1, TOFEk = -1, NaFEk = -1, AGLEk = -1;
    float tk_qin[2][3], tk_exqln[2][2][3], tk_qln[2][9][3], rich_q[2], tof_ql[4];
    bool rich_NaF = false;
    unsigned int cutStatus = 0;
    
    float mmom; double generatedRig = -1; int mtrpar[9];
    // MC branch
    if (isMC) {
        tree->SetBranchAddress("mmom", &mmom);
        tree->SetBranchAddress("generatedRig", &generatedRig);
        tree->SetBranchAddress("mtrpar", mtrpar);
    }
    
    // ISS branch
    float cutoffpi[2] = {-999, -999}, iso_ll[3][3], mcutoffi[4][2];
    double cutOffRig = -1; 
    int btstat_new = -1;
    if(!isMC) {
        tree->SetBranchAddress("cutoffpi", cutoffpi);
        tree->SetBranchAddress("iso_ll", iso_ll);
        tree->SetBranchAddress("cutOffRig", &cutOffRig);
        tree->SetBranchAddress("btstat_new", &btstat_new);
    }
    
    // common branch
    tree->SetBranchAddress("tk_qin", tk_qin);
    tree->SetBranchAddress("rich_q", rich_q);
    tree->SetBranchAddress("tof_ql", tof_ql);
    tree->SetBranchAddress("tk_exqln", tk_exqln);
    tree->SetBranchAddress("tk_qln", tk_qln);
    tree->SetBranchAddress("InnerRig", &InnerRig);
    tree->SetBranchAddress("NaFBeta", &NaFBeta);
    tree->SetBranchAddress("AGLBeta", &AGLBeta);
    tree->SetBranchAddress("TOFBeta", &TOFBeta);
    tree->SetBranchAddress("TOFEk", &TOFEk);
    tree->SetBranchAddress("AGLEk", &AGLEk);
    tree->SetBranchAddress("NaFEk", &NaFEk);
    tree->SetBranchAddress("rich_NaF", &rich_NaF);
    tree->SetBranchAddress("cutStatus", &cutStatus);

    //bin initial
    std::vector<std::array<double, 15>> all_ek_bins;
    std::vector<std::vector<double>> beta_bins;
    for (int mass : config.masses) {
        std::cout << mass << std::endl;
        const auto& ek_bins = Binning::WideBins;
        std::cout << ek_bins[10] << std::endl;
        
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
   
    //safety factor
    constexpr double SAFETY_FACTORS[DETECTOR_COUNT] = {1.06, 1.005, 1.0005}; // TOF, NaF, AGL
    //alpha tuning
    double alphaStep = (config.alphaMax - config.alphaMin) / (config.nAlpha - 1);
    
    HistogramSet histSet;
    ValidHistograms vHist;
    BackgroundHistograms bgHist;
    
    //hist initial
    if (isMC) {
        for (int j = 10; j < 11; ++j) {
            histSet.createMCHistograms(processMass, config.masses, j, all_ek_bins);
        }
    } else {
        histSet.createISSHistograms(config.dataName, config.masses, all_ek_bins);
    }

    // flux for mc reweight-----------
    std::unique_ptr<TFile> fluxFile(TFile::Open(
        (config.charge == 4) ? 
        "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BerylliumFlux.root" :
        "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BoronFlux.root"
    ));
    
    if(!fluxFile || !fluxFile->Get("graph1")) {
        throw std::runtime_error("Cannot access flux file/graph");
    }
    //smooth
    TF1 *sourceFlux_ = SplineFit(static_cast<TGraph*>(fluxFile->Get("graph1")), 
            kSourceFluxFitPoints.data(), 
            kSourceFluxFitPoints.size(), 
            0x38, "b1e1", "f_SourceFlux", 1., 2000.0);
    //----------------------------------

    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessed = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        //if(i>20000)break;

        if (i % 1000000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries 
                     << " (processed: " << nProcessed << ")" << std::endl;
        }
        
        bool basiccut = (cutStatus & BKG_BASIC_MASK) == BKG_BASIC_MASK; 
        //bool basiccut = (cutStatus & 0x01) == 0x01; 
        if(!basiccut) continue;
        
        double q_inner = tk_qin[0][2];
        double q_l1 = tk_qln[0][0][2];
        double q_l1_unbiased = tk_exqln[0][0][2];
        double q_rich = sqrt(rich_q[0]);
        
        bool Qcut = (config.charge == -1) ? 
                    (q_inner > 4.5 && q_inner < 5.5 && q_l1_unbiased > 4.5 && q_l1_unbiased < 5.5) 
                    : (true ? (q_inner > 3.5 && q_inner < 5.5 && q_l1_unbiased > 4.8 && q_l1_unbiased < 5.5)//build be mc template q_inner > 3.5 && q_inner < 4.5
                           : (q_inner > 3.45 && q_inner < 5.45 && q_l1_unbiased > 5.0 && q_l1_unbiased < 5.4)); 
        //select L1 Boron(with Be C contamination) + Inner Be, B frag. to Be + Be contami. sample
        if(!Qcut) continue;

        bool TargetQcut = q_inner > 3.45 && q_inner < 4.45;//q_inner > 3.5 && q_inner < 4.5;
        
        // mc reweight according to flux
        double weight = 1.0;
        if (isMC) {
            if (generatedRig <= 0) continue;
            double flux_weight = sourceFlux_->Eval(generatedRig) / sourceFlux_->Integral(1., 2000.);
            double mc_weight = getMCFunction().Eval(generatedRig) / getMCNorm();
            weight = flux_weight / mc_weight;
            if (weight <= 0.0) continue;
        }
        weight = 1.0; // No Weight change!!
        
        bool detectorValid[DETECTOR_COUNT] = {
            //<0 if no det. beta
            (cutStatus & BKG_TOF_MASK) == BKG_TOF_MASK && isValidBeta(TOFBeta),
            (cutStatus & BKG_RICH_MASK) == BKG_RICH_MASK && isValidBeta(NaFBeta),
            (cutStatus & BKG_RICH_MASK) == BKG_RICH_MASK && isValidBeta(AGLBeta)
        };

        std::array<double, DETECTOR_COUNT> Beta = {TOFBeta, NaFBeta, AGLBeta};
        std::array<double, DETECTOR_COUNT> EkperNuc = {TOFEk, NaFEk, AGLEk};
        std::array<int, 3> BeIsoID = {63, 64, 114};
        
        if (isMC) {
            bool origMC_L1 = (mtrpar[0] == GeneID); // maintain generated particle in L1
            bool origMC_L2 = (mtrpar[1] == GeneID); // maintain generated particle in L2
            
            for (int det = 0; det < DETECTOR_COUNT; det++) {
                if (!detectorValid[det]) continue;
                
                DetectorMeasure mc = DetectorMeasure::get(
                    static_cast<DetectorType>(det), Beta, EkperNuc
                );
                
                //change!!
                //25.5.13 study frag in mc, since we think event fragment between L1L2, almost above tof, so in beta(mass) detect range, event are already Beryllium, so tuning them as be
                bool TempFitTargetQcut = (det == TOF) ? TargetQcut : TargetQcut && q_rich > 3;
                //same as sdiat, require up tof q
                //TempFitTargetQcut = TempFitTargetQcut && (cutStatus & (1<<13)) && q_rich<6;
                // mc smear for rich
                double TunedBeta = (det == TOF) ? mc.beta : GetSmearRichBeta(charge, mc.beta, det == NaF);
                TunedBeta = mc.beta; //No Smear change!! 
                double TunedEk = betaToKineticEnergy(TunedBeta);
               
                // backgroud study: Qcut sample counts(den)
                bgHist.SourceCounts[det][1]->Fill(TunedEk, weight); 
                
                //require B11 in L1 (no TOI influence)
                if(!origMC_L1) continue;
                //---------------!!!---------------change?
                
                bgHist.SourceCounts[det][0]->Fill(TunedEk, weight);
                // pass ber inner q
                if(TargetQcut) {bgHist.FragNucCounts[det]->Fill(TunedEk, weight);}
                // pass ber inner q cut and confirm as ber iso in L2
                for(int iso = 0; iso < 3; iso++){
                    if (TargetQcut && mtrpar[1] == BeIsoID[iso]){
                        bgHist.SourceToIsoCounts[det][iso]->Fill(TunedEk, weight);
                    }
                }
                
                if (TempFitTargetQcut && isValidBeta(TunedBeta)) {
                    for (int j = 10; j < 11; ++j) {
                        double alpha = config.alphaMin + j * alphaStep;
                        auto result = calculateMass(TunedBeta, alpha, InnerRig, charge); //change!!
                        
                        if (result.isValid()) {
                            
                            bgHist.FragNucID_Ek[det]->Fill(mtrpar[1], result.ek, weight);
                            
                            for (size_t k = 0; k < config.masses.size(); ++k) {
                                int hist_idx = k; //if nornal alphat tuing, j * config.masses.size() + k;
                                //if (result.ek != TunedEk) {
                                    //cout << "Wrong Ek!!! " << result.ek << " " << TunedEk  << " Diff: " << fabs(result.ek - TunedEk) << endl;
                                //}
                                histSet.massHist[det][hist_idx]->Fill(result.invMass, result.ek, weight);
                                if (origMC_L2) {
                                    histSet.origMassHist[det][hist_idx]->Fill(result.invMass, result.ek, weight);
                                }
                            }
                            
                            //frag be7 9 10
                            for (size_t k = 0; k < 3; ++k) {
                                int hist_idx = k; //if nornal alphat tuing, j * config.masses.size() + k;
                                if (mtrpar[1] == BeIsoID[k]) {
                                    histSet.pureFragNucMassHist[det][hist_idx]->Fill(result.invMass, result.ek, weight);
                                }
                            }

                        }
                    }
                }
            }
        } else {
            for (size_t j = 0; j < 1; ++j) { //config.masses.size
                const auto& mass_beta_bins = beta_bins[j];
                //cout<<"mass_beta_bins1:"<<mass_beta_bins[1]<<endl;
                
                for (int det = 0; det < DETECTOR_COUNT; det++) {
                    if (!detectorValid[det]) continue;
                    
                    bool TempFitTargetQcut = (det == TOF) ? TargetQcut : TargetQcut && q_rich > 3;
                    //same as sdiat, require up tof q
                    //TempFitTargetQcut = TempFitTargetQcut && (cutStatus & (1<<13)) && q_rich<6;;
                    
                    DetectorMeasure data = DetectorMeasure::get(
                        static_cast<DetectorType>(det), Beta, EkperNuc
                    );
                    
                    int betaBin = findBin(mass_beta_bins, data.beta);
                    
                    if (j == 0) {
                        //vHist.CutoffHist[det]->Fill(cutoffpi[1], data.beta);
                        //vHist.pureCutoffHist[det]->Fill(cutoffpi[1], data.beta);
                    }

                    bool cutoffcut = (betaBin >= 0) ? isBeyondCutoff(mass_beta_bins[betaBin], cutOffRig, SAFETY_FACTORS[det], 5, 10, isMC) : false;
                    
                    if (cutoffcut) {
                        //no truth info in iss, no source[det][0]
                        bgHist.SourceCounts[det][1]->Fill(data.ek);//total counts pass Qcut and cutoff cut

                        if(!TargetQcut) continue;
                        bgHist.FragNucCounts[det]->Fill(data.ek);//counts pass inner frag q cut

                        if(!TempFitTargetQcut) continue; // extra cut for mass tempfit
                        auto result = calculateMass(data.beta, 1.0, InnerRig, charge);
                        
                        if (result.isValid()) {
                            histSet.massHist[det][j]->Fill(result.invMass, data.ek);
                        }
                        
                        if (iso_ll[0][det] >= 0 && iso_ll[1][det] >= 0 && iso_ll[2][det] >= 0) {
                            histSet.estHist[det][j]->Fill(
                                iso_ll[2][det] / (iso_ll[0][det] + iso_ll[1][det] + iso_ll[2][det]), 
                                data.ek
                            );
                        }
                    }
                }
            }
        }
        
        nProcessed++;
    }

    if(isMC) {
        std::cout<<"calculate FragIso/Source Ratios"<<std::endl;
        bgHist.calculateRatios();
    }

    std::cout << "Processed " << nProcessed << "/" << nEntries << " entries" << std::endl;

    auto outputFile = std::make_unique<TFile>(outputFileName, "RECREATE");
    
    histSet.writeToFile(outputFile.get(), isMC);
    bgHist.writeToFile(outputFile.get());
    //vHist.writeToFile(outputFile.get());
    //vHist.drawAndSave();
 
    outputFile->Close();
}

void BuildTempHist(const std::string& isotype) {
    const std::string baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/";
    
    try {
        if (isotope_configs.find(isotype) == isotope_configs.end()) {
            std::cerr << "Unknown isotope type: " << isotype << std::endl;
            return;
        }

        const auto& config = isotope_configs.at(isotype);
        /* 
        // process iss
        {
            //std::string issFile = baseDir + config.dataName + std::to_string(config.masses[0]) + ".root";
            std::string issFile = baseDir + "BerandBforBkg" + ".root";
            std::string outFile = baseDir + config.dataName + "_temp_wide_bkg_4.8to5.5_NoW_NoS.root";
            
            auto file = std::make_unique<TFile>(issFile.c_str());
            if (file && !file->IsZombie()) {
                cout << "Begin ISS processing..." << endl;
                ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), false, config);
                file->Close(); 
            } else {
                std::cerr << "Failed to open ISS file: " << issFile << std::endl;
            }
        }
            */
        // process each mc
        for (int mass : config.masses) {
            std::string mcFile = baseDir + config.mcPrefix + std::to_string(mass) + "_bkg.root";
            std::string outFile = baseDir + config.mcPrefix + std::to_string(mass) + "_temp_wide_bkg_frag_4.8to5.5_NoW_NoS.root";
            
            auto file = std::make_unique<TFile>(mcFile.c_str());
            if (file && !file->IsZombie()) {
                cout << "Begin MC processing for mass " << mass << "..." << endl;
                ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), true, config, mass);
                file->Close(); 
            } else {
                std::cerr << "Failed to open MC file: " << mcFile << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in BuildTempHist: " << e.what() << std::endl;
    }
}

void FillTempHist(const std::string& NucName) {
    BuildTempHist(NucName);  
}