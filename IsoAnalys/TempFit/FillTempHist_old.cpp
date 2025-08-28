/***********************************************************
 *  Mass Template Generation for Isotope Analysis cpp file
 *  
 *  Author: Z.Yan
 *  Date: 2024.10.31
 *  
 *  Purpose: Generate mass templates for isotope separation
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
     float  mmom = 0;
     double innerRig = 0, richBeta = 0, TOFBeta = 0, cutoffRig = 0;
     bool richNaF = false;
     unsigned int cutStatus = 0;
 
     if (isMC) tree->SetBranchAddress("mmom", &mmom);
     tree->SetBranchAddress("InnerRig", &innerRig);
     tree->SetBranchAddress("richBeta", &richBeta);
     tree->SetBranchAddress("TOFBeta", &TOFBeta);
     tree->SetBranchAddress("cutOffRig", &cutoffRig);
     tree->SetBranchAddress("rich_NaF", &richNaF);
     tree->SetBranchAddress("cutStatus", &cutStatus);
 
     // Create histograms
     std::vector<std::unique_ptr<TH2D>> h_inv_ToFMass;
     std::vector<std::unique_ptr<TH2D>> h_inv_NaFMass;
     std::vector<std::unique_ptr<TH2D>> h_inv_AGLMass;
 
     // 获取所有mass对应的动能bins
     std::vector<std::array<double, AMS_Iso::Constants::RIGIDITY_BINS + 1>> all_ek_bins;
     std::vector<std::vector<double>> beta_bins;
         for (int mass : config.masses) {
             cout<<mass<<endl;
             const auto& ek_bins = getKineticEnergyBins(config.charge, mass);
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
         for (int j = 0; j < config.nAlpha; ++j) {
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
         }
     }
     // Process entries
     Long64_t nEntries = tree->GetEntries();
     Long64_t nProcessed = 0;
 
     for (Long64_t i = 0; i < nEntries; ++i) {
         tree->GetEntry(i);
         
         if (i % 1000000 == 0) {
             std::cout << "Processing entry " << i << "/" << nEntries 
                      << " (processed: " << nProcessed << ")" << std::endl;
         }
 
         // Calculate weight for MC
         double weight = 1.0;
         if (isMC) {
             weight = calculateWeight(mmom, config.charge, false);
             if (weight <= 0.0) continue;
         }
 
         // Process conditions using named constants
         bool processTOF = (cutStatus & TRACKER_TOF_MASK) == TRACKER_TOF_MASK && 
                          isValidBeta(TOFBeta);
         bool processNaF = richNaF && (cutStatus & TRACKER_NAF_MASK) == TRACKER_NAF_MASK && 
                          isValidBeta(richBeta);
         bool processAGL = !richNaF && (cutStatus & TRACKER_AGL_MASK) == TRACKER_AGL_MASK && 
                          isValidBeta(richBeta);
 
         // Fill histograms
         if (isMC) {
             if (processTOF) { 
                 for (int j = 0; j < config.nAlpha; ++j) {
                     double alpha = config.alphaMin + j * alphaStep;
                     auto result = calculateMass(TOFBeta, alpha, innerRig, config.charge);
                     if (result.isValid()) {
                         // 为每个mass计算动能并填充
                         for (size_t k = 0; k < config.masses.size(); ++k) {
                             int hist_idx = j * config.masses.size() + k;
                             h_inv_ToFMass[hist_idx]->Fill(result.invMass, result.ek, weight);
                         }
                     }
                 }
             }
             
             if (processNaF) {
                 double smearBeta = GetSmearRichBeta(config.charge, 7, richBeta, true);
                 if (isValidBeta(smearBeta)) {
                     for (int j = 0; j < config.nAlpha; ++j) {
                         double alpha = config.alphaMin + j * alphaStep;
                         auto result = calculateMass(smearBeta, alpha, innerRig, config.charge);
                         if (result.isValid()) {
                             for (size_t k = 0; k < config.masses.size(); ++k) {
                                 int hist_idx = j * config.masses.size() + k;
                                 h_inv_NaFMass[hist_idx]->Fill(result.invMass, result.ek, weight);
                             }
                         }
                     }
                 }
             }
             
             if (processAGL) {
                 double smearBeta = GetSmearRichBeta(config.charge, 7, richBeta, false);
                 if (isValidBeta(smearBeta)) {
                     for (int j = 0; j < config.nAlpha; ++j) {
                         double alpha = config.alphaMin + j * alphaStep;
                         auto result = calculateMass(smearBeta, alpha, innerRig, config.charge);
                         if (result.isValid()) {
                             for (size_t k = 0; k < config.masses.size(); ++k) {
                                 int hist_idx = j * config.masses.size() + k;
                                 h_inv_AGLMass[hist_idx]->Fill(result.invMass, result.ek, weight);
                             }
                         }
                     }
                 }
             }
         } else {
             // ISS数据处理
             int firstMass = config.masses[0];  // 使用第一个质量数进行cutoff判断
 
             for (size_t j = 0; j < config.masses.size(); ++j) {
                 const auto& mass_beta_bins = beta_bins[j];
                 
                 if (processTOF) {
                     int tofBin = findBin(mass_beta_bins, TOFBeta);
                     if (tofBin >= 0 && isBeyondCutoff(mass_beta_bins[tofBin], cutoffRig, 
                             SAFETY_TOF, config.charge, firstMass, false)) {
                         auto result = calculateMass(TOFBeta, 1.0, innerRig, config.charge);
                         if (result.isValid()) {
                             h_inv_ToFMass[j]->Fill(result.invMass, betaToKineticEnergy(TOFBeta));
                         }
                     }
                 }
                 
                 if (processNaF) {
                     int richBin = findBin(mass_beta_bins, richBeta);
                     if (richBin >= 0 && isBeyondCutoff(mass_beta_bins[richBin], cutoffRig, 
                             SAFETY_NAF, config.charge, firstMass, false)) {
                         auto result = calculateMass(richBeta, 1.0, innerRig, config.charge);
                         if (result.isValid()) {
                             h_inv_NaFMass[j]->Fill(result.invMass, betaToKineticEnergy(richBeta));
                         }
                     }
                 }
                 
                 if (processAGL) {
                     int richBin = findBin(mass_beta_bins, richBeta);
                     if (richBin >= 0 && isBeyondCutoff(mass_beta_bins[richBin], cutoffRig, 
                             SAFETY_AGL, config.charge, firstMass, false)) {
                         auto result = calculateMass(richBeta, 1.0, innerRig, config.charge);
                         if (result.isValid()) {
                             h_inv_AGLMass[j]->Fill(result.invMass, betaToKineticEnergy(richBeta));
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
     }
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
         // Process MC for each mass
         for (int mass : config.masses) {
             std::string mcFile = baseDir + config.mcPrefix + std::to_string(mass) + ".root";
             std::string outFile = baseDir + config.mcPrefix + std::to_string(mass) + "_temp.root";
             
             auto file = std::make_unique<TFile>(mcFile.c_str());
             if (file && !file->IsZombie()) {
                 ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), true, config, mass);
             } else {
                 std::cerr << "Failed to open MC file: " << mcFile << std::endl;
             }
         }
         */
         // Process ISS data (using only first mass in filename)
         std::string issFile = baseDir + config.dataName + std::to_string(config.masses[0]) + ".root";
         std::string outFile = baseDir + config.dataName + "_temp.root";
         
         auto file = std::make_unique<TFile>(issFile.c_str());
         if (file && !file->IsZombie()) {
             cout<<"begin iss"<<endl;
             ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), false, config);
         } else {
             std::cerr << "Failed to open ISS file: " << issFile << std::endl;
         }
         
     } catch (const std::exception& e) {
         std::cerr << "Error in BuildTempHist: " << e.what() << std::endl;
     }
 }
 
 void FillTempHist() {
     BuildTempHist("Be");  // 处理铍
 }