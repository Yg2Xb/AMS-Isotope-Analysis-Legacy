/*********************************************************** 
 * Analysis of background events and acceptance from MC fragmentation.
 * Processes MC simulation data for various nuclei and calculates their acceptances.
 * Author: Yan Zixuan  Date: 20250309
 ***********************************************************/
 #include "../Tool.h"
 #include <TFile.h> 
 #include <TTree.h>
 #include <TGraphAsymmErrors.h>
 #include <memory>
 #include <array>
 #include <vector>
 #include <stdexcept>
 
 using namespace AMS_Iso;
 
 // Structure to hold nucleus information
 struct NucleusInfo {
	 std::string name;      // Full name (e.g., "Beryllium")
	 int Z;                 // Atomic number
	 std::vector<int> A;    // Mass numbers of isotopes
	 NucleusInfo(const std::string& n, int z, std::vector<int> a): name(n), Z(z), A(std::move(a)) {}
 };
 
 class BkgAnalyzer {
 public:
	 static constexpr int kNDetectors = 3;    // TOF, NaF, AGL
	 static constexpr int kNIsotopes = 3;     // Be7, Be9, Be10
	 static constexpr int kNHistTypes = 14;   // 扩展后的直方图类型数量
	 static const std::vector<NucleusInfo> kNuclei;
	 
	 BkgAnalyzer(const std::string& inputFile, const std::string& outputDir, const NucleusInfo& nucleus, int mass)
		 : inputFile_(inputFile), outputDir_(outputDir), nucleus_(nucleus), sourceFlux_(nullptr), currentMass_(mass) {
		 Initialize();
	 }
	 
	 ~BkgAnalyzer() { CleanupHistograms(); delete sourceFlux_; }
	 
	 // Main processing function
	 void Process() { 
		 LoadSourceFlux(); 
		 ProcessEvents(); 
		 CalculateAcceptance(); 
		 SaveResults(); 
	 }
 
 private:
	 std::string inputFile_, outputDir_;
	 NucleusInfo nucleus_; 
	 int currentMass_;
	 TF1* sourceFlux_;
	 std::vector<std::vector<std::vector<TH1F*>>> histograms_;    // [detector][isotope][type]
	 std::vector<TH1F*> acceptance_hists_;
	 
	 const std::array<std::string, kNDetectors> detectorNames_{"TOF", "NaF", "AGL"};
	 const std::array<std::string, kNHistTypes> histTypes_{
		 "NumEk", 
		 "genEk", "recEk", "genEk_reweighted", "recEk_reweighted",
		 "genEk_Cut2", "recEk_Cut2", "genEk_reweighted_Cut2", "recEk_reweighted_Cut2",
		 "genEk_Cut3", "recEk_Cut3", "genEk_reweighted_Cut3", "recEk_reweighted_Cut3",
		 "NumEk_reweighted"
	 };
	 
	 // Source flux fit points for spline interpolation
	 std::array<double, 13> kSourceFluxFitPoints{
		 2., 2.6, 3.2, 4, 6, 13.0, 30, 50.0, 100.0, 211.0, 330.0, 643.0, 2000.0
	 };
 
	 // Initialize histograms for all detectors and isotopes
	 void Initialize() {
		 histograms_.resize(kNDetectors);
		 for(auto& det : histograms_) { 
			 det.resize(kNIsotopes);
			 for(auto& iso : det) iso.resize(kNHistTypes, nullptr);
		 }
 
		 for(int det = 0; det < kNDetectors; ++det) {
			 for(int iso = 0; iso < kNIsotopes; ++iso) {
				 int mass = IsotopeData[3].mass_[iso];  // Z=4 (Be) is at index 3
				 const auto& ekBins = Binning::KineticEnergyBins[3][iso];
				 for(int type = 0; type < kNHistTypes; ++type) {
					 std::string name = Form("h_%s_%s_Be%d", histTypes_[type].c_str(), 
						 detectorNames_[det].c_str(), mass);
					 std::string xtitle = (histTypes_[type].find("rec") != std::string::npos) ? 
						 detectorNames_[det] + "Ek" : "geneEk";
					 
					 histograms_[det][iso][type] = new TH1F(name.c_str(), 
						 Form("%s;%s;Events", name.c_str(), xtitle.c_str()),
						 Constants::RIGIDITY_BINS, ekBins.data());
					 histograms_[det][iso][type]->Sumw2();
				 }
			 }
		 }
	 }
// Load source flux from data file
void LoadSourceFlux() {
    std::unique_ptr<TFile> fluxFile(TFile::Open(("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_" + 
        nucleus_.name + "Flux.root").c_str()));
    if(!fluxFile || !fluxFile->Get("graph1")) throw std::runtime_error("Cannot access flux file/graph");
    
    sourceFlux_ = SplineFit(static_cast<TGraph*>(fluxFile->Get("graph1")), 
        kSourceFluxFitPoints.data(), 
        kSourceFluxFitPoints.size(), 
        0x38, "b1e1", "f_SourceFlux", 1., 2000.0);
}

// Calculate weight for each event based on flux/MC ratio
double CalculateWeight(double rigidity) const {
    if (rigidity <= 0) return 0.0;
    
    // Calculate flux weight
    double flux_weight = sourceFlux_->Eval(rigidity) / sourceFlux_->Integral(1., 2000.);
    
    // Calculate MC weight
    double mc_weight = getMCFunction().Eval(rigidity) / getMCNorm();
    
    // Return flux/MC ratio
    return flux_weight / mc_weight;
}

// Calculate generated number for each bin
void CalculateBinGeneNum(double N_gen, int isoIndex, bool is_reweight) {
    std::unique_ptr<TF1> f_flux(is_reweight ? sourceFlux_ : 
        new TF1("f_flux", "TMath::Power(x,-1)", Acc::_xlow, Acc::_xup));
    double baseFluxIntegral = f_flux->Integral(Acc::_xlow, Acc::_xup);
    
    for(int det = 0; det < kNDetectors; ++det) {
        TH1F* hist = histograms_[det][isoIndex][is_reweight ? 13 : 0]; // 注意这里改成13，对应新的NumEk_reweighted位置
        const double* bins = hist->GetXaxis()->GetXbins()->GetArray();
        for(int bin = 0; bin < hist->GetNbinsX(); ++bin) {
            int mass = IsotopeData[3].mass_[isoIndex];  // Z=4 (Be) is at index 3
            double Rlow = betaToRigidity(kineticEnergyToBeta(bins[bin]), 4, mass, false);
            double Rup = betaToRigidity(kineticEnergyToBeta(bins[bin + 1]), 4, mass, false);
            hist->SetBinContent(bin + 1, N_gen * f_flux->Integral(Rlow, Rup) / baseFluxIntegral);
        }
    }
    if(is_reweight) f_flux.release();  // Don't delete sourceFlux_
}

// Process all events from the input file
void ProcessEvents() {
    std::unique_ptr<TFile> file(TFile::Open(inputFile_.c_str()));
    if(!file || !file->Get("hMCnum") || !file->Get("saveTree")) 
        throw std::runtime_error("Cannot access input file/contents");
    
    double N_gen = ((TH1F*)file->Get("hMCnum"))->GetBinContent(1);
    for(int iso = 0; iso < kNIsotopes; ++iso) {
        CalculateBinGeneNum(N_gen, iso, false);
        CalculateBinGeneNum(N_gen, iso, true);
    }

    struct EventData {
        double InnerRig, generatedRig, generatedEk;
        double ToFBeta, NaFBeta, AGLBeta, ToFEk, NaFEk, AGLEk;
        int DetectorAndBkgCut[3][3];
    } data;
    
    TTree* tree = (TTree*)file->Get("saveTree");
    tree->SetBranchAddress("InnerRig", &data.InnerRig);
    tree->SetBranchAddress("generatedRig", &data.generatedRig);
    tree->SetBranchAddress("generatedEk", &data.generatedEk);
    tree->SetBranchAddress("ToFEk", &data.ToFEk);
    tree->SetBranchAddress("NaFEk", &data.NaFEk);
    tree->SetBranchAddress("AGLEk", &data.AGLEk);
    tree->SetBranchAddress("DetectorAndBkgCut", data.DetectorAndBkgCut);

    // 新增的分支
    Int_t cutStatus;
    int mtrpar[9];
    float mevmom1[23], mparc[2];
    bool rich_NaF;
    tree->SetBranchAddress("cutStatus", &cutStatus);
    tree->SetBranchAddress("mtrpar", &mtrpar);
    tree->SetBranchAddress("mevmom1", &mevmom1);
    tree->SetBranchAddress("mparc", &mparc);
    tree->SetBranchAddress("rich_NaF", &rich_NaF);

    const std::array<double*, 3> recoEks{&data.ToFEk, &data.NaFEk, &data.AGLEk};
    for(Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);

        if(entry % 100000 == 1) cout<<entry<<endl;

        if(mevmom1[0]!=-1000) continue;

        double weight = CalculateWeight(data.generatedRig);
        
        for(int det = 0; det < kNDetectors; ++det) {
            for(int iso = 0; iso < kNIsotopes; ++iso) {
                int mass = IsotopeData[3].mass_[iso];
                
                // Cut1: 原有的DetectorAndBkgCut判断
                if(data.DetectorAndBkgCut[det][iso] == 1) {
                    histograms_[det][iso][1]->Fill(data.generatedEk);
                    histograms_[det][iso][2]->Fill(*recoEks[det]);
                    histograms_[det][iso][3]->Fill(data.generatedEk, weight);
                    histograms_[det][iso][4]->Fill(*recoEks[det], weight);
                }
                
                // 检查基本探测器cuts
                bool trackerCut = (cutStatus & 0x01) == 0x01;  // Bit 0
                bool tofCut = (cutStatus & 0x02) == 0x02;      // Bit 1
                bool richCut = (cutStatus & 0x04) == 0x04;     // Bit 2
                
                // 根据探测器类型确定cut条件
                bool detectorCut = false;
                if(det == 0) {  // TOF
                    detectorCut = trackerCut && tofCut;
                } else if(det == 1) {  // NaF
                    detectorCut = trackerCut && richCut && rich_NaF;
                } else {  // AGL
                    detectorCut = trackerCut && richCut && !rich_NaF;
                }

                // Cut2: detector cuts && mtrpar要求
                bool cut2Pass = detectorCut && (
                    (mass == 7 && mtrpar[0] == 63) ||
                    (mass == 9 && mtrpar[0] == 64) ||
                    (mass == 10 && mtrpar[0] == 114)
                );
                if(cut2Pass) {
                    histograms_[det][iso][5]->Fill(data.generatedEk);
                    histograms_[det][iso][6]->Fill(*recoEks[det]);
                    histograms_[det][iso][7]->Fill(data.generatedEk, weight);
                    histograms_[det][iso][8]->Fill(*recoEks[det], weight);
                }

                // Cut3: detector cuts && mevmom1要求 && mparc要求(只对Be7)
                bool cut3Pass = detectorCut && mevmom1[0] == -1000;
                if(mass == 7) cut3Pass = cut3Pass && (mparc[0] == 4);
                if(cut3Pass) {
                    histograms_[det][iso][9]->Fill(data.generatedEk);
                    histograms_[det][iso][10]->Fill(*recoEks[det]);
                    histograms_[det][iso][11]->Fill(data.generatedEk, weight);
                    histograms_[det][iso][12]->Fill(*recoEks[det], weight);
                }
                
            }
        }
    }
}
// Calculate acceptance histograms
void CalculateAcceptance() {
    acceptance_hists_.reserve(kNDetectors * kNIsotopes * 4);
    const double scale_factor = 3.9 * 3.9 * TMath::Pi();
    
    for(int det = 0; det < kNDetectors; ++det) {
        for(int iso = 0; iso < kNIsotopes; ++iso) {
            std::array<std::pair<int,int>, 4> hist_pairs{{{1,0}, {2,0}, {3,13}, {4,13}}};
            
            for(const auto& [num, den] : hist_pairs) {
                TH1F* h_num = (TH1F*)histograms_[det][iso][num]->Clone();
                std::string acc_type = (num == 1 || num == 3) ? "gen" : "rec";
                std::string weight_tag = (num > 2) ? "reweighted" : "";
                
                int mass = IsotopeData[3].mass_[iso];  // Z=4 (Be) is at index 3
                std::string xtitle = histograms_[det][iso][num]->GetXaxis()->GetTitle();
                
                std::string name = Form("h_acc_%s_%s_%s_Be%d", acc_type.c_str(), weight_tag.c_str(),
                    detectorNames_[det].c_str(), mass);
                
                h_num->SetName(name.c_str());
                h_num->SetTitle(Form("%s %s %s Be%d Acceptance;%s;Acceptance (m^{2} sr)", 
                    weight_tag.c_str(), acc_type.c_str(), detectorNames_[det].c_str(), mass, xtitle.c_str()));
                
                h_num->Divide(histograms_[det][iso][den]);
                h_num->Scale(scale_factor);
                acceptance_hists_.push_back(h_num);
            }
        }
    }
}

// Save results to output file
void SaveResults() {
    std::unique_ptr<TFile> outFile(TFile::Open((outputDir_ + "/" + nucleus_.name + 
        std::to_string(currentMass_) + "_bkg_analysis.root").c_str(), "RECREATE"));
    if(!outFile) throw std::runtime_error("Cannot create output file");
    cout << outFile->GetName() << endl;

    for(const auto& det : histograms_)
        for(const auto& iso : det)
            for(const auto& hist : iso)
                if(hist) hist->Write();
                
    for(const auto& hist : acceptance_hists_) 
        if(hist) hist->Write();
        
    if(sourceFlux_) sourceFlux_->Write();
}

// Cleanup allocated memory
void CleanupHistograms() {
    for(auto& det : histograms_)
        for(auto& iso : det)
            for(auto& hist : iso) { 
                delete hist; 
                hist = nullptr; 
            }
            
    for(auto& hist : acceptance_hists_) { 
        delete hist; 
        hist = nullptr; 
    }
}
};

// Initialize supported nuclei configurations
const std::vector<NucleusInfo> BkgAnalyzer::kNuclei = {
    {"Beryllium", 4, {9, 10}},  
    {"Boron", 5, {10, 11}}, 
    {"Carbon", 6, {12}},
    {"Nitrogen", 7, {14, 15}}, 
    {"Oxygen", 8, {16}}
};

// Main function to run the analysis
void cal_bkghist() {
    try {
        for(const auto& nucleus : BkgAnalyzer::kNuclei) {
            for(const auto& mass : nucleus.A) {
                BkgAnalyzer(("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MC" + nucleus.name + 
                    std::to_string(mass) + ".root").c_str(),
                    "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg", nucleus, mass).Process();
                std::cout << "Completed analysis for " << nucleus.name << mass << std::endl;
            }
        }
    } catch(const std::exception& e) { 
        std::cerr << "Error: " << e.what() << std::endl; 
    }
}