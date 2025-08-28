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
 
 struct NucleusInfo {
	 std::string name;
	 int Z;
	 std::vector<int> A;
	 NucleusInfo(const std::string& n, int z, std::vector<int> a): name(n), Z(z), A(std::move(a)) {}
 };
 
 class BkgAnalyzer {
	 public:
		 static constexpr int kNDetectors = 3;    // TOF, NaF, AGL
		 static constexpr int kNIsotopes = 3;     // Be7, Be9, Be10
		 static const std::vector<NucleusInfo> kNuclei;
 
		 BkgAnalyzer(const std::string& inputFile, const std::string& outputDir, const NucleusInfo& nucleus, int mass)
			 : inputFile_(inputFile), outputDir_(outputDir), nucleus_(nucleus), sourceFlux_(nullptr), currentMass_(mass) {
				 Initialize();
			 }
 
		 ~BkgAnalyzer() { CleanupHistograms(); delete sourceFlux_; }
 
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
 
		 // 新的直方图存储结构
		 std::vector<std::vector<TH1F*>> h_num_;  // [iso][2] 每个同位素的reweight前后的总数
		 std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>> h_bkg_study_;  // [iso][2][2][3][3] iso/reweight/sel/bkg/det
		 std::vector<TH1F*> acceptance_hists_;   
 
		 const std::array<std::string, kNDetectors> detectorNames_{"TOF", "NaF", "AGL"};
 
		 std::array<double, 13> kSourceFluxFitPoints{
			 2., 2.6, 3.2, 4, 6, 13.0, 30, 50.0, 100.0, 211.0, 330.0, 643.0, 2000.0
		 };
 
		 void Initialize() {
			 const std::array<std::string, 2> selTypes{"Complete", "Loose"};
			 const std::array<std::string, 3> bkgTypes{"Bkg", "NoBkg", "StrictBkg"};
 
			 // 初始化h_num_，每个同位素两个(原始和reweight)
			 h_num_.resize(kNIsotopes);
			 for(int iso = 0; iso < kNIsotopes; ++iso) {
				 h_num_[iso].resize(2, nullptr);
				 int mass = IsotopeData[3].mass_[iso];  // Z=4 (Be) is at index 3
				 const auto& ekBins = Binning::KineticEnergyBins[3][iso];
 
				 for(int rw = 0; rw < 2; ++rw) {
					 std::string rwTag = rw ? "reweight" : "orig";
					 std::string name = Form("h_num_%s_Be%d", rwTag.c_str(), mass);
					 std::string title = Form("%s Be%d Number;Ek;Events", rwTag.c_str(), mass);
 
					 h_num_[iso][rw] = new TH1F(name.c_str(), title.c_str(),
							 Constants::RIGIDITY_BINS, ekBins.data());
					 h_num_[iso][rw]->Sumw2();
				 }
			 }
 
			 // 初始化h_bkg_study_
			 h_bkg_study_.resize(kNIsotopes);
			 for(int iso = 0; iso < kNIsotopes; ++iso) {
				 int mass = IsotopeData[3].mass_[iso];
				 const auto& ekBins = Binning::KineticEnergyBins[3][iso];
 
				 h_bkg_study_[iso].resize(2);  // orig/reweight
				 for(auto& v1 : h_bkg_study_[iso]) {
					 v1.resize(2);  // Complete/Loose
					 for(auto& v2 : v1) {
						 v2.resize(3);  // Bkg/NoBkg/StrictBkg
						 for(auto& v3 : v2) {
							 v3.resize(3, nullptr);  // TOF/NaF/AGL
						 }
					 }
				 }
 
				 for(int rw = 0; rw < 2; ++rw) {
					 std::string rwTag = rw ? "reweight" : "orig";
					 for(int sel = 0; sel < 2; ++sel) {
						 for(int bkg = 0; bkg < 3; ++bkg) {
							 for(int det = 0; det < 3; ++det) {
								 std::string name = Form("h_bkg_%s_%s_%s_%s_Be%d",
										 rwTag.c_str(), selTypes[sel].c_str(), 
										 bkgTypes[bkg].c_str(), detectorNames_[det].c_str(), mass);
 
								 std::string title = Form("%s %s %s %s Be%d;%sEk;Events",
										 rwTag.c_str(), selTypes[sel].c_str(), bkgTypes[bkg].c_str(),
										 detectorNames_[det].c_str(), mass,
										 sel == 0 ? detectorNames_[det].c_str() : "gene");  // Complete用detector Ek，Loose用gene Ek
 
								 h_bkg_study_[iso][rw][sel][bkg][det] = new TH1F(name.c_str(),
										 title.c_str(), Constants::RIGIDITY_BINS, ekBins.data());
								 h_bkg_study_[iso][rw][sel][bkg][det]->Sumw2();
							 }
						 }
					 }
				 }
			 }
		 }
 
		 void LoadSourceFlux() {
			 std::unique_ptr<TFile> fluxFile(TFile::Open(("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_" + 
							 nucleus_.name + "Flux.root").c_str()));
			 if(!fluxFile || !fluxFile->Get("graph1")) throw std::runtime_error("Cannot access flux file/graph");
 
			 sourceFlux_ = SplineFit(static_cast<TGraph*>(fluxFile->Get("graph1")), 
					 kSourceFluxFitPoints.data(), 
					 kSourceFluxFitPoints.size(), 
					 0x38, "b1e1", "f_SourceFlux", 1., 2000.0);
		 }
 
		 double CalculateWeight(double rigidity) const {
			 if (rigidity <= 0) return 0.0;
			 double flux_weight = sourceFlux_->Eval(rigidity) / sourceFlux_->Integral(1., 2000.);
			 double mc_weight = getMCFunction().Eval(rigidity) / getMCNorm();
			 return flux_weight / mc_weight;
		 }
 
		 void CalculateBinGeneNum(double N_gen, int isoIndex, bool is_reweight) {
			 std::unique_ptr<TF1> f_flux(is_reweight ? sourceFlux_ : 
					 new TF1("f_flux", "TMath::Power(x,-1)", Acc::_xlow, Acc::_xup));
			 double baseFluxIntegral = f_flux->Integral(Acc::_xlow, Acc::_xup);
 
			 TH1F* hist = h_num_[isoIndex][is_reweight ? 1 : 0];
			 const double* bins = hist->GetXaxis()->GetXbins()->GetArray();
			 for(int bin = 0; bin < hist->GetNbinsX(); ++bin) {
				 int mass = IsotopeData[3].mass_[isoIndex];
				 double Rlow = betaToRigidity(kineticEnergyToBeta(bins[bin]), 4, mass, false);
				 double Rup = betaToRigidity(kineticEnergyToBeta(bins[bin + 1]), 4, mass, false);
				 hist->SetBinContent(bin + 1, N_gen * f_flux->Integral(Rlow, Rup) / baseFluxIntegral);
			 }
			 if(is_reweight) f_flux.release();  // Don't delete sourceFlux_
		 }
 
		 void ProcessEvents() {
			 std::unique_ptr<TFile> file(TFile::Open(inputFile_.c_str()));
			 if(!file || !file->Get("hMCnum") || !file->Get("saveTree")) 
				 throw std::runtime_error("Cannot access input file/contents");
 
			 double N_gen = ((TH1F*)file->Get("hMCnum"))->GetBinContent(1);
 
			 TTree* tree = (TTree*)file->Get("saveTree");
	
			 // 为每个同位素计算bin中的生成数
			 for(int iso = 0; iso < kNIsotopes; ++iso) {
				 CalculateBinGeneNum(N_gen, iso, false);  // orig
				 CalculateBinGeneNum(N_gen, iso, true);   // reweight
			 }
			 double InnerRig, generatedRig, generatedEk;
			 double TOFEk, NaFEk, AGLEk;
			 bool rich_NaF;
			 int BkgStudyCut[2][3][3];  // [sel][bkg][det]
			 int MCInterCut[3];  // [iso]
			 int mtrpar[9];
 
			 tree->SetBranchAddress("InnerRig", &InnerRig);
			 tree->SetBranchAddress("generatedRig", &generatedRig);
			 tree->SetBranchAddress("generatedEk", &generatedEk);
			 tree->SetBranchAddress("TOFEk", &TOFEk);
			 tree->SetBranchAddress("NaFEk", &NaFEk);
			 tree->SetBranchAddress("AGLEk", &AGLEk);
			 tree->SetBranchAddress("rich_NaF", &rich_NaF);
			 tree->SetBranchAddress("BkgStudyCut", BkgStudyCut);
			 tree->SetBranchAddress("MCInterCut", MCInterCut);
			 tree->SetBranchAddress("mtrpar", mtrpar);
 
			 const std::array<double*, 3> recoEks{&TOFEk, &NaFEk, &AGLEk};
			 int isoCode[3] = {63, 64, 114};
 
			 for(Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
				 tree->GetEntry(entry);
 
				 if(entry % 100000 == 0) cout << "Processing entry " << entry << endl;
				 //if(MCInterCut[0] == 0 && MCInterCut[1] == 0 && MCInterCut[2] == 0) continue;
				 if(mtrpar[7]!=63 && mtrpar[7]!=64 && mtrpar[7]!=114) continue;
 
				 double weight = CalculateWeight(generatedRig);
 
				 // 对每个同位素处理
				 for(int iso = 0; iso < kNIsotopes; ++iso) {
					 int mass = IsotopeData[3].mass_[iso];  // Z=4 (Be) is at index 3
 
					 // 对每个探测器处理
					 for(int det = 0; det < kNDetectors; ++det) {
						 // 对每个选择类型处理 (Complete/Loose)
						 for(int sel = 0; sel < 2; ++sel) {
							 // 对每个背景类型处理 (Bkg/NoBkg/StrictBkg)
							 for(int bkg = 0; bkg < 3; ++bkg) {
								 if(BkgStudyCut[sel][bkg][det] == 1) {
									 // 填充原始直方图
									 //if(MCInterCut[iso]) {  // mc truth info
									 if(mtrpar[7] == isoCode[iso]) {  // mc truth info
										 if(generatedEk > 0) {
											 h_bkg_study_[iso][0][sel][bkg][det]->Fill(generatedEk);
											 h_bkg_study_[iso][1][sel][bkg][det]->Fill(generatedEk, weight);
										 }
									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }
 
		 void CalculateAcceptance() {
			 const double scale_factor = 3.9 * 3.9 * TMath::Pi();
			 // 对每个同位素计算接收度
			 for(int iso = 0; iso < kNIsotopes; ++iso) {
				 int mass = IsotopeData[3].mass_[iso];
				 // 对orig和reweight分别计算
				 for(int rw = 0; rw < 2; rw++) {
					 std::string rwTag = rw ? "reweight" : "orig";
					 // 对Complete和Loose分别计算
					 for(int sel = 0; sel < 2; sel++) {
						 std::string selTag = sel == 0 ? "Complete" : "Loose";
						 // 对不同背景cut计算
						 for(int bkg = 0; bkg < 3; bkg++) {
							 std::string bkgTag = bkg == 0 ? "Bkg" : (bkg == 1 ? "NoBkg" : "StrictBkg");
							 // 对每个探测器计算
							 for(int det = 0; det < 3; det++) {
								 // 创建接收度直方图
								 std::string name = Form("h_acc_%s_%s_%s_%s_Be%d", 
										 rwTag.c_str(), selTag.c_str(), bkgTag.c_str(), 
										 detectorNames_[det].c_str(), mass);
								 std::string title = Form("%s %s %s %s Be%d Acceptance;%sEk;Acceptance (m^{2} sr)",
										 rwTag.c_str(), selTag.c_str(), bkgTag.c_str(),
										 detectorNames_[det].c_str(), mass,
										 sel == 0 ? detectorNames_[det].c_str() : "gene");
								 TH1F* h_acc = (TH1F*)h_bkg_study_[iso][rw][sel][bkg][det]->Clone(name.c_str());
								 h_acc->SetTitle(title.c_str());
								 // 除以对应的num直方图并乘以比例因子
								 h_acc->Divide(h_num_[iso][rw]);
								 h_acc->Scale(scale_factor);
								 acceptance_hists_.push_back(h_acc);
							 }
						 }
					 }
				 }
			 }
		 }
 
		 void SaveResults() {
			 std::unique_ptr<TFile> outFile(TFile::Open((outputDir_ + "/" + nucleus_.name + 
							 std::to_string(currentMass_) + "_bkg_analysis_mtrpar7.root").c_str(), "RECREATE"));
			 if(!outFile) throw std::runtime_error("Cannot create output file");
			 std::cout << "Saving results to: " << outFile->GetName() << std::endl;
 
			 // Save h_num_
			 for(const auto& iso_hists : h_num_) {
				 for(const auto& hist : iso_hists) {
					 if(hist) hist->Write();
				 }
			 }
 
			 // Save h_bkg_study_
			 for(const auto& iso_hists : h_bkg_study_) {        // [iso]
				 for(const auto& rw_hists : iso_hists) {        // [reweight]
					 for(const auto& sel_hists : rw_hists) {    // [sel]
						 for(const auto& bkg_hists : sel_hists) { // [bkg]
							 for(const auto& hist : bkg_hists) {  // [det]
								 if(hist) hist->Write();
							 }
						 }
					 }
				 }
			 }
 
			 // Save acceptance histograms
			 for(const auto& hist : acceptance_hists_) {
				 if(hist) hist->Write();
			 }
 
			 if(sourceFlux_) sourceFlux_->Write();
		 }
 
		 void CleanupHistograms() {
			 // Cleanup h_num_
			 for(auto& iso_hists : h_num_) {
				 for(auto& hist : iso_hists) {
					 delete hist;
					 hist = nullptr;
				 }
			 }
 
			 // Cleanup h_bkg_study_
			 for(auto& iso_hists : h_bkg_study_) {        // [iso]
				 for(auto& rw_hists : iso_hists) {        // [reweight]
					 for(auto& sel_hists : rw_hists) {    // [sel]
						 for(auto& bkg_hists : sel_hists) { // [bkg]
							 for(auto& hist : bkg_hists) {  // [det]
								 delete hist;
								 hist = nullptr;
							 }
						 }
					 }
				 }
			 }
 
			 // Cleanup acceptance histograms
			 for(auto& hist : acceptance_hists_) {
				 delete hist;
				 hist = nullptr;
			 }
		 }
 };
 
 const std::vector<NucleusInfo> BkgAnalyzer::kNuclei = {
	 {"Beryllium", 4, {7, 9 ,10}}, 
	 //{"Boron", 5, {10, 11}}, 
	 //{"Carbon", 6, {12}},
	 //{"Nitrogen", 7, {14, 15}}, 
	 //{"Oxygen", 8, {16}}
 };
 
 // Main function to run the analysis
 void cal_bkghist_par8() {
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
 