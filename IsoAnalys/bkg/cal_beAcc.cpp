/*********************************************************** 
 * Analysis of background events and acceptance from MC fragmentation.
 * Special processor for Beryllium isotopes with mixed acceptance calculation.
 * Author: Yan Zixuan  Date: 20250327
 ***********************************************************/
 #include "../Tool.h"
 #include <TFile.h> 
 #include <TTree.h>
 #include <TGraphAsymmErrors.h>
 #include <memory>
 #include <array>
 #include <vector>
 #include <map>
 #include <stdexcept>
 
 using namespace AMS_Iso;
 
 class BeBkgAnalyzer {
 public:
 static constexpr int kNDetectors = 3;    // TOF, NaF, AGL
    
 // 现在这些是非static成员变量
 const std::array<std::string, kNDetectors> detectorNames_;
 const std::array<std::string, 2> selTypes;
 const std::array<std::string, 3> bkgTypes;
 const std::map<int, double> isotopeMixWeights;

 // 在构造函数的初始化列表中初始化这些const成员
 BeBkgAnalyzer(const std::string& be7File, const std::string& be9File, 
			   const std::string& be10File, const std::string& outputDir)
	 : detectorNames_{{"TOF", "NaF", "AGL"}}
	 , selTypes{{"Complete", "Loose"}}
	 , bkgTypes{{"Bkg", "NoBkg", "StrictBkg"}}
	 , isotopeMixWeights{{7, 0.1}, {9, 0.3}, {10, 0.6}}
	 , outputDir_(outputDir)
 {
	 inputFiles_[7] = be7File;
	 inputFiles_[9] = be9File;
	 inputFiles_[10] = be10File;
	 Initialize();
 }
 
	 ~BeBkgAnalyzer() { 
		 CleanupHistograms(); 
		 for(auto& flux : sourceFluxes_) {
			 delete flux.second;
		 }
	 }
 
	 void Process() {
		 LoadSourceFluxes();
		 ProcessFiles();
		 CalculateAcceptances();
		 MixHistograms();
		 CalculateMixedAcceptance();
		 FitAcceptances();
		 SaveResults();
	 }
 
 private:
 template<typename T, typename Func>
 void ForEachHistogram(T& container, Func f) {
	 if constexpr (std::is_same_v<T, std::vector<TH1F*>>) {
		 // 处理一维数组 (h_num_mix_)
		 for (auto& hist : container) {
			 if (hist) f(hist);
		 }
	 }
	 else if constexpr (std::is_same_v<T, std::map<int, std::vector<TH1F*>>>) {
		 // 处理h_num_
		 for (auto& [_, histVec] : container) {
			 for (auto& hist : histVec) {
				 if (hist) f(hist);
			 }
		 }
	 }
	 else {
		 // 处理多维数组 (h_bkg_study_, h_acc_, h_bkg_study_mix_, h_acc_mix_)
		 for (auto& v1 : container) {
			 for (auto& v2 : v1) {
				 for (auto& v3 : v2) {
					 for (auto& hist : v3) {
						 if (hist) f(hist);
					 }
				 }
			 }
		 }
	 }
 }

	 std::map<int, std::string> inputFiles_;  // mass -> file path
	 std::string outputDir_;
	 
	 // Histograms for each isotope
	 std::map<int, std::vector<TH1F*>> h_num_;  // mass -> [rw]
	 std::map<int, std::vector<std::vector<std::vector<std::vector<TH1F*>>>>> h_bkg_study_;  // mass -> [rw][sel][bkg][det]
	 std::map<int, std::vector<std::vector<std::vector<std::vector<TH1F*>>>>> h_acc_;  // mass -> [rw][sel][bkg][det]
	 
	 // Mixed histograms
	 std::vector<TH1F*> h_num_mix_;  // [rw]
	 std::vector<std::vector<std::vector<std::vector<TH1F*>>>> h_bkg_study_mix_;  // [rw][sel][bkg][det]
	 std::vector<std::vector<std::vector<std::vector<TH1F*>>>> h_acc_mix_;  // [rw][sel][bkg][det]
	 
	 // Source flux functions
	 std::map<int, TF1*> sourceFluxes_;  // mass -> flux function
 
	 std::array<double, 13> kSourceFluxFitPoints{
		 2., 2.6, 3.2, 4, 6, 13.0, 30, 50.0, 100.0, 211.0, 330.0, 643.0, 2000.0
	 };
 
	 void Initialize() {
		 // Initialize histograms for each isotope
		 for(const auto& [mass, _] : inputFiles_) {
			 InitializeHistogramsForIsotope(mass);
		 }
		 
		 // Initialize mixed histograms
		 InitializeMixedHistograms();
	 }
 
	 void InitializeHistogramsForIsotope(int mass) {
        int massIndex = mass == 7 ? 0 : (mass == 9 ? 1 : 2);
		 const auto& ekBins = Binning::KineticEnergyBins[3][massIndex];  // Be7 starts at index 0
 
		 // Initialize h_num_
		 h_num_[mass].resize(2, nullptr);
		 for(int rw = 0; rw < 2; ++rw) {
			 std::string rwTag = rw ? "reweight" : "orig";
			 std::string name = Form("h_num_%s_Be%d", rwTag.c_str(), mass);
			 std::string title = Form("%s Be%d Number;Ek;Events", rwTag.c_str(), mass);
 
			 h_num_[mass][rw] = new TH1F(name.c_str(), title.c_str(),
					 Constants::RIGIDITY_BINS, ekBins.data());
			 h_num_[mass][rw]->Sumw2();
		 }
 
		 // Initialize h_bkg_study_
		 h_bkg_study_[mass].resize(2);  // orig/reweight
		 for(auto& v1 : h_bkg_study_[mass]) {
			 v1.resize(2);  // Complete/Loose
			 for(auto& v2 : v1) {
				 v2.resize(3);  // Bkg/NoBkg/StrictBkg
				 for(auto& v3 : v2) {
					 v3.resize(3, nullptr);  // TOF/NaF/AGL
				 }
			 }
		 }
 
		 // Initialize all background study histograms
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
								 sel == 0 ? detectorNames_[det].c_str() : "gene");
 
						 h_bkg_study_[mass][rw][sel][bkg][det] = new TH1F(name.c_str(),
								 title.c_str(), Constants::RIGIDITY_BINS, ekBins.data());
						 h_bkg_study_[mass][rw][sel][bkg][det]->Sumw2();
					 }
				 }
			 }
		 }
	 }
	 void InitializeMixedHistograms() {
        const auto& ekBins = Binning::KineticEnergyBins[3][0];  // Use Be7 binning for mixed histograms

        // Initialize mixed h_num_
        h_num_mix_.resize(2, nullptr);
        for(int rw = 0; rw < 2; ++rw) {
            std::string rwTag = rw ? "reweight" : "orig";
            std::string name = Form("h_num_%s_BeTotal", rwTag.c_str());
            std::string title = Form("%s BeTotal Number;Ek;Events", rwTag.c_str());

            h_num_mix_[rw] = new TH1F(name.c_str(), title.c_str(),
                    Constants::RIGIDITY_BINS, ekBins.data());
            h_num_mix_[rw]->Sumw2();
        }

        // Initialize mixed h_bkg_study_
        h_bkg_study_mix_.resize(2);  // orig/reweight
        for(auto& v1 : h_bkg_study_mix_) {
            v1.resize(2);  // Complete/Loose
            for(auto& v2 : v1) {
                v2.resize(3);  // Bkg/NoBkg/StrictBkg
                for(auto& v3 : v2) {
                    v3.resize(3, nullptr);  // TOF/NaF/AGL
                }
            }
        }

        // Initialize all mixed background study histograms
        for(int rw = 0; rw < 2; ++rw) {
            std::string rwTag = rw ? "reweight" : "orig";
            for(int sel = 0; sel < 2; ++sel) {
                for(int bkg = 0; bkg < 3; ++bkg) {
                    for(int det = 0; det < 3; ++det) {
                        std::string name = Form("h_bkg_%s_%s_%s_%s_BeTotal",
                                rwTag.c_str(), selTypes[sel].c_str(), 
                                bkgTypes[bkg].c_str(), detectorNames_[det].c_str());

                        std::string title = Form("%s %s %s %s BeTotal;%sEk;Events",
                                rwTag.c_str(), selTypes[sel].c_str(), bkgTypes[bkg].c_str(),
                                detectorNames_[det].c_str(),
                                sel == 0 ? detectorNames_[det].c_str() : "gene");

                        h_bkg_study_mix_[rw][sel][bkg][det] = new TH1F(name.c_str(),
                                title.c_str(), Constants::RIGIDITY_BINS, ekBins.data());
                        h_bkg_study_mix_[rw][sel][bkg][det]->Sumw2();
                    }
                }
            }
        }
    }

    void LoadSourceFluxes() {
        for(const auto& [mass, _] : inputFiles_) {
            std::unique_ptr<TFile> fluxFile(TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BerylliumFlux.root"));
            if(!fluxFile || !fluxFile->Get("graph1")) {
                throw std::runtime_error("Cannot access flux file/graph for Be" + std::to_string(mass));
            }

            sourceFluxes_[mass] = SplineFit(static_cast<TGraph*>(fluxFile->Get("graph1")), 
                    kSourceFluxFitPoints.data(), 
                    kSourceFluxFitPoints.size(), 
                    0x38, "b1e1", ("f_SourceFlux_Be" + std::to_string(mass)).c_str(), 1., 2000.0);
        }
    }

    double CalculateWeight(double rigidity, int mass) const {
        if (rigidity <= 0) return 0.0;
        auto fluxIt = sourceFluxes_.find(mass);
        if(fluxIt == sourceFluxes_.end()) return 0.0;
        
        double flux_weight = fluxIt->second->Eval(rigidity) / fluxIt->second->Integral(1., 2000.);
        double mc_weight = getMCFunction().Eval(rigidity) / getMCNorm();
        return flux_weight / mc_weight;
    }

    void CalculateBinGeneNum(double N_gen, int mass, bool is_reweight) {
        auto fluxIt = sourceFluxes_.find(mass);
        if(fluxIt == sourceFluxes_.end()) return;

        std::unique_ptr<TF1> f_flux(is_reweight ? fluxIt->second : 
                new TF1("f_flux", "TMath::Power(x,-1)", Acc::_xlow, Acc::_xup));
        double baseFluxIntegral = f_flux->Integral(Acc::_xlow, Acc::_xup);

        TH1F* hist = h_num_[mass][is_reweight ? 1 : 0];
        const double* bins = hist->GetXaxis()->GetXbins()->GetArray();
        for(int bin = 0; bin < hist->GetNbinsX(); ++bin) {
            double Rlow = betaToRigidity(kineticEnergyToBeta(bins[bin]), 4, mass, false);
            double Rup = betaToRigidity(kineticEnergyToBeta(bins[bin + 1]), 4, mass, false);
            hist->SetBinContent(bin + 1, N_gen * f_flux->Integral(Rlow, Rup) / baseFluxIntegral);
        }
        if(is_reweight) f_flux.release();  // Don't delete sourceFlux_
    }
	void ProcessFiles() {
        for(const auto& [mass, inputFile] : inputFiles_) {
            std::unique_ptr<TFile> file(TFile::Open(inputFile.c_str()));
            if(!file || !file->Get("hMCnum") || !file->Get("saveTree")) {
                throw std::runtime_error("Cannot access input file/contents for Be" + std::to_string(mass));
            }

            double N_gen = ((TH1F*)file->Get("hMCnum"))->GetBinContent(1);
            
            // Calculate bin generation numbers
            CalculateBinGeneNum(N_gen, mass, false);  // orig
            CalculateBinGeneNum(N_gen, mass, true);   // reweight

            ProcessEvents(file.get(), mass);
        }
    }

    void ProcessEvents(TFile* file, int mass) {
        TTree* tree = (TTree*)file->Get("saveTree");
        
        double InnerRig, generatedRig, generatedEk;
        double TOFEk, NaFEk, AGLEk;
        bool rich_NaF;
        int BkgStudyCut[2][3][3];  // [sel][bkg][det]

        tree->SetBranchAddress("InnerRig", &InnerRig);
        tree->SetBranchAddress("generatedRig", &generatedRig);
        tree->SetBranchAddress("generatedEk", &generatedEk);
        tree->SetBranchAddress("TOFEk", &TOFEk);
        tree->SetBranchAddress("NaFEk", &NaFEk);
        tree->SetBranchAddress("AGLEk", &AGLEk);
        tree->SetBranchAddress("rich_NaF", &rich_NaF);
        tree->SetBranchAddress("BkgStudyCut", BkgStudyCut);

        const std::array<double*, 3> recoEks{&TOFEk, &NaFEk, &AGLEk};

        for(Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);

            if(entry % 400000 == 0) cout << "Processing Be" << mass << " entry " << entry << endl;

            double weight = CalculateWeight(generatedRig, mass);

            // 对每个探测器处理
            for(int det = 0; det < kNDetectors; ++det) {
                // 对每个选择类型处理 (Complete/Loose)
                for(int sel = 0; sel < 2; ++sel) {
                    // 对每个背景类型处理 (Bkg/NoBkg/StrictBkg)
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        if(BkgStudyCut[sel][bkg][det] == 1) {
                            if(generatedEk > 0) {
                                h_bkg_study_[mass][0][sel][bkg][det]->Fill(generatedEk);
                                h_bkg_study_[mass][1][sel][bkg][det]->Fill(generatedEk, weight);
                            }
                        }
                    }
                }
            }
        }
    }

    void CalculateAcceptances() {
        const double scale_factor = 3.9 * 3.9 * TMath::Pi();
        
        // Initialize acceptance histograms for each isotope
        for(const auto& [mass, _] : inputFiles_) {
            h_acc_[mass].resize(2);  // orig/reweight
            for(auto& v1 : h_acc_[mass]) {
                v1.resize(2);  // Complete/Loose
                for(auto& v2 : v1) {
                    v2.resize(3);  // Bkg/NoBkg/StrictBkg
                    for(auto& v3 : v2) {
                        v3.resize(3, nullptr);  // TOF/NaF/AGL
                    }
                }
            }
        }

        // Calculate acceptance for each isotope
        for(const auto& [mass, _] : inputFiles_) {
            for(int rw = 0; rw < 2; ++rw) {
                for(int sel = 0; sel < 2; ++sel) {
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        for(int det = 0; det < 3; ++det) {
                            std::string name = Form("h_acc_%s_%s_%s_%s_Be%d",
                                    (rw ? "reweight" : "orig"), selTypes[sel].c_str(),
                                    bkgTypes[bkg].c_str(), detectorNames_[det].c_str(), mass);
                            
                            h_acc_[mass][rw][sel][bkg][det] = 
                                (TH1F*)h_bkg_study_[mass][rw][sel][bkg][det]->Clone(name.c_str());
                            h_acc_[mass][rw][sel][bkg][det]->Divide(h_num_[mass][rw]);
                            h_acc_[mass][rw][sel][bkg][det]->Scale(scale_factor);
                        }
                    }
                }
            }
        }
    }

	void MixHistograms() {
		for(int rw = 0; rw < 2; ++rw) {
			// 清零目标直方图
			h_num_mix_[rw]->Reset();
			
			// 加权相加
			for(const auto& [mass, weight] : isotopeMixWeights) {
				h_num_mix_[rw]->Add(h_num_[mass][rw], weight);
			}
		}
		
		// 同样处理 h_bkg_study
		for(int rw = 0; rw < 2; ++rw) {
			for(int sel = 0; sel < 2; ++sel) {
				for(int bkg = 0; bkg < 3; ++bkg) {
					for(int det = 0; det < 3; ++det) {
						auto hist = h_bkg_study_mix_[rw][sel][bkg][det];
						hist->Reset();
						for(const auto& [mass, weight] : isotopeMixWeights) {
							hist->Add(h_bkg_study_[mass][rw][sel][bkg][det], weight);
						}
					}
				}
			}
		}
	}

	void CalculateMixedAcceptance() {
        const double scale_factor = 3.9 * 3.9 * TMath::Pi();
        
        // Initialize mixed acceptance histograms
        h_acc_mix_.resize(2);  // orig/reweight
        for(auto& v1 : h_acc_mix_) {
            v1.resize(2);  // Complete/Loose
            for(auto& v2 : v1) {
                v2.resize(3);  // Bkg/NoBkg/StrictBkg
                for(auto& v3 : v2) {
                    v3.resize(3, nullptr);  // TOF/NaF/AGL
                }
            }
        }

        // Calculate mixed acceptance
        for(int rw = 0; rw < 2; ++rw) {
            for(int sel = 0; sel < 2; ++sel) {
                for(int bkg = 0; bkg < 3; ++bkg) {
                    for(int det = 0; det < 3; ++det) {
                        std::string name = Form("h_acc_%s_%s_%s_%s_BeTotal",
                                (rw ? "reweight" : "orig"), selTypes[sel].c_str(),
                                bkgTypes[bkg].c_str(), detectorNames_[det].c_str());
                        
                        h_acc_mix_[rw][sel][bkg][det] = 
                            (TH1F*)h_bkg_study_mix_[rw][sel][bkg][det]->Clone(name.c_str());
                        h_acc_mix_[rw][sel][bkg][det]->Divide(h_num_mix_[rw]);
                        h_acc_mix_[rw][sel][bkg][det]->Scale(scale_factor);
                    }
                }
            }
        }
    }

    void FitAcceptances() {
        // Fit individual isotope acceptances
        for(const auto& [mass, accHists] : h_acc_) {
            for(int rw = 0; rw < 2; ++rw) {
                for(int sel = 0; sel < 2; ++sel) {
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        for(int det = 0; det < 3; ++det) {
                            auto& xpoints = (mass == 7) ? Acc::xpointsEkfor7[det] :
                                          (mass == 9) ? Acc::xpointsEkfor9[det] :
                                                       Acc::xpointsEkfor10[det];
                            
                            std::string fitName = Form("f_acc_%s_%s_%s_%s_Be%d",
                                    (rw ? "reweight" : "orig"), selTypes[sel].c_str(),
                                    bkgTypes[bkg].c_str(), detectorNames_[det].c_str(), mass);
                            
                            auto f_acc = SplineFit(h_acc_[mass][rw][sel][bkg][det], 
                                                 xpoints.data(), xpoints.size(),
                                                 0x38, "b2e2", fitName.c_str(), 
                                                 xpoints[0], 650);
                        }
                    }
                }
            }
        }

        // Fit mixed acceptances
        for(int rw = 0; rw < 2; ++rw) {
            for(int sel = 0; sel < 2; ++sel) {
                for(int bkg = 0; bkg < 3; ++bkg) {
                    for(int det = 0; det < 3; ++det) {
                        // Use Be7 xpoints for mixed acceptance
                        auto& xpoints = Acc::xpointsEkfor10[det];
                        
                        std::string fitName = Form("f_acc_%s_%s_%s_%s_BeTotal",
                                (rw ? "reweight" : "orig"), selTypes[sel].c_str(),
                                bkgTypes[bkg].c_str(), detectorNames_[det].c_str());
                        
                        auto f_acc = SplineFit(h_acc_mix_[rw][sel][bkg][det], 
                                             xpoints.data(), xpoints.size(),
                                             0x38, "b2e2", fitName.c_str(), 
                                             xpoints[0], 650);
                        f_acc->Write();
                    }
                }
            }
        }
    }

    void SaveResults() {
        std::unique_ptr<TFile> outFile(TFile::Open((outputDir_ + "/BeTotal_analysis.root").c_str(), "RECREATE"));
        if(!outFile) throw std::runtime_error("Cannot create output file");
        std::cout << "Saving results to: " << outFile->GetName() << std::endl;

        // 定义写入操作
        auto writeHist = [](TH1F* h) { h->Write(); };

        // 保存所有直方图
        ForEachHistogram(h_num_, writeHist);
        for (auto& [mass, hists] : h_bkg_study_) {
            ForEachHistogram(hists, writeHist);
        }
        for (auto& [mass, hists] : h_acc_) {
            ForEachHistogram(hists, writeHist);
        }
        
        // 保存混合直方图
        ForEachHistogram(h_num_mix_, writeHist);
        ForEachHistogram(h_bkg_study_mix_, writeHist);
        ForEachHistogram(h_acc_mix_, writeHist);

        // 保存源通量函数
        for(const auto& [_, flux] : sourceFluxes_) {
            if(flux) flux->Write();
        }
    }

    void CleanupHistograms() {
        // 定义清理操作
        auto cleanupHist = [](TH1F*& h) { 
            delete h; 
            h = nullptr; 
        };

        // 清理所有直方图
        ForEachHistogram(h_num_, cleanupHist);
        for (auto& [mass, hists] : h_bkg_study_) {
            ForEachHistogram(hists, cleanupHist);
        }
        for (auto& [mass, hists] : h_acc_) {
            ForEachHistogram(hists, cleanupHist);
        }
        
        // 清理混合直方图
        ForEachHistogram(h_num_mix_, cleanupHist);
        ForEachHistogram(h_bkg_study_mix_, cleanupHist);
        ForEachHistogram(h_acc_mix_, cleanupHist);
    }
};

// Main function to run the analysis
void cal_beAcc() {
    try {
        BeBkgAnalyzer analyzer(
            "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MCBeryllium7.root",
            "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MCBeryllium9.root",
            "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MCBeryllium10.root",
            "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg"
        );
        analyzer.Process();
        std::cout << "Analysis completed successfully" << std::endl;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}