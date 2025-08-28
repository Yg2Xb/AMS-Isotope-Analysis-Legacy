/*********************************************************** 
 * Calculation of mixed acceptance for Beryllium isotopes.
 * Author: Yan Zixuan  Date: 20250328
 ***********************************************************/
#include "../Tool.h"
#include <TFile.h>
#include <TH1F.h>
#include <memory>
#include <map>
#include <array>
#include <stdexcept>

using namespace AMS_Iso;

class BeTotalAccCalculator {
	public:
		static constexpr std::array<int, 3> kMasses{7, 9, 10};
		static constexpr std::array<double, 3> kWeights{0.6, 0.3, 0.1};
		static constexpr int kNSel = 2;
		static constexpr int kNBkg = 3;
		static constexpr int kNDet = 3;

	private:
		const std::array<std::string, 3> detNames_{"TOF", "NaF", "AGL"};
		const std::array<std::string, 2> selTypes_{"Complete", "Loose"};
		const std::array<std::string, 3> bkgTypes_{"Bkg", "NoBkg", "StrictBkg"};

		// Input histograms
		std::map<int, std::vector<std::vector<TH1F*>>> h_input_;  // [mass][sel*3+bkg][det]
		std::map<int, double> n_events_;
		std::map<int, TH1F*> h_num_;  // Generated number histograms

		// Mixed histograms
		std::vector<std::vector<TH1F*>> h_mix_;  // [sel*3+bkg][det]
		TH1F* h_num_mix_;
		std::vector<std::vector<TH1F*>> h_acc_mix_;  // [sel*3+bkg][det]

	public:
		BeTotalAccCalculator() {
			Initialize();
			LoadHistograms();
		}

		~BeTotalAccCalculator() {
			Cleanup();
		}

		void Process() {
			CalculateBinnedNumbers();
			MixHistograms();
			CalculateAcceptance();
			SaveResults();
		}

	private:
		void Initialize() {
			for(int mass : kMasses) {
				h_input_[mass].resize(6);  // sel*3+bkg
				for(auto& det_vec : h_input_[mass]) {
					det_vec.resize(3, nullptr);  // det
				}
			}

			h_mix_.resize(6);
			for(auto& det_vec : h_mix_) {
				det_vec.resize(3, nullptr);
			}

			h_acc_mix_.resize(6);
			for(auto& det_vec : h_acc_mix_) {
				det_vec.resize(3, nullptr);
			}
		}

        void LoadHistograms() {
			const auto& ekBins = Binning::KineticEnergyBins[3][2];

            for (int mass : kMasses) {
                std::cout << "Processing Be" << mass << "..." << std::endl;
        
                // 打开文件
                std::string filename = Form("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/be%dacc.root", mass);
                std::unique_ptr<TFile> filein(TFile::Open(filename.c_str()));
                if (!filein) {
                    throw std::runtime_error("Cannot open filein: " + filename);
                }
        
                // 加载事件总数
                TH1F* h_num = (TH1F*)filein->Get("hMCnum");
                if (!h_num) {
                    throw std::runtime_error("Cannot find hMCnum in " + filename);
                }
                n_events_[mass] = h_num->GetBinContent(1);
        
                // 加载事件直方图
                for (int sel = 0; sel < kNSel; ++sel) {
                    for (int bkg = 0; bkg < kNBkg; ++bkg) {
                        for (int det = 0; det < kNDet; ++det) {
                            std::string histname = Form("Event_Ek_%s_Berlium_%d", 
                                                      detNames_[det].c_str(), sel * 3 + bkg);
                            
                            // 从文件中获取直方图
                            TH1F* h = (TH1F*)filein->Get(histname.c_str());
                            if (!h) {
                                throw std::runtime_error("Cannot find histogram: " + histname);
                            }
                            h_input_[mass][sel*3+bkg][det] = new TH1F((histname + "_" + std::to_string(mass)).c_str(),
                                    h->GetTitle(), Constants::RIGIDITY_BINS, ekBins.data());
                            h->Copy(*h_input_[mass][sel*3+bkg][det]);
                            h_input_[mass][sel*3+bkg][det]->SetDirectory(0);
                        }
                    }
                }
            }
        }

		void CalculateBinnedNumbers() {
			// Use Be7 binning for all
			const auto& ekBins = Binning::KineticEnergyBins[3][2];

			for(int mass : kMasses) {
				h_num_[mass] = new TH1F(Form("h_num_Be%d", mass), 
						Form("Be%d Number;Ek;Events", mass),
						Constants::RIGIDITY_BINS, ekBins.data());
				h_num_[mass]->Sumw2();

				std::unique_ptr<TF1> f_flux(new TF1("f_flux", "TMath::Power(x,-1)", 
							Acc::_xlow, Acc::_xup));
				double baseFluxIntegral = f_flux->Integral(Acc::_xlow, Acc::_xup);

				const double* bins = h_num_[mass]->GetXaxis()->GetXbins()->GetArray();
				for(int bin = 0; bin < h_num_[mass]->GetNbinsX(); ++bin) {
					double Rlow = betaToRigidity(kineticEnergyToBeta(bins[bin]), 4, mass, false);
					double Rup = betaToRigidity(kineticEnergyToBeta(bins[bin + 1]), 4, mass, false);
					h_num_[mass]->SetBinContent(bin + 1, 
							n_events_[mass] * f_flux->Integral(Rlow, Rup) / baseFluxIntegral);
				}
			}
		}

		void MixHistograms() {
			const auto& ekBins = Binning::KineticEnergyBins[3][2];
			// Mix generated numbers
			h_num_mix_ = new TH1F("h_num_BeTotal", "BeTotal Number;Ek;Events",
					Constants::RIGIDITY_BINS, ekBins.data());
			h_num_mix_->Sumw2();


			for(size_t i = 0; i < kMasses.size(); ++i) {
				h_num_mix_->Add(h_num_[kMasses[i]], kWeights[i]);
                printf("B%d num in bin22=%f, weight=%f, weightedN = %f\n",kMasses[i], h_num_[kMasses[i]]->GetBinContent(22), kWeights[i], h_num_[kMasses[i]]->GetBinContent(22)*kWeights[i]);
			}
            printf("Betotal weighted num in bin22:%f\n", h_num_mix_->GetBinContent(22));

			// Mix event histograms
			for(int sel = 0; sel < kNSel; ++sel) {
				for(int bkg = 0; bkg < kNBkg; ++bkg) {
					for(int det = 0; det < kNDet; ++det) {
						std::string name = Form("h_mix_%s_%s_%s_BeTotal",
								selTypes_[sel].c_str(), bkgTypes_[bkg].c_str(), 
								detNames_[det].c_str());

						h_mix_[sel*3+bkg][det] = new TH1F(name.c_str(), name.c_str(),
								Constants::RIGIDITY_BINS, ekBins.data());
						h_mix_[sel*3+bkg][det]->Sumw2();

						for(size_t i = 0; i < kMasses.size(); ++i) {
							h_mix_[sel*3+bkg][det]->Add(h_input_[kMasses[i]][sel*3+bkg][det], kWeights[i]);
                            if(sel == 0 && bkg == 0) printf("det%d, B%d events in bin22=%f, weight=%f, weightedN = %f\n",det, kMasses[i], h_input_[kMasses[i]][sel*3+bkg][det]->GetBinContent(22), kWeights[i], h_input_[kMasses[i]][sel*3+bkg][det]->GetBinContent(22)*kWeights[i] );
						}
                        if(sel == 0 && bkg == 0) printf("Betotal weighted events in bin22:%f\n", h_mix_[sel*3+bkg][det]->GetBinContent(22));
					}
				}
			}
		}

		void CalculateAcceptance() {
			const double scale_factor = 3.9 * 3.9 * TMath::Pi();

			for(int sel = 0; sel < kNSel; ++sel) {
				for(int bkg = 0; bkg < kNBkg; ++bkg) {
					for(int det = 0; det < kNDet; ++det) {
						std::string name = Form("h_acc_orig_%s_%s_%s_BeTotal",
								selTypes_[sel].c_str(), bkgTypes_[bkg].c_str(),
								detNames_[det].c_str());

						h_acc_mix_[sel*3+bkg][det] = (TH1F*)h_mix_[sel*3+bkg][det]->Clone(
								name.c_str());
						h_acc_mix_[sel*3+bkg][det]->Divide(h_num_mix_);
						h_acc_mix_[sel*3+bkg][det]->Scale(scale_factor);
					}
				}
			}
		}

		void SaveResults() {
			TFile* outFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/BeTotal_analysis.root", 
					"UPDATE");
			if(!outFile) throw std::runtime_error("Cannot open output file");

			// Delete existing BeTotal histograms
			TList* list = outFile->GetListOfKeys();
			TIter next(list);
			TKey* key;
			std::vector<std::string> toDelete;
			while((key = (TKey*)next())) {
				std::string name = key->GetName();
				if(name.find("BeTotal") != std::string::npos) {
					toDelete.push_back(name);
				}
			}

			std::cout << "Found " << toDelete.size() << " existing BeTotal histograms to delete" << std::endl;

			for(const auto& name : toDelete) {
				outFile->Delete((name + ";*").c_str());
			}

			// Save new acceptances with spline fits
			for(int sel = 0; sel < kNSel; ++sel) {
				for(int bkg = 0; bkg < kNBkg; ++bkg) {
					for(int det = 0; det < kNDet; ++det) {
						auto hist = h_acc_mix_[sel*3+bkg][det];

						// Fit acceptance
						auto& xpoints = Acc::xpointsEkfor7[det];  // Use Be7 points
						std::string fitName = Form("f_acc_orig_%s_%s_%s_BeTotal",
								selTypes_[sel].c_str(), bkgTypes_[bkg].c_str(),
								detNames_[det].c_str());

						auto f_acc = SplineFit(hist, xpoints.data(), xpoints.size(),
								0x38, "b2e2", fitName.c_str(), xpoints[0], 650);

						hist->Write();
					}
				}
			}

			outFile->Close();
		}

		void Cleanup() {
			delete h_num_mix_;
			for(auto& [_, hist] : h_num_) {
				delete hist;
			}
			for(auto& v1 : h_mix_) {
				for(auto& hist : v1) {
					delete hist;
				}
			}
			for(auto& v1 : h_acc_mix_) {
				for(auto& hist : v1) {
					delete hist;
				}
			}
			for(auto& [_, v1] : h_input_) {
				for(auto& v2 : v1) {
					for(auto& hist : v2) {
						delete hist;
					}
				}
			}
		}
};

void cal_BeTotalAcc() {
	try {
		BeTotalAccCalculator calculator;
		calculator.Process();
		std::cout << "BeTotal acceptance calculation completed successfully" << std::endl;
	} catch(const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}
}
