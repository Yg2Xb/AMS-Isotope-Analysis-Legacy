/***********************************************************
 *  File: acc_bkg.cpp
 *
 *  Analysis of background events from MC fragmentation.
 *  Processes events from MC simulation (e.g., O16 -> Be)
 *  and generates histograms for different isotopes and detectors.
 *
 *  Author: [Your name]
 *  Date: 20250308
 ***********************************************************/

#include "../Tool.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <memory>
#include <array>
#include <string>
#include <stdexcept>

using namespace AMS_Iso;

class BkgAnalyzer {
	public:
		// Constants for array dimensions
		static constexpr int kNDetectors = 3;    // TOF, NaF, AGL
		static constexpr int kNIsotopes = 3;     // Be7, Be9, Be10
		static constexpr int kNHistTypes = 5;    // Different histogram types

		// Constructor takes input file path and output directory
		BkgAnalyzer(const std::string& inputFile, const std::string& outputDir)
			: inputFile_(inputFile)
			  , outputDir_(outputDir)
			  , sourceFlux_(nullptr) {
				  Initialize();
			  }

		// Destructor ensures proper cleanup
		~BkgAnalyzer() {
			for (auto& det : histograms_) {
				for (auto& iso : det) {
					for (auto& hist : iso) {
						delete hist;
					}
				}
			}
			delete sourceFlux_;
		}

		// Main processing function
		void Process() {
			LoadSourceFlux();
			ProcessEvents();
			SaveResults();
		}

	private:
		// Member variables
		std::string inputFile_;
		std::string outputDir_;
		TF1* sourceFlux_;
		// Source flux fit points
		std::array<double, 13> kSourceFluxFitPoints{
			2., 2.6, 3.2, 4, 6, 13.0, 30, 50.0, 100.0, 211.0, 330.0, 643.0, 2000.0
		};

		// 3D array of histograms [detector][isotope][type]
		std::array<std::array<std::array<TH1F*, kNHistTypes>, kNIsotopes>, kNDetectors> histograms_;

		// Detector names for convenient access
		const std::array<std::string, kNDetectors> detectorNames_{"TOF", "NaF", "AGL"};

		// Histogram type names for convenient access
		const std::array<std::string, kNHistTypes> histTypes_{
			"NumEk", "genEk", "recEk", "genEk_weighted", "recEk_weighted"
		};

		void Initialize() {
			// Create histograms for each detector, isotope, and type
			for (int det = 0; det < kNDetectors; ++det) {
				for (int iso = 0; iso < kNIsotopes; ++iso) {
					// Get mass number for current isotope
					int mass = IsotopeData[3].mass_[iso];  // Z=4 (Be) is at index 3

					// Get kinetic energy bins for this isotope
					const auto& ekBins = Binning::KineticEnergyBins[3][iso];

					for (int type = 0; type < kNHistTypes; ++type) {
						// Create histogram name and title
						std::string name = Form("h_%s_%s_Be%d", 
							histTypes_[type].c_str(), 
							detectorNames_[det].c_str(), 
							mass);
						std::string xtitle = (type == 2 || type == 4) ? 
								detectorNames_[det] + "Ek" : "geneEk";
						std::string title = Form("%s;%s;Events", name.c_str(), xtitle.c_str());
						
						histograms_[det][iso][type] = new TH1F(name.c_str(), title.c_str(),
								Constants::RIGIDITY_BINS, ekBins.data());

						// Enable error calculation
						histograms_[det][iso][type]->Sumw2();
					}
				}
			}
		}

		void LoadSourceFlux() {
			// Open Source flux data file
			std::unique_ptr<TFile> fluxFile(TFile::Open(
						"/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_OxygenFlux.root"));
			if (!fluxFile || fluxFile->IsZombie()) {
				throw std::runtime_error("Cannot open Source flux file");
			}

			// Get the graph and create spline fit
			TGraphAsymmErrors* graph = (TGraphAsymmErrors*)fluxFile->Get("graph1");
			if (!graph) {
				throw std::runtime_error("Cannot find Source flux graph");
			}

			// Create spline fit of the Source flux using predefined points
			sourceFlux_ = SplineFit(static_cast<TGraph*>(graph), 
					kSourceFluxFitPoints.data(),
					kSourceFluxFitPoints.size(),
					0x38,
					"b1e1",
					"f_SourceFlux",
					1.,
					2000.0);
			if (!sourceFlux_) {
				throw std::runtime_error("Failed to create Source flux fit");
			}
		}

		double getSourceFluxNorm() const {
			static const double norm = sourceFlux_->Integral(1., 2000.);
			return norm;
		}

		double calculateWeight(double rigidity) const {
			if (rigidity <= 0) {
				return 0.0;
			}

			// Use existing MC functions from Tool.h
			double mcWeight = getMCFunction().Eval(rigidity) / getMCNorm();
			// Calculate Source flux weight
			double fluxWeight = sourceFlux_->Eval(rigidity) / getSourceFluxNorm();

			return fluxWeight / mcWeight;
		}

		void CalculateBinGeneNum(double N_gen, int isoIndex) {
			int mass = IsotopeData[3].mass_[isoIndex];

			// Create flux function for background calculation
			TF1 f_flux("f_flux", "TMath::Power(x,-1)", Acc::_xlow, Acc::_xup);
			double baseFluxIntegral = f_flux.Integral(Acc::_xlow, Acc::_xup);

			// Calculate background for each detector
			for (int det = 0; det < kNDetectors; ++det) {
				TH1F* hist = histograms_[det][isoIndex][0];  // bgEk histogram
				const double* bins = hist->GetXaxis()->GetXbins()->GetArray();

				// Fill each bin
				for (int bin = 0; bin < hist->GetNbinsX(); ++bin) {
					double Rlow = betaToRigidity(kineticEnergyToBeta(bins[bin]), 4, mass, false);
					double Rup = betaToRigidity(kineticEnergyToBeta(bins[bin + 1]), 4, mass, false);

					double N_gen_bin = N_gen * (f_flux.Integral(Rlow, Rup) / baseFluxIntegral);
					hist->SetBinContent(bin + 1, N_gen_bin);
					//hist->SetBinContent(bin + 1, N_gen_bin / (3.9 * 3.9 * M_PI));
				}
			}
		}

		void ProcessEvents() {
			// Open input file
			std::unique_ptr<TFile> file(TFile::Open(inputFile_.c_str()));
			if (!file || file->IsZombie()) {
				throw std::runtime_error("Cannot open input file: " + inputFile_);
			}

			// Get total MC events
			TH1F* hMCnum = (TH1F*)file->Get("hMCnum");
			if (!hMCnum) {
				throw std::runtime_error("Cannot find hMCnum histogram");
			}
			double N_gen = hMCnum->GetBinContent(1);

			// Calculate background histograms
			for (int iso = 0; iso < kNIsotopes; ++iso) {
				CalculateBinGeneNum(N_gen, iso);
			}

			// Get tree
			TTree* tree = (TTree*)file->Get("saveTree");
			if (!tree) {
				throw std::runtime_error("Cannot find saveTree");
			}

			// Set up branch addresses
			double InnerRig, generatedRig, generatedEk;
			double ToFBeta, NaFBeta, AGLBeta;
			double ToFEk, NaFEk, AGLEk;
			int DetectorAndBkgCut[3][3];
			tree->SetBranchAddress("InnerRig", &InnerRig);
			tree->SetBranchAddress("generatedRig", &generatedRig);
			tree->SetBranchAddress("generatedEk", &generatedEk);
			tree->SetBranchAddress("ToFBeta", &ToFBeta);
			tree->SetBranchAddress("NaFBeta", &NaFBeta);
			tree->SetBranchAddress("AGLBeta", &AGLBeta);
			tree->SetBranchAddress("ToFEk", &ToFEk);
			tree->SetBranchAddress("NaFEk", &NaFEk);
			tree->SetBranchAddress("AGLEk", &AGLEk);
			tree->SetBranchAddress("DetectorAndBkgCut", DetectorAndBkgCut);

			// Process all events
			Long64_t nEntries = tree->GetEntries();
			for (Long64_t entry = 0; entry < nEntries; ++entry) {
				tree->GetEntry(entry);

				// Arrays for convenient access to beta and Ek values
				const std::array<double, 3> betas{ToFBeta, NaFBeta, AGLBeta};
				const std::array<double, 3> recoEks{ToFEk, NaFEk, AGLEk};

				// Fill histograms for each detector and isotope
				for (int det = 0; det < kNDetectors; ++det) {
					for (int iso = 0; iso < kNIsotopes; ++iso) {
						if (DetectorAndBkgCut[det][iso] != 1) continue;

						// Calculate weight using Tool.h functions and Source flux
						double weight = calculateWeight(generatedRig);

						// Fill unweighted histograms
						histograms_[det][iso][1]->Fill(generatedEk);  // genEk
						histograms_[det][iso][2]->Fill(recoEks[det]); // recEk

						// Fill weighted histograms
						histograms_[det][iso][3]->Fill(generatedEk, weight);  // genEk_weighted
						histograms_[det][iso][4]->Fill(recoEks[det], weight); // recEk_weighted
					}
				}
			}
		}

		void SaveResults() {
			// Create output file
			std::string outFileName = outputDir_ + "/Oxygen16_bkg_analysis.root";
			std::unique_ptr<TFile> outFile(TFile::Open(outFileName.c_str(), "RECREATE"));
			if (!outFile || outFile->IsZombie()) {
				throw std::runtime_error("Cannot create output file: " + outFileName);
			}

			// Save all histograms
			for (const auto& det : histograms_) {
				for (const auto& iso : det) {
					for (const auto& hist : iso) {
						hist->Write();
					}
				}
			}

			// Save Source flux function
			sourceFlux_->Write();
		}
};

// Main function
void cal_bkghist() {
	try {
		std::string inputFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/MCO16.root";
		std::string outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg";

		BkgAnalyzer analyzer(inputFile, outputDir);
		analyzer.Process();

		std::cout << "Background analysis completed successfully" << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}
}
