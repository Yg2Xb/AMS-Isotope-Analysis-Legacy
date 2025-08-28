/***********************************************************
 *  Mass Template Generation for Isotope Analysis cpp file
 *  Author: Z.Yan
 *  Date: 2024.10.31
 *  Purpose: select pure be10 data, study data and mc consisitancy
 ***********************************************************/

#include "FillTempHist.h"
using namespace AMS_Iso;

struct DetectorConfig {
    const char* name;
    double beta_min;
    double beta_max; 
    double cutoff_max;
    double safety;
    unsigned int mask;
};

const std::vector<DetectorConfig> DETECTOR_CONFIGS = {
    {"ToF", 0.6, 0.9, 8, 1.01, TRACKER_TOF_MASK},
    {"NaF", 0.75, 1.001, 20, 1.005, TRACKER_NAF_MASK},
    {"AGL", 0.954, 1.001, 25, 1.0005, TRACKER_AGL_MASK}
};

class PureStudyHists {
public:
    TH2F *h_all, *h_pure, *h_invMass_rrc;
    TH1F *h_RRc_all, *h_RRc_pure, *h_RRc_highMass;

    PureStudyHists(const DetectorConfig& config) {
        h_all = new TH2F(Form("h%s", config.name), 
                        Form("Pass %s Selection;cutoff rigidity[GV];%s measured #beta", config.name, config.name), 
                        250, 0, config.cutoff_max, 300, config.beta_min, config.beta_max);
        
        h_pure = new TH2F(Form("h%s_pure", config.name),
                         Form("Pass %s and Pure Be10 Selection;cutoff rigidity[GV];%s measured #beta", config.name, config.name),
                         250, 0, config.cutoff_max, 300, config.beta_min, config.beta_max);
        
        h_RRc_all = new TH1F(Form("h%s_RRc_all", config.name),
                            Form("%s R/Rc Distribution;InnerRig/Cutoff Rig;Counts", config.name),
                            100, 0, 2);
        
        h_RRc_pure = new TH1F(Form("h%s_RRc_pure", config.name),
                             Form("%s R/Rc Pure Selection;InnerRig/Cutoff Rig;Counts", config.name),
                             100, 0, 2);
        
        h_RRc_highMass = new TH1F(Form("h%s_RRc_highMass", config.name),
                                 Form("%s R/Rc High Mass;InnerRig/Cutoff Rig;Counts", config.name),
                                 100, 0, 2);

        h_invMass_rrc = new TH2F(Form("h%s_invMass_rrc", config.name),
                                Form("%s 1/Mass vs R/Rc;1/Mass;InnerRig/Cutoff Rig", config.name),
                                100, 0.06, 0.22, 100, 0, 2);
        
        SetHistStyle();
    }

    void SetHistStyle() {
        const std::vector<std::pair<TH1*, int>> hist_colors = {
            {h_RRc_all, kBlue}, {h_RRc_pure, kRed}, {h_RRc_highMass, kGreen+2}
        };
        
        for(const auto& [hist, color] : hist_colors) {
            hist->SetLineColor(color);
            hist->SetMarkerColor(color);
            hist->SetMarkerStyle(20);
            hist->SetMarkerSize(0.8);
        }

		h_RRc_highMass->SetMarkerStyle(22);
    }

    void Write() const {
        for(auto* hist : {h_RRc_all, h_RRc_pure, h_RRc_highMass}) {
            hist->Write();
        }
        for(auto* hist : {h_all, h_pure, h_invMass_rrc}) {
            hist->Write();
        }
    }

    void SavePlots(const char* detector) const {
        SaveBetaCutoffPlots(detector);
        SaveRRcDistributions(detector);
        SaveInvMassRRcPlot(detector);
    }

private:
    void SaveBetaCutoffPlots(const char* detector) const {
        TCanvas c("c1", Form("%s Beta vs Cutoff", detector), 1200, 900);
        c.SetRightMargin(0.16);
        c.SetLogz();
        
        h_all->Draw("colz");
        DrawTheoryLines(&c, h_all->GetYaxis()->GetXmin(), h_all->GetYaxis()->GetXmax());
        c.SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/cf30/%s_beta_cutoff_cf30.png", detector));
        
        h_pure->GetZaxis()->SetRangeUser(1, h_all->GetMaximum());
        h_pure->Draw("colz");
        DrawTheoryLines(&c, h_pure->GetYaxis()->GetXmin(), h_pure->GetYaxis()->GetXmax());
        c.SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/cf30/%s_beta_cutoff_pure_cf30.png", detector));
    }

    void SaveRRcDistributions(const char* detector) const {
        TCanvas c("c2", Form("%s R/Rc Distributions", detector), 800, 600);
        
        h_RRc_pure->Draw("P");
        h_RRc_highMass->Draw("P SAME");
        
        TLegend leg(0.65, 0.25, 0.85, 0.45);
        leg.AddEntry(h_RRc_all, "All events", "p");
        leg.AddEntry(h_RRc_pure, "Pure Be10", "p");
        leg.AddEntry(h_RRc_highMass, "1/Mass > 0.12", "p");
        leg.Draw();
        
        c.SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/cf30/%s_RRc_distributions_cf30.png", detector));
    }

    void SaveInvMassRRcPlot(const char* detector) const {
        TCanvas c("c3", Form("%s InvMass vs R/Rc", detector), 800, 600);
        c.SetLogz(0);
		c.SetGrid();
        
        h_invMass_rrc->Draw("colz");
		c.SetRightMargin(0.16);
        c.SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/cf30/%s_invMass_rrc_cf30.png", detector));
    }
};

class MassAnalyzer {
public:
    MassAnalyzer(TTree* tree, const IsotopeConfig& config) 
        : tree_(tree), config_(config) {
        SetupBranches();
        SetupHistograms();
    }

    void Process() {
        const Long64_t nEntries = tree_->GetEntries();
        Long64_t nProcessed = 0;

        for (Long64_t i = 0; i < nEntries; ++i) {
            if (i % 100000 == 0) {
                std::cout << "Processing " << i << "/" << nEntries 
                         << " (processed: " << nProcessed << ")" << std::endl;
            }

            tree_->GetEntry(i);
            ProcessEntry();
            nProcessed++;
        }

        std::cout << "Processed " << nProcessed << "/" << nEntries << " entries" << std::endl;
    }

    void SaveResults(const char* outputFileName) {
        auto outputFile = std::make_unique<TFile>(outputFileName, "RECREATE");
        
        for (size_t i = 0; i < h_inv_mass_.size(); ++i) {
            for (auto& hist : h_inv_mass_[i]) {
                hist->GetYaxis()->SetRangeUser(0.2, 10);
                hist->Write();
            }
        }

        for (const auto& study : detector_studies_) {
            study.Write();
            study.SavePlots(DETECTOR_CONFIGS[&study - &detector_studies_[0]].name);
        }

        SaveMassDistributionsPDF();
    }

private:
    void SetupBranches() {
        tree_->SetBranchAddress("InnerRig", &InnerRig_);
        tree_->SetBranchAddress("richBeta", &richBeta_);
        tree_->SetBranchAddress("tofBeta", &tofBeta_);
        tree_->SetBranchAddress("cutOffRig", &cutOffRig_);
        tree_->SetBranchAddress("cutOffRig25", &cutOffRig25_);
        tree_->SetBranchAddress("cutOffRig40", &cutOffRig40_);
        tree_->SetBranchAddress("rich_NaF", &richNaF_);
        tree_->SetBranchAddress("cutStatus", &cutStatus_);
        tree_->SetBranchAddress("ntrack", &ntrack_);
        tree_->SetBranchAddress("btstat", &btstat_);
    }

    void SetupHistograms() {
        SetupMassTemplates();
        for (const auto& config : DETECTOR_CONFIGS) {
            detector_studies_.emplace_back(config);
        }
    }

    void SetupMassTemplates() {
        h_inv_mass_.resize(config_.masses.size());
        for (size_t j = 0; j < config_.masses.size(); ++j) {
            const auto& ek_bins = getKineticEnergyBins(config_.charge, config_.masses[j]);
            h_inv_mass_[j].resize(3);  // ToF, NaF, AGL

            const std::vector<std::string> types = {"ToF", "NaF", "AGL"};
            for (size_t i = 0; i < types.size(); ++i) {
                h_inv_mass_[j][i] = std::make_unique<TH2D>(
                    Form("h_inv_%sMass_%sUseMass%d", types[i].c_str(), 
                         config_.dataName.c_str(), config_.masses[j]),
                    Form("inv_%sMass", types[i].c_str()),
                    200, 0.0, 0.5, ek_bins.size()-1, ek_bins.data());
            }
        }
    }

    void ProcessEntry() {
        for (size_t i = 0; i < DETECTOR_CONFIGS.size(); ++i) {
            const auto& config = DETECTOR_CONFIGS[i];
            if (ShouldProcessDetector(config)) {
                ProcessDetector(i);
            }
        }
    }

	bool ShouldProcessDetector(const DetectorConfig& config) {
		if (strcmp(config.name, "ToF") == 0) {
			return (cutStatus_ & TRACKER_TOF_MASK) == TRACKER_TOF_MASK && 
				   isValidBeta(tofBeta_);
		}
		else if (strcmp(config.name, "NaF") == 0) {
			return richNaF_ && 
				   (cutStatus_ & TRACKER_NAF_MASK) == TRACKER_NAF_MASK && 
				   isValidBeta(richBeta_);
		}
		else if (strcmp(config.name, "AGL") == 0) {
			return !richNaF_ && 
				   (cutStatus_ & TRACKER_AGL_MASK) == TRACKER_AGL_MASK && 
				   isValidBeta(richBeta_);
		}
		return false;
	}

    void ProcessDetector(size_t detector_idx) {
        const auto& config = DETECTOR_CONFIGS[detector_idx];
        const double beta = (detector_idx == 0) ? tofBeta_ : richBeta_;
        const double rrc = InnerRig_ / cutOffRig25_;

        auto& study = detector_studies_[detector_idx];
        study.h_all->Fill(cutOffRig25_, beta);
        study.h_RRc_all->Fill(rrc);

        auto result = calculateMass(beta, 1.0, InnerRig_, config_.charge);
        if (result.isValid() && selectPure3(beta, cutOffRig_, cutOffRig25_, cutOffRig40_, config.safety, 
                                         config_.charge, InnerRig_, false)) {
            study.h_pure->Fill(cutOffRig25_, beta);
            study.h_RRc_pure->Fill(rrc);
            study.h_invMass_rrc->Fill(result.invMass, rrc);

            if (result.invMass > 0.12) {
                study.h_RRc_highMass->Fill(rrc);
            }

            for (size_t j = 2; j < config_.masses.size(); ++j) {
                h_inv_mass_[j][detector_idx]->Fill(result.invMass, 
                                                 betaToKineticEnergy(beta));
            }
        }
    }

    void SaveMassDistributionsPDF() {
        TCanvas c("cMass", "Mass Distributions", 1200, 900);
        TString pdfName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/cf30/mass_distributions_cf30.pdf";
        c.Print(pdfName + "[");
        int beginBin[] = {6, 11, 20};
        int NBin[] = {11, 21, 22};
        for (size_t det_idx = 0; det_idx < DETECTOR_CONFIGS.size(); ++det_idx) {
            const auto& hist = h_inv_mass_[2][det_idx];
            for (int i = beginBin[det_idx]; i < beginBin[det_idx] + NBin[det_idx]; i += NBin[det_idx]) {
                auto* projX = hist->ProjectionX(
                    Form("px_%s_%d", DETECTOR_CONFIGS[det_idx].name, i), i, i+NBin[det_idx]-1);

                double ek_min = hist->GetYaxis()->GetBinLowEdge(i);
                double ek_max = hist->GetYaxis()->GetBinUpEdge(i+NBin[det_idx]-1);

                projX->GetXaxis()->SetTitle("1/mass");
                projX->GetYaxis()->SetTitle("Counts");
                projX->GetXaxis()->SetRangeUser(0.07, 0.21);
                projX->SetTitle(Form("%s %.2f-%.2f GeV/n", 
                                   DETECTOR_CONFIGS[det_idx].name, ek_min, ek_max));
                
                projX->Draw("E");
                c.Print(pdfName);
                delete projX;
            }
        }
        c.Print(pdfName + "]");
    }

    TTree* tree_;
    const IsotopeConfig& config_;
    
    // Branch variables
    double InnerRig_, richBeta_, tofBeta_, cutOffRig_, cutOffRig25_, cutOffRig40_;
    bool richNaF_;
    unsigned int cutStatus_;
    int ntrack_, btstat_;

    // Histograms
    std::vector<std::vector<std::unique_ptr<TH2D>>> h_inv_mass_;
    std::vector<PureStudyHists> detector_studies_;
};

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

    MassAnalyzer analyzer(tree, config);
    analyzer.Process();
    analyzer.SaveResults(outputFileName);
}

void BuildTempHist(const std::string& isotype) {
    const std::string baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/";
    try {
        if (isotope_configs.find(isotype) == isotope_configs.end()) {
            std::cerr << "Unknown isotope type: " << isotype << std::endl;
            return;
        }

        const auto& config = isotope_configs.at(isotype);
        std::string issFile = baseDir + config.dataName + 
                             std::to_string(config.masses[2]) + "_ntrk.root";
        std::string outFile = baseDir + config.dataName + "_temp_pure10_cf30.root";

        auto file = std::make_unique<TFile>(issFile.c_str());
        if (file && !file->IsZombie()) {
            ProcessTreeWithMass(file.get(), "saveTree", outFile.c_str(), false, config);
        } else {
            std::cerr << "Failed to open ISS file: " << issFile << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in BuildTempHist: " << e.what() << std::endl;
    }
}

void tempValid() {
    BuildTempHist("Be");
}