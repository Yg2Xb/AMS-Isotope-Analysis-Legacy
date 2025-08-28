#include "../Tool.h"
using namespace AMS_Iso;

class AcceptanceAnalyzer {
	public:
		AcceptanceAnalyzer(int z) : Z_(z) {
			const auto& isotope = IsotopeData[Z_ - 1];
			h_bgEk_ = new TH1F("h_bgEk", "N_gen", 
					Constants::RIGIDITY_BINS,
					Binning::KineticEnergyBins[Z_ - 1][0].data());;
		}

		~AcceptanceAnalyzer() {
			if (h_bgEk_) delete h_bgEk_;
		}

		void ProcessIsotope(int mass, const std::string& inputDir, const std::string& outputDir) {
			const auto& isotope = IsotopeData[Z_ - 1];
			int massIndex = GetMassIndex(mass);
			if (massIndex < 0) {
				std::cerr << "Invalid mass number " << mass << " for Z=" << Z_ << std::endl;
				return;
			}

			std::string infile = inputDir + ConvertElementName(isotope.name_, false) + std::to_string(mass) + ".root";
			std::string outfile = outputDir + ConvertElementName(isotope.name_, false) + std::to_string(mass) + "_acc";

			std::unique_ptr<TFile> file(TFile::Open(infile.c_str()));
			if (!file || file->IsZombie()) {
				std::cerr << "Could not open " << infile << std::endl;
				return;
			}

			TH1F* hMCnum = (TH1F*)file->Get("hMCnum");
			if (!hMCnum) {
				std::cerr << "Could not find hMCnum in " << infile << std::endl;
				return;
			}
			cout<<"break?"<<endl;
			h_bgEk_->SetBins(Constants::RIGIDITY_BINS,
					Binning::KineticEnergyBins[Z_ - 1][massIndex].data());
			cout<<"check bg "<<h_bgEk_->GetBinLowEdge(10)<<endl;
			h_bgEk_->Sumw2();
			cout<<"break?"<<endl;

			CalculateBackground(hMCnum->GetBinContent(1), mass);
			cout<<"break?"<<endl;
			ProcessDetectors(file.get(), outfile, isotope.color_[massIndex]);
		}

	private:
		int Z_;
		TH1F* h_bgEk_;

		int GetMassIndex(int mass) const {
			const auto& isotope = IsotopeData[Z_ - 1];
			for (size_t i = 0; i < isotope.mass_.size(); ++i) {
				if (isotope.mass_[i] == mass) return i;
			}
			return -1;
		}

		void CalculateBackground(double N_gen, int mass) {
			TF1 f_flux("f_flux", "TMath::Power(x,-1)", Acc::_xlow, Acc::_xup);
			double base_fluxIntegral = f_flux.Integral(Acc::_xlow, Acc::_xup);

			const auto& bins = h_bgEk_->GetXaxis()->GetXbins()->GetArray();
			for (int j = 0; j < h_bgEk_->GetNbinsX(); ++j) {
				double Rlow = betaToRigidity(kineticEnergyToBeta(bins[j]), Z_, mass, false);
				double Rup = betaToRigidity(kineticEnergyToBeta(bins[j + 1]), Z_, mass, false);
				//double Rlow = Binning::RigidityBins[j];
				cout<<Rlow<<endl;
				//double Rup = Binning::RigidityBins[j+1];
				double N_gen_bin = N_gen * (f_flux.Integral(Rlow, Rup) / base_fluxIntegral);
				h_bgEk_->SetBinContent(j + 1, N_gen_bin/(3.9 * 3.9 * TMath::Pi()));
			}
		}

		void ProcessDetectors(TFile* file, const std::string& outfile, int color) {
			std::unique_ptr<TFile> outRoot(TFile::Open((outfile + ".root").c_str(), "RECREATE"));

			for (size_t detIdx = 0; detIdx < Detector::BetaTypes.size(); ++detIdx) {
				const auto& det = Detector::BetaTypes[detIdx];
				std::vector<TH1D*> cuts_hists;

				const auto& cuts = (det.name_ == "TOF") ? 
					Detector::TOFCuts : Detector::RICHCuts;

				ProcessDetectorCuts(file, det, cuts, color, cuts_hists, outRoot.get());

				if (!cuts_hists.empty()) {
					DrawDetectorAcceptance(cuts_hists, det, outfile);
					for (auto hist : cuts_hists) delete hist;
				}
			}
		}

		void ProcessDetectorCuts(TFile* file, const BetaConfig& det,
				const std::array<std::string, 10>& cuts,
				int color, std::vector<TH1D*>& cuts_hists,
				TFile* outRoot) {
			const auto& isotope = IsotopeData[Z_ - 1];
			std::string prefix = isotope.name_;

			std::string filename = file->GetName();
			bool isBe7 = (filename.find("Be7") != std::string::npos);
			bool isBe9 = (filename.find("Be9") != std::string::npos);
			bool isBe10 = (filename.find("Be10") != std::string::npos);

			for (size_t i = 0; i < cuts.size(); ++i) {
				std::string detName = det.name_ ;
				std::string histName = Form("Event_Ek_%s_%s_%zu", 
						detName.c_str(), prefix.c_str(), i);
				TH1D* hist = (TH1D*)file->Get(histName.c_str());
				if (!hist) {cout<<"No hist!"<<endl; continue;}

				TH1D* acc = (TH1D*)hist->Clone(Form("acc_%s", histName.c_str()));
				acc->Divide(hist, h_bgEk_, 1, 1, "B");

				size_t detIndex = (&det - &Detector::BetaTypes[0]);
				// 根据是否为Be10选择不同的xpoints
				auto& xpoints = isBe7 ? Acc::xpointsEkfor7[detIndex] : 
                                      isBe9 ? Acc::xpointsEkfor9[detIndex] : 
                                              Acc::xpointsEkfor10[detIndex];
                cout<<"begin fit:"<<isBe7<<" "<<filename<<" "<<xpoints[0]<<endl;
				auto f_acc = SplineFit(acc, xpoints.data(), 
						xpoints.size(), 0x38, "b2e2", 
						"f_acc", xpoints[0], 650);
				f_acc->SetLineColor(color);
				f_acc->SetLineWidth(2);

				cuts_hists.push_back(acc);
				outRoot->cd();
				acc->Write();
				int last = (det.name_ == "TOF") ? 8 : 9;
				if(i == last)
				{
					f_acc->SetName(Form("%sAccSplineFit", det.name_.c_str()));
					f_acc->SetTitle(Form("%sAccSplineFit;Ek/n[GeV/n];Acceptance(m_{-2} s_{-1})", det.name_.c_str()));
					f_acc->Write();
				} 
			}
		}

		void DrawDetectorAcceptance(const std::vector<TH1D*>& hists, 
				const BetaConfig& det,
				const std::string& outfile) {
			std::unique_ptr<TCanvas> c(new TCanvas("c", "", 800, 600));
			c->SetLogx();

			for (size_t i = 0; i < hists.size(); ++i) {
				hists[i]->SetLineColor(colors[i]);
				hists[i]->SetMarkerColor(colors[i]);
				hists[i]->SetMarkerStyle(20);
				hists[i]->GetXaxis()->SetTitle("Kinetic Energy [GeV/n]");
				hists[i]->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");

				if (i == 0) {
					hists[i]->Draw("P");
					hists[i]->GetYaxis()->SetRangeUser(0,0.82);
				} else {
					hists[i]->Draw("P SAME");
				}
			}

			c->SaveAs((outfile + "_" + det.name_ + ".pdf").c_str());
		}
};

void DrawCombinedAcceptance(int Z, const std::string& outputDir, int last) {
	const auto& isotope = IsotopeData[Z - 1];
	std::unique_ptr<TCanvas> c(new TCanvas("c", "", 800, 600));
	c->SetLogx();

	std::unique_ptr<TLegend> leg(new TLegend(0.24, 0.37, 0.36, 0.53));
	leg->SetBorderSize(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.04);

	const std::map<std::string, int> lastCuts = {
		{"TOF", 8}, {"NaF", 9}, {"AGL", 9}
	};

	std::vector<TH1D*> histograms;
	bool first = true;

	for (size_t i = 0; i < isotope.mass_.size(); ++i) {
		int mass = isotope.mass_[i];
		for (const auto& det : Detector::BetaTypes) {
			if(last <= 7 && (det.name_ != "TOF")) continue;
			std::string filename = outputDir + ConvertElementName(isotope.name_, false) + 
				std::to_string(mass) + "_acc.root";
			std::cout << "Looking for file: " << filename << std::endl;
			TFile *f = TFile::Open(filename.c_str());
			if (!f) continue;

			std::string detName = det.name_ ;
			std::string histname = Form("acc_Event_Ek_%s_%s_%d", 
					detName.c_str(), isotope.name_.c_str(), lastCuts.at(det.name_));

			std::cout << "Looking for histogram: " << histname << std::endl;
			TH1D* hist = (TH1D*)f->Get(histname.c_str());
			if (!hist || hist->GetEntries() == 0) {cout<<"No combined hist!"<<endl; continue;}

			hist = (TH1D*)hist->Clone();
			histograms.push_back(hist);
			cout<<hist->GetBinLowEdge(10)<<endl;

			hist->SetLineColor(isotope.color_[i]);
			hist->SetMarkerColor(isotope.color_[i]);
			hist->SetMarkerStyle(20);
			hist->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");
			hist->GetXaxis()->SetTitle("generated Ek [GeV/n]");
			hist->SetTitle("MC Acceptance");
			hist->GetXaxis()->SetRangeUser(0.2, 200);
			hist->SetMinimum(0.002);

			if (first) {
				hist->Draw("P");
				first = false;
			} else {
				hist->Draw("P SAME");
			}

			if (det.name_ == "TOF") {
				leg->AddEntry(hist, Form("%s%d", ConvertElementName(isotope.name_, false).c_str(), mass), "p");
			}
		}
	}

	leg->Draw();
	c->Update();
	c->Print((outputDir + ConvertElementName(isotope.name_, false) + Form("_Compare%d.pdf",last)).c_str());

	TFile *fcanvas = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Acc/acc_canvas.root","recreate");
	fcanvas->cd();
	c->Write();
	fcanvas->Close();
	
	for (auto hist : histograms) delete hist;
}

void acc() {
	cout<<"begin cal acc"<<endl;
	const std::string inputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/";
	const std::string outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Acc/";

	AcceptanceAnalyzer analyzer(4);
	for (int mass : IsotopeData[3].mass_) {
		analyzer.ProcessIsotope(mass, inputDir, outputDir);
	}
    /*
	for(int i=0;i<8;i++){
		DrawCombinedAcceptance(4, outputDir, i);
	}
        */
	DrawCombinedAcceptance(4, outputDir, 10);
}
