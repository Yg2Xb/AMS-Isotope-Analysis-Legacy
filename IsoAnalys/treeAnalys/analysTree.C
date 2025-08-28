#define analysTree_cxx
#include "/afs/cern.ch/work/z/zuhao/public/yanzx/IsoAnalys/treeAnalys/analysTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <vector>
#include <TString.h>
#include <memory>

struct HistogramSet {
	std::vector<std::vector<TH2D*>> h_rich_beta_theta;
	std::vector<std::vector<TH2D*>> h_rich_beta_phi;

	TH1F* h_tof_chist;
	TH1F* h_tof_chisc;
	TH1F* h_rich_pb;
	TH1F* h_rich_hit;
	TH1F* h_rich_usedm;
	TH1F* h_rich_npe0;
	TH1F* h_rich_npe2;
	TH1F* h_rich_npe_ratio;
	TH1F* h_rich_pmt;
};

void SaveAllHistograms(const HistogramSet& histset, const std::string& sampleName) {
	std::string outpath = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/analysTree/";
	std::string pdfname = outpath + "analysis_" + sampleName + ".pdf";

	TCanvas* c = new TCanvas("c", "", 800, 600);
	c->SetRightMargin(0.2);

	c->Print((pdfname + "[").c_str());
	c->SetLogz(); 
	const char* det_names[2] = {"NaF", "AGL"};
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 4; j++) {
			c->cd();
			TH2D* h_theta = histset.h_rich_beta_theta[i][j];
			h_theta->Draw("COLZ");

			if (strcmp(det_names[i], "AGL") == 0) {
				h_theta->SetTitle(TString(h_theta->GetTitle()) + " InnerRig > 100GV");
				h_theta->GetXaxis()->SetRangeUser(2.7, 3.15);
				h_theta->GetYaxis()->SetRangeUser(0.997, 1.003);
			} else {
				h_theta->RebinX(2);
				h_theta->RebinY(2);
				h_theta->SetTitle(TString(h_theta->GetTitle()) + " InnerRig > 50GV");
				h_theta->GetXaxis()->SetRangeUser(2.7, 3.15);
				h_theta->GetYaxis()->SetRangeUser(0.99, 1.02);
			}
			c->Print(pdfname.c_str());

			c->cd();
			TH2D* h_phi = histset.h_rich_beta_phi[i][j];
			h_phi->Draw("COLZ");

			if (strcmp(det_names[i], "AGL") == 0) {
				h_phi->SetTitle(TString(h_phi->GetTitle()) + " InnerRig > 100GV");
				h_phi->GetYaxis()->SetRangeUser(0.997, 1.003);
			} else {
				h_phi->RebinX(2);
				h_phi->RebinY(2);
				h_phi->SetTitle(TString(h_phi->GetTitle()) + " InnerRig > 50GV");
				h_phi->GetYaxis()->SetRangeUser(0.993, 1.007);
			}
			c->Print(pdfname.c_str());
		}
	}
	c->SetLogz(0); 

	TH1F* h1d_array[] = {
		histset.h_tof_chist,
		histset.h_tof_chisc,
		histset.h_rich_pb,
		histset.h_rich_hit,
		histset.h_rich_usedm,
		histset.h_rich_npe0,
		histset.h_rich_npe2,
		histset.h_rich_npe_ratio,
		histset.h_rich_pmt
	};
	int lg = 0;
	for(auto hist : h1d_array) {
		if(lg<2)
		{
			c->SetLogy();
		}
		else{
			c->SetLogy(0);
		}
		c->cd();
		hist->Draw();
		c->Print(pdfname.c_str());
		lg++;
	}

	c->Print((pdfname + "]").c_str());
	delete c;
}

void InitHistograms(HistogramSet& histset, const std::string& sampleName) {
	const char* det_names[2] = {"NaF", "AGL"};
	const char* cut_names[4] = {"NoCut", "Geo", "Basic", "Charge"};

	histset.h_rich_beta_theta.resize(2);
	histset.h_rich_beta_phi.resize(2);

	for(int i = 0; i < 2; i++) {
		histset.h_rich_beta_theta[i].resize(4);
		histset.h_rich_beta_phi[i].resize(4);

		for(int j = 0; j < 4; j++) {
			histset.h_rich_beta_theta[i][j] = new TH2D(
					Form("RBeta_Theta_%s_%s_%s", det_names[i], cut_names[j], sampleName.c_str()),
					Form("%s Cut%s %s;rich_theta;1/beta_{%s};Events", 
						det_names[i], cut_names[j], sampleName.c_str(), det_names[i]),
					1000, 2.7, 3.2,
					2000, 0.95, 1.05
					);

			histset.h_rich_beta_phi[i][j] = new TH2D(
					Form("RBeta_Phi_%s_%s_%s", det_names[i], cut_names[j], sampleName.c_str()),
					Form("%s Cut%s %s;rich_phi;1/beta_{%s};Events", 
						det_names[i], cut_names[j], sampleName.c_str(), det_names[i]),
					1000, -3.2, 3.2,
					2000, 0.95, 1.05
					);
		}
	}

	histset.h_tof_chist = new TH1F(Form("tof_chist_%s", sampleName.c_str()),
			"TOF Time Chi2;#chi^{2}_{t};Events", 2000, 0, 20);
	histset.h_tof_chisc = new TH1F(Form("tof_chisc_%s", sampleName.c_str()),
			"TOF Coordinate Chi2;#chi^{2}_{c};Events", 2000, 0, 20);
	histset.h_rich_pb = new TH1F(Form("rich_pb_%s", sampleName.c_str()),
			"RICH Kolmogrov Probability;rich_pb;Events", 1000, 0, 1);
	histset.h_rich_hit = new TH1F(Form("rich_hit_%s", sampleName.c_str()),
			"RICH Hits;rich_hit;Events", 100, 0, 100);
	histset.h_rich_usedm = new TH1F(Form("rich_usedm_%s", sampleName.c_str()),
			"RICH Reflected Hits;rich_usedm;Events", 50, 0, 50);
	histset.h_rich_npe0 = new TH1F(Form("rich_npe0_%s", sampleName.c_str()),
			"RICH Number of Photoelectrons in Ring;rich_npe[0];Events", 1000, 0, 100);
	histset.h_rich_npe2 = new TH1F(Form("rich_npe2_%s", sampleName.c_str()),
			"RICH Total Number of Photoelectrons;rich_npe[2];Events", 1000, 0, 100);
	histset.h_rich_npe_ratio = new TH1F(Form("rich_npe_ratio_%s", sampleName.c_str()),
			"RICH NPE[0]/NPE[2];rich_npe[0]/rich_npe[2];Events", 1000, 0, 1);
	histset.h_rich_pmt = new TH1F(Form("rich_pmt_%s", sampleName.c_str()),
			"RICH PMTs;PMTs;Events", 100, 0, 100);
}

void ProcessFile(const std::string& filename) {
	TFile *f = TFile::Open(filename.c_str());
	if(!f || f->IsZombie()) {
		std::cout << "Cannot open file: " << filename << std::endl;
		return;
	}

	TTree* tree = (TTree*)f->Get("saveTree");
	if(!tree) {
		std::cout << "Cannot find saveTree in file: " << filename << std::endl;
		f->Close();
		return;
	}

	size_t lastSlash = filename.find_last_of('/');
	size_t lastDot = filename.find_last_of('.');
	std::string sampleName = filename.substr(lastSlash + 1, lastDot - lastSlash - 1);

	analysTree t(tree);

	HistogramSet histset;
	InitHistograms(histset, sampleName);

	std::string outpath = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/analysTree/";
	TFile* outfile = new TFile((outpath + "analysis_" + sampleName + ".root").c_str(), "RECREATE");

	Long64_t nentries = tree->GetEntriesFast();
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		t.GetEntry(jentry);
		if (jentry % 100000 == 0) {
			std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;
		}

		bool trackerCut = (t.cutStatus & 1);
		bool richCut = (t.cutStatus & 4);
		bool cutoffCut = (t.cutStatus & 8);

		// Fill 1D histograms only if passing tracker cut and having positive rigidity
		if(trackerCut && t.InnerRig > 0) {
			// TOF quality histograms
			if(t.tofBeta > 0){
				histset.h_tof_chist->Fill(t.tof_chist);  // Track fit chi-square
				histset.h_tof_chisc->Fill(t.tof_chisc);  // Charge fit chi-square
			}
			// RICH quality histograms
			if(t.richBeta > 0){
				histset.h_rich_pb->Fill(t.rich_pb);      // Photoelectron background
				histset.h_rich_hit->Fill(t.rich_hit);    // Number of hits
				histset.h_rich_usedm->Fill(t.rich_usedm);// Number of used PMTs in ring reconstruction
				histset.h_rich_pmt->Fill(t.rich_pmt);    // Total number of fired PMTs

				// RICH photoelectron histograms
				histset.h_rich_npe0->Fill(t.rich_npe[0]); // Direct photoelectrons
				histset.h_rich_npe2->Fill(t.rich_npe[2]); // Total photoelectrons
				if(t.rich_npe[2] > 0) {
					// Ratio of direct to total photoelectrons
					histset.h_rich_npe_ratio->Fill(t.rich_npe[0]/t.rich_npe[2]);
				}
			}
		}

		if(!trackerCut || !cutoffCut || (t.richBeta <= 0)) continue;

		double invRichBeta = 1. / t.richBeta;
		int det_idx = t.rich_NaF ? 0 : 1;
		double rig_cut = t.rich_NaF ? 50 : 100;

		if(t.InnerRig > rig_cut) {
			// Fill for different RICH cut levels
			for(int cut_level = 0; cut_level < 4; cut_level++) {
				bool pass_cut = (t.richBeta > 0);
				for(int bit = 0; bit < cut_level; bit++) {
					if(!(t.cutStatus & (1 << (9 + bit)))) {
						pass_cut = false;
						break;
					}
				}
				if(cut_level==3 && (pass_cut != richCut))
					continue;
				if(pass_cut) {
					histset.h_rich_beta_theta[det_idx][cut_level]->Fill(t.rich_theta, invRichBeta);
					histset.h_rich_beta_phi[det_idx][cut_level]->Fill(t.rich_phi, invRichBeta);
				}
			}
		}
	}

	// Save histograms to ROOT file
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 4; j++) {
			histset.h_rich_beta_theta[i][j]->Write();
			histset.h_rich_beta_phi[i][j]->Write();
		}
	}

	histset.h_tof_chist->Write();
	histset.h_tof_chisc->Write();
	histset.h_rich_pb->Write();
	histset.h_rich_hit->Write();
	histset.h_rich_usedm->Write();
	histset.h_rich_npe0->Write();
	histset.h_rich_npe2->Write();
	histset.h_rich_npe_ratio->Write();
	histset.h_rich_pmt->Write();

	outfile->Close();
	delete outfile;

	// Save histograms to PDF
	SaveAllHistograms(histset, sampleName);

	// Cleanup
	for(auto& vec1 : histset.h_rich_beta_theta) {
		for(auto& hist : vec1) delete hist;
	}
	for(auto& vec1 : histset.h_rich_beta_phi) {
		for(auto& hist : vec1) delete hist;
	}
	delete histset.h_tof_chist;
	delete histset.h_tof_chisc;
	delete histset.h_rich_pb;
	delete histset.h_rich_hit;
	delete histset.h_rich_usedm;
	delete histset.h_rich_npe0;
	delete histset.h_rich_npe2;
	delete histset.h_rich_npe_ratio;
	delete histset.h_rich_pmt;

	f->Close();
	delete f;
}

void AnalyzeIsotope(const std::string& isotope) {
	std::vector<std::string> filenames;
	std::string baseDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/";

	if(isotope == "Li") {
		filenames = {
			baseDir + "Li6.root",
			baseDir + "Li7.root",
			baseDir + "Lit.root"
		};
	}
	else if(isotope == "Be") {
		filenames = {
			baseDir + "Be7.root",
			baseDir + "Be9.root",
			baseDir + "Be10.root",
			baseDir + "Ber.root"
		};
	}
	else if(isotope == "B") {
		filenames = {
			baseDir + "B10.root",
			baseDir + "B11.root",
			baseDir + "Bor.root"
		};
	}
	else {
		std::cout << "Unknown isotope: " << isotope << std::endl;
		return;
	}

	for(const auto& filename : filenames) {
		ProcessFile(filename);
	}
}

void analysTree() {
	AnalyzeIsotope("Li");
}
