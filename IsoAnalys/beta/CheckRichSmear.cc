#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPad.h>
#include <TString.h> // For Form
#include <TLine.h>   // For pull plot line at 1
#include <iostream>
#include <vector>
#include <string>
#include <map> // For mapping index to mass

#include "/afs/cern.ch/work/z/zuhao/public/yanzx/IsoAnalys/TempFit/richTune.cc"
#include "../Tool.h"
using namespace AMS_Iso; // Assuming isValidBeta is here

// --- Configuration ---
const double charge = 4.0; // Set charge to Beryllium (Z=4)
const char* iss_file_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root";
const std::vector<std::string> mc_file_paths = {
    "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7.root",
    "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be9.root",
    "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be10.root"
};
const std::vector<std::string> mc_labels = {"MC Be7", "MC Be9", "MC Be10"};
const std::vector<int> mc_colors = {kRed, kBlue, kGreen + 2};
const char* output_dir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/richBeta/";
const char* output_pdf_name = "RICH_Smear_Comparison_Ber.pdf"; // Single output PDF
const char* treeName = "saveTree";
// --------------------

void NormalizeHist(TH1F* h) {
    if(h->Integral() > 0) {
        double binWidth = h->GetBinWidth(1);
        h->Scale(1.0/(h->Integral() * binWidth));
    }
}

// --- Function to Process a Single Tree ---
void ProcessTree(TFile* file, const char* treeName,
                 TH1F* h_AGL_Target, TH1F* h_NaF_Target, // Histograms to fill (ISS or MC NoSmear)
                 TH1F* h_AGL_Smear_Target, TH1F* h_NaF_Smear_Target, // MC Smear histograms (nullptr for ISS)
                 int ia, // Pass MC index (-1 or invalid for ISS)
                 bool isMC) {
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file: " << (file ? file->GetName() : "nullptr") << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Error: Unable to retrieve tree '" << treeName << "' from file: " << file->GetName() << std::endl;
        return;
    }
   
    std::vector<std::vector<double>> beta_bins;
    int masses[3] = {7,9,10};
    
    for (int i = 0; i < 3; i++) {
        cout << "Mass: " << masses[i] << endl;
        std::vector<double> Rbins_beta;
        for (int j = 0; j < Binning::RigidityBins.size(); j++) {
            double rig = Binning::RigidityBins[j];
            double beta = rigidityToBeta(rig, 4, masses[i], false);
            Rbins_beta.push_back(beta);
        }
        
        cout << "Finished mass " << masses[i] << " with " << Rbins_beta.size() << " beta values" << endl;
        beta_bins.push_back(Rbins_beta);
    }

    // Setup branches
    double InnerRig = -1, TOFBeta = -1, NaFBeta = -1, AGLBeta = -1, richBeta = -1, cutOffRig = -1;
    bool rich_NaF = false;
    unsigned int cutStatus = 0;
    double generatedRig = -1, MC_weight_27 = 1.0; // Default weight to 1 for ISS

    // Set branch addresses common to both ISS and MC
    tree->SetBranchAddress("InnerRig", &InnerRig);
    tree->SetBranchAddress("TOFBeta", &TOFBeta);
    tree->SetBranchAddress("NaFBeta", &NaFBeta);
    tree->SetBranchAddress("AGLBeta", &AGLBeta);
    tree->SetBranchAddress("richBeta", &richBeta);
    tree->SetBranchAddress("cutOffRig", &cutOffRig);
    tree->SetBranchAddress("rich_NaF", &rich_NaF);
    tree->SetBranchAddress("cutStatus", &cutStatus);

    double current_mc_mass = 0; // Will be set only for MC
    if (isMC) {
        if (tree->GetBranch("generatedRig")) tree->SetBranchAddress("generatedRig", &generatedRig);
        if (tree->GetBranch("MC_weight_27")) tree->SetBranchAddress("MC_weight_27", &MC_weight_27);
    }

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        //if(i>500000) continue;
        if(i%1000000 == 0) cout<<i<<endl;
        tree->GetEntry(i);

        double weight = 1.;//(isMC && tree->GetBranch("MC_weight_27")) ? MC_weight_27 : 1.0;
        if (weight <= 0) continue; // Skip events with non-positive weights if necessary

        // Pre-selection: Check bit 0
        if(!((cutStatus & 0x01) == 0x01)) continue;

        // Define selection criteria
        bool isAGLEvent = ((cutStatus & 0x04) == 0x04) && AGLBeta > 0; // Bits 0, 2, 6
        bool isNaFEvent = ((cutStatus & 0x04) == 0x04) && NaFBeta > 0; // Bits 0, 2, 5

        // Rigidity cuts
        bool passAGLRigCut = InnerRig > 200;
        bool passNaFRigCut = InnerRig > 100;

        if (isMC) {
            if (!h_AGL_Smear_Target || !h_NaF_Smear_Target) {
                 std::cerr << "Error: Smear target histograms not provided for MC processing." << std::endl;
                 return;
            }
            if (passAGLRigCut && isAGLEvent) {
                // No Smear - Apply weight
                h_AGL_Target->Fill(1.0 / AGLBeta, weight);
                // Smear - Apply weight
                // *** Corrected: Pass current_mc_mass ***
                double smearedBeta = GetSmearRichBeta(charge, ia, AGLBeta, false);
                if (smearedBeta > 0) {
                   h_AGL_Smear_Target->Fill(1.0 / smearedBeta, weight);
                }
            }
            if (passNaFRigCut && isNaFEvent) {
                // No Smear - Apply weight
                h_NaF_Target->Fill(1.0 / NaFBeta, weight);
                // Smear - Apply weight
                // *** Corrected: Pass current_mc_mass ***
                double smearedBeta = GetSmearRichBeta(charge, ia,  NaFBeta, true);
                 if (smearedBeta > 0) {
                    h_NaF_Smear_Target->Fill(1.0 / smearedBeta, weight);
                 }
            }
        } 
        else if(richBeta > 0){ // is ISS Data (weight is 1.0)
            int richbin = findBin(beta_bins[0], richBeta);
            if(false && AGLBeta <= 1 && InnerRig > 200 && !rich_NaF && !(beta_bins[0][richbin] > 1.0005 * rigidityToBeta(cutOffRig, 4, 7, false)))
            {
                cout<<"richBeta: "<<AGLBeta<<endl;
                cout<<"betalowedge: "<<beta_bins[0][richbin]<<endl;
                cout<<"cf beta: "<<rigidityToBeta(cutOffRig, 4, 7, false)<<endl;
                cout<<"cf beta*safe: "<<1.0005 * rigidityToBeta(cutOffRig, 4, 7, false)<<endl;
                cout<<"cutoffrig:"<<cutOffRig<<endl;
                cout<<"inner rig: "<<InnerRig<<endl;
                cout<<""<<endl;
                
            }
            if (passAGLRigCut && isAGLEvent) {
                h_AGL_Target->Fill(1.0 / AGLBeta); // Weight is 1 here
            }
            if (passNaFRigCut && isNaFEvent) {
                h_NaF_Target->Fill(1.0 / NaFBeta); // Weight is 1 here
            }
        }
    }
    std::cout << "Finished processing " << nEntries << " entries from " << file->GetName() << std::endl;
}

// --- Function to Draw a SINGLE Pull Plot (ISS vs One MC) ---
void DrawSinglePullPlot(TPad* pad, TH1F* h_ISS, TH1F* h_MC_Smear, int color, const char* radiatorName, bool isNaF, bool is1) {
    pad->cd();
    pad->SetGridy();
    pad->SetBottomMargin(0.35);
    pad->SetTopMargin(0.1); // Adjusted slightly
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.1);

    int nBins = h_ISS->GetNbinsX();
    TGraphErrors* pull = new TGraphErrors(); // Create the graph for this single comparison
    pull->SetName(Form("pull_%s_%s", radiatorName, h_MC_Smear->GetName())); // Unique name
    pull->SetMarkerStyle(20);
    pull->SetMarkerSize(0.8);
    pull->SetMarkerColor(color);
    pull->SetLineColor(color);

    int pointsAdded = 0;
    for (int bin = 1; bin <= nBins; bin++) {
        double x = h_ISS->GetBinCenter(bin);
        double ISS_val = h_ISS->GetBinContent(bin);
        if(ISS_val<=1) continue;
        double MC_val = h_MC_Smear->GetBinContent(bin);
        // Optional: Add error calculation if needed
        double ISS_err = h_ISS->GetBinError(bin);
        double MC_err = h_MC_Smear->GetBinError(bin);

        if (MC_val > 1e-12 && ISS_val > 1e-12) { // Avoid division by zero and meaningless points
            double pull_val = ISS_val / MC_val;
            cout<<x<<" "<<pull_val<<endl;
            double pull_err = pull_val * TMath::Sqrt( TMath::Power(ISS_err/ISS_val, 2) + TMath::Power(MC_err/MC_val, 2) ); // If errors needed
            pull->SetPoint(pointsAdded, x, pull_val);
            pull->SetPointError(pointsAdded, 0, pull_err); // Add errors if calculated
            pointsAdded++;
        }
    }

    // Define common X-axis range for pull plots
    double xl = isNaF ? 0.995 : 0.997;
    double xu = isNaF ? 1.005 : 1.003;

        pull->Draw(is1 ? "AP" : "Psame"); // Draw axis with points
        pull->SetTitle(""); // Clear default title

        // Configure Axes
        pull->GetXaxis()->SetLimits(xl, xu);
        pull->GetYaxis()->SetRangeUser(0.0, 3.0); // Adjust Y range if needed
        pull->GetXaxis()->SetTitle(Form("1/#beta_{%s}", radiatorName));
        pull->GetYaxis()->SetTitle("ISS / MC");

        // Adjust label and title sizes
        pull->GetXaxis()->SetLabelSize(0.12);
        pull->GetYaxis()->SetLabelSize(0.10);
        pull->GetYaxis()->SetNdivisions(505);
        pull->GetXaxis()->SetTitleSize(0.14);
        pull->GetYaxis()->SetTitleSize(0.14);
        pull->GetYaxis()->SetTitleOffset(0.5);
        pull->GetXaxis()->SetTitleOffset(1.0);

        // Add line at y=1
        TLine *line = new TLine(xl, 1.0, xu, 1.0);
        line->SetLineStyle(kDashed);
        line->SetLineColor(kGray + 1);
        line->Draw("SAME");
}


// --- Main Function ---
void CheckRichSmear() {
    gStyle->SetOptStat(0);

    // --- Create Histograms ---
    int nBinsagl = 60;
    int nBinsnaf = 50;
    double agl_min = 1.0 - 0.003, agl_max = 1.0 + 0.003;
    double naf_min = 1.0 - 0.005,  naf_max = 1.0 + 0.005;

    TH1F* h_ISS_AGL = new TH1F("h_ISS_AGL", ";1/#beta_{AGL};Normalized Counts", nBinsagl, agl_min, agl_max);
    TH1F* h_ISS_NaF = new TH1F("h_ISS_NaF", ";1/#beta_{NaF};Normalized Counts", nBinsnaf, naf_min, naf_max);
    // Style ISS histograms
    h_ISS_AGL->SetMarkerStyle(20); h_ISS_AGL->SetMarkerColor(kBlack); h_ISS_AGL->SetLineColor(kBlack);
    h_ISS_NaF->SetMarkerStyle(20); h_ISS_NaF->SetMarkerColor(kBlack); h_ISS_NaF->SetLineColor(kBlack);

    size_t nMC = mc_file_paths.size();
    std::vector<TH1F*> h_MC_AGL_NoSmear(nMC), h_MC_NaF_NoSmear(nMC);
    std::vector<TH1F*> h_MC_AGL_Smear(nMC), h_MC_NaF_Smear(nMC);

    for (size_t i = 0; i < nMC; ++i) {
        const char* label = mc_labels[i].c_str();
        int color = mc_colors[i];
        // No Smear
        h_MC_AGL_NoSmear[i] = new TH1F(Form("h_MC_AGL_NoSmear_%s", label), Form("%s AGL No Smear;1/#beta_{AGL};Counts", label), nBinsagl, agl_min, agl_max);
        h_MC_NaF_NoSmear[i] = new TH1F(Form("h_MC_NaF_NoSmear_%s", label), Form("%s NaF No Smear;1/#beta_{NaF};Counts", label), nBinsnaf, naf_min, naf_max);
        // Smear
        h_MC_AGL_Smear[i] = new TH1F(Form("h_MC_AGL_Smear_%s", label), Form("%s AGL Smear;1/#beta_{AGL};Counts", label), nBinsagl, agl_min, agl_max);
        h_MC_NaF_Smear[i] = new TH1F(Form("h_MC_NaF_Smear_%s", label), Form("%s NaF Smear;1/#beta_{NaF};Counts", label), nBinsnaf, naf_min, naf_max);

        // Style MC histograms
        h_MC_AGL_NoSmear[i]->SetLineColor(kBlue); h_MC_AGL_NoSmear[i]->SetLineWidth(3);
        h_MC_NaF_NoSmear[i]->SetLineColor(kBlue); h_MC_NaF_NoSmear[i]->SetLineWidth(3);
        h_MC_AGL_Smear[i]->SetLineColor(kRed); h_MC_AGL_Smear[i]->SetLineWidth(3); 
        h_MC_NaF_Smear[i]->SetLineColor(kRed); h_MC_NaF_Smear[i]->SetLineWidth(3); 
    }

    // --- Process Files ---
    TFile* file_ISS = TFile::Open(iss_file_path);
    if (file_ISS && !file_ISS->IsZombie()) {
        ProcessTree(file_ISS, treeName, h_ISS_AGL, h_ISS_NaF, nullptr, nullptr, -1, false); // mc_index = -1 for ISS
        file_ISS->Close();
        delete file_ISS;
    } else {
        std::cerr << "Error: Failed to open ISS file: " << iss_file_path << std::endl; return;
    }

    for (size_t i = 0; i < nMC; ++i) {
        TFile* file_MC = TFile::Open(mc_file_paths[i].c_str());
         if (file_MC && !file_MC->IsZombie()) {
            ProcessTree(file_MC, treeName,
                        h_MC_AGL_NoSmear[i], h_MC_NaF_NoSmear[i],
                        h_MC_AGL_Smear[i], h_MC_NaF_Smear[i],
                        i, true); // Pass index i
            file_MC->Close();
            delete file_MC;
         } else {
             std::cerr << "Error: Failed to open MC file: " << mc_file_paths[i] << std::endl; // Continue processing others
         }
    }

    // --- Normalize Histograms ---
    // Normalize ISS (using sum of weights, which is just entries here)
    // 3. 在ProcessTree中统一使用
    NormalizeHist(h_ISS_AGL);
    NormalizeHist(h_ISS_NaF);
    for(size_t i = 0; i < nMC; ++i) {
        NormalizeHist(h_MC_AGL_NoSmear[i]);
        NormalizeHist(h_MC_NaF_NoSmear[i]);
        NormalizeHist(h_MC_AGL_Smear[i]);
        NormalizeHist(h_MC_NaF_Smear[i]);
    }

    // --- Create Plots (Multi-page PDF) ---
    TCanvas* c = new TCanvas("c", "RICH Smear Comparison", 800, 800);
    TString pdfPath = Form("%s%s", output_dir, output_pdf_name);

    // Open the PDF file
    c->Print(pdfPath + "[");

    for (size_t i = 0; i < nMC; ++i) {
        const char* mcLabel = mc_labels[i].c_str();
        int mcColor = mc_colors[i];

        // --- AGL Plot for MC mass i ---
        c->Clear();
        c->Divide(1, 2, 0, 0);

        // Top Pad (Main Comparison)
        TPad* pad1_agl = (TPad*)c->cd(1);
        pad1_agl->SetPad(0, 0.3, 1, 1);
        pad1_agl->SetBottomMargin(0.02); 
        pad1_agl->SetTopMargin(0.1); 
        pad1_agl->SetLeftMargin(0.15); 
        pad1_agl->SetRightMargin(0.1);
        pad1_agl->SetLogy(0);

        // Draw ISS first for axes
        h_ISS_AGL->Draw("PE1");
        pad1_agl->Update(); // Ensure axes are set before getting range
        double ymax_agl = TMath::Max(h_ISS_AGL->GetMaximum(), TMath::Max(h_MC_AGL_NoSmear[i]->GetMaximum(), h_MC_AGL_Smear[i]->GetMaximum())) * 1.05;
        double ymin_agl = TMath::Min(h_ISS_AGL->GetMinimum(1e-7), TMath::Min(h_MC_AGL_NoSmear[i]->GetMinimum(1e-7), h_MC_AGL_Smear[i]->GetMinimum(1e-7))) * 0.7;
        if (ymin_agl <= 0) ymin_agl = 1e-5; // Ensure positive minimum for log scale
        h_ISS_AGL->GetYaxis()->SetRangeUser(ymin_agl, ymax_agl);
        h_ISS_AGL->GetXaxis()->SetLabelSize(0); h_ISS_AGL->GetXaxis()->SetTitle("");
        h_ISS_AGL->GetYaxis()->SetTitleOffset(1.2); h_ISS_AGL->GetYaxis()->SetTitleSize(0.05); h_ISS_AGL->GetYaxis()->SetLabelSize(0.04);
        h_ISS_AGL->SetTitle(Form("AGL Comparison: ISS vs %s", mcLabel)); // Add title to top plot

        // Draw MC histograms

        h_MC_AGL_NoSmear[i]->Draw("HIST SAME"); // Solid line
        h_MC_AGL_Smear[i]->Draw("HIST SAME");   // Dashed line (style set earlier)
        // Redraw ISS points on top
        h_ISS_AGL->Draw("PE1 SAME");

        // Legend for this plot
        TLegend* leg_AGL = new TLegend(0.65, 0.65, 0.88, 0.88); // Adjusted position
        leg_AGL->AddEntry(h_ISS_AGL, "ISS", "EP");
        leg_AGL->AddEntry(h_MC_AGL_NoSmear[i], Form("%s (No Smear)", mcLabel), "l");
        leg_AGL->AddEntry(h_MC_AGL_Smear[i], Form("%s (Smear)", mcLabel), "l");
        leg_AGL->SetBorderSize(0); leg_AGL->SetFillStyle(0);
        leg_AGL->Draw();

        // Bottom Pad (Pull Plot)
        TPad* pad2_agl = (TPad*)c->cd(2);
        pad2_agl->SetPad(0, 0, 1, 0.3);
        // Pass the SMEARED MC histogram for this mass
        DrawSinglePullPlot(pad2_agl, h_ISS_AGL, h_MC_AGL_Smear[i], kRed, "AGL", false, true);
        DrawSinglePullPlot(pad2_agl, h_ISS_AGL, h_MC_AGL_NoSmear[i], kBlue, "AGL", false, false);

        // Print this page to the PDF
        c->Print(pdfPath);

        // --- NaF Plot for MC mass i ---
        c->Clear();
        c->Divide(1, 2, 0, 0);

        // Top Pad (Main Comparison)
        TPad* pad1_naf = (TPad*)c->cd(1);
        pad1_naf->SetPad(0, 0.3, 1, 1);
        pad1_naf->SetBottomMargin(0.02); 
        pad1_naf->SetTopMargin(0.1); 
        pad1_naf->SetLeftMargin(0.15); 
        pad1_naf->SetRightMargin(0.1);

        h_ISS_NaF->Draw("PE1");
        pad1_naf->Update();
        double ymax_naf = TMath::Max(h_ISS_NaF->GetMaximum(), TMath::Max(h_MC_NaF_NoSmear[i]->GetMaximum(), h_MC_NaF_Smear[i]->GetMaximum())) * 1.05;
        double ymin_naf = TMath::Min(h_ISS_NaF->GetMinimum(1e-7), TMath::Min(h_MC_NaF_NoSmear[i]->GetMinimum(1e-7), h_MC_NaF_Smear[i]->GetMinimum(1e-7))) * 0.7;
         if (ymin_naf <= 0) ymin_naf = 1e-5;
        h_ISS_NaF->GetYaxis()->SetRangeUser(ymin_naf, ymax_naf);
        h_ISS_NaF->GetXaxis()->SetLabelSize(0); h_ISS_NaF->GetXaxis()->SetTitle("");
        h_ISS_NaF->GetYaxis()->SetTitleOffset(1.2); h_ISS_NaF->GetYaxis()->SetTitleSize(0.05); h_ISS_NaF->GetYaxis()->SetLabelSize(0.04);
        h_ISS_NaF->SetTitle(Form("NaF Comparison: ISS vs %s", mcLabel)); // Add title

        h_MC_NaF_NoSmear[i]->Draw("HIST SAME");
        h_MC_NaF_Smear[i]->Draw("HIST SAME");
        h_ISS_NaF->Draw("PE1 SAME");

        TLegend* leg_NaF = new TLegend(0.65, 0.65, 0.88, 0.88);
        leg_NaF->AddEntry(h_ISS_NaF, "ISS (Ber7)", "EP");
        leg_NaF->AddEntry(h_MC_NaF_NoSmear[i], Form("%s (No Smear)", mcLabel), "l");
        leg_NaF->AddEntry(h_MC_NaF_Smear[i], Form("%s (Smear)", mcLabel), "l");
        leg_NaF->SetBorderSize(0); leg_NaF->SetFillStyle(0);
        leg_NaF->Draw();

        // Bottom Pad (Pull Plot)
        TPad* pad2_naf = (TPad*)c->cd(2);
        pad2_naf->SetPad(0, 0, 1, 0.3);
        DrawSinglePullPlot(pad2_naf, h_ISS_NaF, h_MC_NaF_Smear[i], kRed, "NaF", true, true);
        DrawSinglePullPlot(pad2_naf, h_ISS_NaF, h_MC_NaF_NoSmear[i], kBlue, "NaF", true, false);

        // Print this page to the PDF
        c->Print(pdfPath);
    }

    // Close the PDF file
    c->Print(pdfPath + "]");

    // --- Cleanup ---
    delete h_ISS_AGL;
    delete h_ISS_NaF;
    for(auto h : h_MC_AGL_NoSmear) delete h;
    for(auto h : h_MC_NaF_NoSmear) delete h;
    for(auto h : h_MC_AGL_Smear) delete h;
    for(auto h : h_MC_NaF_Smear) delete h;
    delete c;

    std::cout << "Processing and plotting complete. Output saved to: " << pdfPath.Data() << std::endl;
}
