#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

void PlotTofComparison(const char* file1Path, const char* file2Path, const char* outputPath, 
                       const char* particle1, const char* particle2, const char* type) {
    
    TFile *f1 = TFile::Open(file1Path), *f2 = TFile::Open(file2Path);
    TH1D* h1[2], *h2[2];
    
    // Get histograms, h1 = iss, h2 = mc
    for(int i = 0; i < 2; i++) {
        h1[i] = (TH1D*)f1->Get(Form("TOF_%s_%s_%s", "NaF", type, particle1))->Clone();
        h2[i] = (TH1D*)f2->Get(Form("TOF_%s_%s_%s", "NaF", type, particle2))->Clone();
    }

    // For Mean plots, get additional sigma values
    TH1D *h1_sig = nullptr, *h2_sig = nullptr;
    if(strcmp(type, "Mean") == 0) {
        h1_sig = (TH1D*)f1->Get(Form("TOF_NaF_Sig_%s", particle1))->Clone();
        h2_sig = (TH1D*)f2->Get(Form("TOF_NaF_Sig_%s", particle2))->Clone();
        for(int i = 1; i <= h1[0]->GetNbinsX(); i++) {
            //h1[0]->SetBinError(i, h1_sig->GetBinContent(i)); 
            //h2[0]->SetBinError(i, h2_sig->GetBinContent(i));
        }
    }
        
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetLogx(); c->SetGrid(); c->SetLeftMargin(.22);

    h1[0]->SetTitle("");
    h1[0]->GetXaxis()->SetRangeUser(0.51,12); 
    h1[0]->GetYaxis()->SetRangeUser(strcmp(type, "Sig") == 0 ? 0.008 : -0.03, 
                                   strcmp(type, "Sig") == 0 ? 0.014 : 0.01);
    h1[0]->GetYaxis()->SetTitle(Form("Mean of 1/Beta_{TOF}-1/Beta_{NaF}", type));
    h1[0]->SetMarkerSize(.8);
    h1[0]->GetYaxis()->SetLabelSize(0.06);
    h1[0]->GetYaxis()->SetTitleSize(0.06);
    h1[0]->GetXaxis()->SetLabelSize(0.06);
    h1[0]->GetXaxis()->SetTitleSize(0.06);
    h1[0]->GetXaxis()->SetLabelOffset(-0.01);
    h1[0]->GetYaxis()->SetTitleOffset(1.6);
    h1[0]->GetYaxis()->SetNdivisions(508);
    h1[0]->GetXaxis()->SetTitle("NaF Ek/n [GeV/n]");

    h2[0]->SetMarkerColor(kBlue);
    h2[0]->SetLineColor(kBlue);
    h2[0]->SetMarkerSize(.8);
    if(strcmp(type, "Mean") == 0){
        h1[0]->SetLineStyle(2);
        h2[0]->SetLineStyle(2);
    }

    h1[0]->Draw("e[]");
    h2[0]->Draw("e[]SAME");
    if(strcmp(type, "Mean") == 0) {
        h1[1]->SetLineStyle(1);
        h2[1]->SetLineStyle(1);
        h1[1]->SetMarkerColorAlpha(kBlue, 0);
        h2[1]->SetMarkerColorAlpha(kBlue, 0);
        //h1[1]->Draw("e");
        //h2[1]->Draw("esame");
    }

    // Draw additional error bars for Mean plots
    TLegend *leg2 = new TLegend(0.55, 0.2, 0.8, 0.35);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(h2[0], "MC Mean", "p");
    leg2->AddEntry(h1[0], "ISS Mean", "p");
    leg2->Draw("same");
    TLegend *leg = new TLegend(0.7, 0.2, 0.82, 0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h2[0], "MC mean", "P");
    leg->AddEntry(h1[0], "ISS mean", "P");
    //leg->Draw("same");




    c->SaveAs(outputPath);
    delete c; delete leg;
    for(int i = 0; i < 2; i++) { delete h1[i]; delete h2[i]; }
    delete h1_sig; delete h2_sig;
    f1->Close(); f2->Close();
}

void CompareBeIsotopes() {
    const char* path = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/";
    const char* masses[] = {"7", "9", "10"};
    for(auto mass : masses) {
        for(const char* type : {"Mean"}) {
            PlotTofComparison(Form("%sFitTOFinEkBin_Ber%s.root", path, mass),
                            Form("%sFitTOFinEkBin_Be%s.root", path, mass),
                            Form("%sBe%sCP_%s.png", path, mass, type), 
                            Form("Ber%s", mass), Form("Be%s", mass), type);
        }
    }
}

void TofBorCP() { CompareBeIsotopes(); }