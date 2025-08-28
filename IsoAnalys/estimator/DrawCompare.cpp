#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>

const double Be_bins[] = {0.08,0.13,0.17,0.21,0.27,0.33,0.41,0.49,0.59,0.70,0.82,0.96,1.11,1.28,1.47,1.68,1.91,2.16,2.44,2.73,3.06,3.41,3.79,4.20,4.65,5.14,5.64,6.18,6.78,7.42,8.12,8.86,9.66,10.51,11.45,12.45,13.50,14.65};
const int nBins = sizeof(Be_bins)/sizeof(Be_bins[0]);

void FillNewHist(TH1F* src, TH1F* dst, double xMin, double xMax, bool useErrors = true, double maxlimit = 0.99) {
    if (!src || !dst) return;
    for (int i = 1; i <= src->GetNbinsX(); ++i) {
        double binCenter = src->GetBinCenter(i);
        if (binCenter >= xMin && binCenter <= xMax) {
            int dstBin = dst->FindBin(binCenter);
            if (dstBin > 0 && dstBin <= dst->GetNbinsX() && src->GetBinContent(i) > 0 && src->GetBinContent(i) < maxlimit) {
                dst->SetBinContent(dstBin, src->GetBinContent(i));
                if (useErrors) dst->SetBinError(dstBin, src->GetBinError(i));
            }
        }
    }
}

TH1F* CalculateBe10(TH1F* h7, TH1F* h9) {
    TH1F* h10 = (TH1F*)h7->Clone("h10");
    for (int i = 1; i <= h7->GetNbinsX(); ++i) {
        h10->SetBinContent(i, 1.0 - h7->GetBinContent(i) - h9->GetBinContent(i));
        h10->SetBinError(i, sqrt(pow(h7->GetBinError(i),2) + pow(h9->GetBinError(i),2)));
    }
    return h10;
}

std::pair<double,double> GetYRange(const std::vector<TH1F*>& hists) {
    double minVal = 1e9, maxVal = -1e9;
    for (auto h : hists) {
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double val = h->GetBinContent(i);
            if (val > 0 && val < 100) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
            }
        }
    }
    return {minVal*0.5, maxVal*1.5};
}

void DrawCompare() {
    gStyle->SetOptStat(0);
    struct Detector {const char* name; double minE; double maxE; int color;} 
    detectors[] = {{"ToF",0.4,1.5,kBlue}, {"NaF",1.0,4.0,kGreen+2}, {"AGL",2.2,15,kOrange+1}};
    const int nDet = sizeof(detectors)/sizeof(detectors[0]);
    const string methodname[] = {"Yan 1/Mass", "Lu 1/Mass", "Lu estimator"};
    
    struct Method {const char* name; const char* file; bool isTH2D; int marker; int color;} 
    methods[] = {{"Yan Mass","/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/result_%sMass_Be_compareLu.root",false,20,kBlack},
                 {"Lu Mass","/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/RMassTempFit.root",true,21,kRed},
                 {"Lu Est","/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/EstTempFit.root",true,22,kBlue}};
    const int nMethod = sizeof(methods)/sizeof(methods[0]);

    const char* vars[] = {"Be7", "Be9", "Be10", "chi2", "sgf"};
    const int nVar = sizeof(vars)/sizeof(vars[0]);
    std::map<std::string, std::vector<std::vector<TH1F*>>> allHists;
    
    // First collect all histograms
    for (int iVar = 0; iVar < nVar; ++iVar) {
        const char* var = vars[iVar];
        std::vector<std::vector<TH1F*>> methodHists;
        for (int iMethod = 0; iMethod < nMethod; ++iMethod) {
            std::vector<TH1F*> detectorHists;
            for (int iDet = 0; iDet < nDet; ++iDet) {
                TH1F* h = new TH1F(Form("h_%s_%s_%s",methods[iMethod].name,detectors[iDet].name,var), 
                                  var, nBins-1, Be_bins);
                if (methods[iMethod].isTH2D) {
                    TFile* f = TFile::Open(methods[iMethod].file);
                    if (strcmp(var,"Be10") == 0) {
                        TH2D* h2d7 = (TH2D*)f->Get("hBe7Frac");
                        TH2D* h2d9 = (TH2D*)f->Get("hBe9Frac");
                        int detIdx = strcmp(detectors[iDet].name,"ToF")==0 ? 1 : 
                                   strcmp(detectors[iDet].name,"NaF")==0 ? 2 : 3;
                        TH1F* h7 = (TH1F*)h2d7->ProjectionX("temp7",detIdx,detIdx);
                        TH1F* h9 = (TH1F*)h2d9->ProjectionX("temp9",detIdx,detIdx);
                        TH1F* h10 = CalculateBe10(h7,h9);
                        FillNewHist(h10, h, detectors[iDet].minE, detectors[iDet].maxE, true, 0.9);
                        delete h7; delete h9; delete h10;
                    } else {
                        TH2D* h2d = (TH2D*)f->Get(strcmp(var,"chi2")==0 ? "hChi2NDF" : 
                                                 strcmp(var,"sgf")==0 ? "hBeSgf" : Form("h%sFrac",var));
                        int detIdx = strcmp(detectors[iDet].name,"ToF")==0 ? 1 : 
                                   strcmp(detectors[iDet].name,"NaF")==0 ? 2 : 3;
                        TH1F* h1d = (TH1F*)h2d->ProjectionX("temp",detIdx,detIdx);
                        bool NoErr = (strcmp(var,"chi2")==0 || strcmp(var,"sgf")==0 );
                        double maxlimit = NoErr ? 100 : 0.95;
                        FillNewHist(h1d, h, detectors[iDet].minE, detectors[iDet].maxE, !NoErr, maxlimit);
                        delete h1d;
                    }
                    f->Close();
                } else {
                    TFile* f = TFile::Open(Form(methods[iMethod].file,detectors[iDet].name));
                    if (strcmp(var,"Be10") == 0) {
                        TH1F* h7 = (TH1F*)f->Get("h_best_Be7_frac");
                        TH1F* h9 = (TH1F*)f->Get("h_best_Be9_frac");
                        TH1F* h10 = CalculateBe10(h7,h9);
                        FillNewHist(h10, h, detectors[iDet].minE, detectors[iDet].maxE, true, 0.9);
                        delete h10;
                    } else {
                        TH1F* h1d = (TH1F*)f->Get(strcmp(var,"chi2")==0 ? "h_best_chi2" :
                                                 strcmp(var,"sgf")==0 ? "h_best_sgf" : Form("h_best_%s_frac",var));
                        bool NoErr = (strcmp(var,"chi2")==0 || strcmp(var,"sgf")==0 );
                        double maxlimit = NoErr ? 100 : 0.95;
                        FillNewHist(h1d, h, detectors[iDet].minE, detectors[iDet].maxE, !NoErr, maxlimit);
                    }
                    f->Close();
                }
                detectorHists.push_back(h);
            }
            methodHists.push_back(detectorHists);
        }
        allHists[vars[iVar]] = methodHists;
    }

    // Draw detector-specific plots
    for (const char* var : vars) {
        for (int iDet = 0; iDet < nDet; ++iDet) {
            TCanvas* c = new TCanvas(); c->SetLogx();
            TLegend* leg = new TLegend(0.15,0.72,0.3,0.88);
            leg->SetFillStyle(0); leg->SetBorderSize(1); //eg->SetTextSize(0.035);
            
            TH1F* hTemp = new TH1F("hTemp","",nBins-1, Be_bins);
            hTemp->GetXaxis()->SetRangeUser(detectors[iDet].minE, detectors[iDet].maxE);
            hTemp->GetXaxis()->SetTitle("Ek [GeV/n]");
            if (strcmp(var,"sgf") == 0) {
                hTemp->SetTitle(Form("%s significance",detectors[iDet].name));
                hTemp->GetYaxis()->SetTitle("Significance");
            } else if (strcmp(var,"chi2") == 0){
                hTemp->GetYaxis()->SetTitle("Chi2/NDF");
            }
            else {
                hTemp->GetYaxis()->SetTitle(Form("%s %s", var, "Fraction"));
            }
            hTemp->GetXaxis()->SetNdivisions(520);
            
            std::vector<TH1F*> plotHists;
            for (const auto& methodHists : allHists[var]) plotHists.push_back(methodHists[iDet]);
            auto yRange = GetYRange(plotHists);
            hTemp->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
            hTemp->Draw();

            for (int iMethod = 0; iMethod < nMethod; ++iMethod) {
                TH1F* h = allHists[var][iMethod][iDet];
                h->SetMarkerStyle(methods[iMethod].marker);
                h->SetMarkerColor(methods[iMethod].color);
                h->SetLineColor(methods[iMethod].color);
                h->SetMarkerSize(1.2);
                h->GetXaxis()->SetNdivisions(520);
                
                if (strcmp(var,"chi2")==0 || strcmp(var,"sgf")==0) h->Draw("P SAME");
                else h->Draw("PE SAME");
                leg->AddEntry(h, methods[iMethod].name, "P");
            }
            leg->Draw();
            c->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/%s_%s.png",
                          var,detectors[iDet].name));
            delete c;
            delete hTemp;
        }

        // Draw method-specific plots
        for (int iMethod = 0; iMethod < nMethod; ++iMethod) {
            TCanvas* c = new TCanvas(); c->SetLogx();
            TLegend* leg = new TLegend(0.15,0.72,0.3,0.88);
            leg->SetFillStyle(0); leg->SetBorderSize(1); leg->SetTextSize(0.035);
            
            TH1F* hTemp = new TH1F("hTemp","",nBins-1, Be_bins);
            hTemp->GetXaxis()->SetRangeUser(0.3, 15.0);
            hTemp->GetXaxis()->SetTitle("Ek [GeV/n]");
            if (strcmp(var,"sgf") == 0) {
                hTemp->SetTitle(Form("%s significance", methodname[iMethod].c_str()));
                hTemp->GetYaxis()->SetTitle("Significance");
            } else if (strcmp(var,"chi2") == 0){
                hTemp->GetYaxis()->SetTitle("Chi2/NDF");
            }
            else {
                hTemp->GetYaxis()->SetTitle(Form("%s %s", var, "Fraction"));
            }
            hTemp->GetXaxis()->SetNdivisions(520);
            
            auto yRange = GetYRange({allHists[var][iMethod][0], 
                                   allHists[var][iMethod][1], 
                                   allHists[var][iMethod][2]});
            hTemp->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
            hTemp->Draw();

            for (int iDet = 0; iDet < nDet; ++iDet) {
                TH1F* h = allHists[var][iMethod][iDet];
                h->SetMarkerStyle(20);
                h->SetMarkerColor(detectors[iDet].color);
                h->SetLineColor(detectors[iDet].color);
                h->SetMarkerSize(1.2);
                
                if (strcmp(var,"chi2")==0 || strcmp(var,"sgf")==0) h->Draw("P SAME");
                else h->Draw("PE SAME");
                leg->AddEntry(h, detectors[iDet].name, "P");
            }
            leg->Draw();
            c->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/%s_%s.png",
                          var,methods[iMethod].name));
            delete c;
            delete hTemp;
        }
    }

    // Cleanup
    for (auto& var : vars) {
        for (auto& methodHists : allHists[var]) {
            for (auto& h : methodHists) delete h;
        }
    }
}