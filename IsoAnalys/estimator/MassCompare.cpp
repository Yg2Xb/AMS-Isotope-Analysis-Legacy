#include <TFile.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPDF.h>
#include <TPaveText.h>
#include <TF1.h>
#include <iostream>
#include <vector>
#include <utility>

// 两次高斯拟合函数
pair<double, double> doubleFitGaussian(TH1* hist, double mass) {
    double center = 1.0/(mass*0.931);
    TF1* f1 = new TF1("f1", "gaus", center-0.01, center+0.01);
    hist->Fit(f1, "QRN");
    double mean1 = f1->GetParameter(1);
    double sigma1 = f1->GetParameter(2);
    TF1* f2 = new TF1("f2", "gaus", mean1-2.*sigma1, mean1+2.*sigma1);
    hist->Fit(f2, "R");
    cout<<"range: "<<mean1<<"+-"<<2*sigma1<<endl;
    cout << "Fit Chi2/NDF = " << f2->GetChisquare() << "/" << f2->GetNDF() 
     << " = " << f2->GetChisquare()/f2->GetNDF() << endl;
    double mean = f2->GetParameter(1);
    double sigma = f2->GetParameter(2);
    f2->Draw("same");
    delete f1; delete f2;
    return {mean, sigma};
}

void setHistStyle(TH1* hist, const char* title, const char* xtitle, const char* ytitle, 
                 int color = 1, int marker = 20) {
    hist->SetMarkerStyle(marker);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
}

void MassCompare() {
    gStyle->SetOptStat(0);
    
    const vector<double> Be_bins = {
        0.08,0.13,0.17,0.21,0.27,0.33,0.41,0.49,0.59,0.70,0.82,0.96,1.11,1.28,1.47,1.68,1.91,
        2.16,2.44,2.73,3.06,3.41,3.79,4.20,4.65,5.14,5.64,6.18,6.78,7.42,8.12,8.86,9.66,10.51,
        11.45,12.45,13.50,14.65
    };

    struct DetInfo {
        const char* name;
        const char* yan_hist;
        const char* lu_hist;
        double minE;
        double maxE;
    } dets[] = {
        {"ToF", "h_inv_ToFMass_alpha_20", "hist011", 0.4, 1.5},
        {"NaF", "h_inv_NaFMass_alpha_20", "hist021", 1.0, 4.0},
        {"AGL", "h_inv_AGLMass_alpha_20", "hist031", 2.2, 15.0}
    };

    const double masses[] = {7.0, 9.0, 10.0};
    const int colors[] = {kRed, kBlue, kGreen+2, kOrange+1, kMagenta, kBlack};
    const char* isotopes[] = {"Be7", "Be9", "Be10"};

    // 创建直方图
    vector<TH1D*> h_yan_mean, h_lu_mean, h_yan_sigma, h_lu_sigma;
    vector<TH1D*> h_yan_counts(3), h_lu_counts(3);  // 每个探测器一个counts直方图

    // 初始化直方图
    for(int i = 0; i < 3; i++) {
        h_yan_mean.push_back(new TH1D(Form("h_yan_mean_Be%d", int(masses[i])), 
            Form("Mean vs Ek/n; Ek/n(GeV/n);1/Mass Mean(1/GeV)"), Be_bins.size()-1, &Be_bins[0]));
        h_lu_mean.push_back(new TH1D(Form("h_lu_mean_Be%d", int(masses[i])), 
            Form("Mean vs Ek/n; Ek/n(GeV/n);1/Mass Mean(1/GeV)"), Be_bins.size()-1, &Be_bins[0]));
        h_yan_sigma.push_back(new TH1D(Form("h_yan_sigma_Be%d", int(masses[i])), 
            Form("Sigma vs Ek/n;Ek/n(GeV/n);1/Mass Sigma (1/GeV)"), Be_bins.size()-1, &Be_bins[0]));
        h_lu_sigma.push_back(new TH1D(Form("h_lu_sigma_Be%d", int(masses[i])), 
            Form("Sigma vs Ek/n;Ek/n(GeV/n);1/Mass Sigma (1/GeV)"), Be_bins.size()-1, &Be_bins[0]));
            
        setHistStyle(h_yan_mean[i], "", "", "", colors[i], 20);
        setHistStyle(h_lu_mean[i], "", "", "", colors[i+3], 21);
        setHistStyle(h_yan_sigma[i], "", "", "", colors[i], 20);
        setHistStyle(h_lu_sigma[i], "", "", "", colors[i+3], 21);
    }

    // 为每个探测器创建counts直方图
    for(int i = 0; i < 3; i++) {
        h_yan_counts[i] = new TH1D(Form("h_yan_counts_%s", dets[i].name), 
            Form("%s ISS Counts;Ek/n (GeV/n);Counts", dets[i].name), Be_bins.size()-1, &Be_bins[0]);
        h_lu_counts[i] = new TH1D(Form("h_lu_counts_%s", dets[i].name), 
            Form("%s ISS Counts;Ek/n (GeV/n);Counts", dets[i].name), Be_bins.size()-1, &Be_bins[0]);
        setHistStyle(h_yan_counts[i], "", "", "", kBlack, 20);
        setHistStyle(h_lu_counts[i], "", "", "", kRed, 21);
    }
// Open files
    TFile *f_yan_7 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7_temp.root");
    TFile *f_yan_9 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be9_temp.root");
    TFile *f_yan_10 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be10_temp.root");
    TFile *f_yan_iss = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber_temp.root");
    TFile *f_lu_7 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be7.root");
    TFile *f_lu_9 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be9.root");
    TFile *f_lu_10 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be10.root");
    TFile *f_lu_iss = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_iss.root");

    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/MassCompare.pdf[");

    // Loop over detectors for PDF generation and data collection
    for (int idet = 0; idet < 3; ++idet) {
        const auto& det = dets[idet];
        
        // Get histograms
        TH2F* h_yan_7 = (TH2F*)f_yan_7->Get(det.yan_hist);
        TH2F* h_yan_9 = (TH2F*)f_yan_9->Get(det.yan_hist);
        TH2F* h_yan_10 = (TH2F*)f_yan_10->Get(det.yan_hist);
        TH2F* h_yan_iss = (TH2F*)f_yan_iss->Get(det.yan_hist);
        TH2D* h_lu_7 = (TH2D*)f_lu_7->Get(det.lu_hist);
        TH2D* h_lu_9 = (TH2D*)f_lu_9->Get(det.lu_hist);
        TH2D* h_lu_10 = (TH2D*)f_lu_10->Get(det.lu_hist);
        TH2D* h_lu_iss = (TH2D*)f_lu_iss->Get(det.lu_hist);

        // Loop over Ek/n bins
        for (int binE = 1; binE <= h_yan_7->GetNbinsX(); ++binE) {
            double Ek = h_yan_7->GetYaxis()->GetBinCenter(binE);
            if (Ek < det.minE || Ek > det.maxE) continue;
            double Ek2 = h_yan_7->GetYaxis()->GetBinCenter(binE+1);

            // Template canvas
            TCanvas* c1 = new TCanvas("c1", "", 800, 600);
            TLegend* leg1 = new TLegend(0.65, 0.65, 0.85, 0.85);
            leg1->SetFillStyle(0);
            leg1->SetBorderSize(1);

            // Process templates
            TH1D* temp_hists[6];
            temp_hists[0] = h_yan_7->ProjectionX(Form("yan_7_%s_%.2f", det.name, Ek), binE, binE);
            temp_hists[1] = h_yan_9->ProjectionX(Form("yan_9_%s_%.2f", det.name, Ek), binE, binE);
            temp_hists[2] = h_yan_10->ProjectionX(Form("yan_10_%s_%.2f", det.name, Ek), binE, binE);
            
            int luBinE = h_lu_7->GetXaxis()->FindBin(Ek);
            temp_hists[3] = h_lu_7->ProjectionY(Form("lu_7_%s_%.2f", det.name, Ek), luBinE, luBinE);
            temp_hists[4] = h_lu_9->ProjectionY(Form("lu_9_%s_%.2f", det.name, Ek), luBinE, luBinE);
            temp_hists[5] = h_lu_10->ProjectionY(Form("lu_10_%s_%.2f", det.name, Ek), luBinE, luBinE);

            // Normalize and process templates
            double maxY = 0;
            for (int i = 0; i < 6; ++i) {
                double integral = temp_hists[i]->Integral(temp_hists[i]->FindBin(0.06), temp_hists[i]->FindBin(0.22));
                if (integral > 0) temp_hists[i]->Scale(1.0/integral);
                maxY = std::max(maxY, temp_hists[i]->GetMaximum());
                
                if (i < 3) {
                    auto [mean_yan, sigma_yan] = doubleFitGaussian(temp_hists[i], masses[i]);
                    auto [mean_lu, sigma_lu] = doubleFitGaussian(temp_hists[i+3], masses[i]);
                    
                    h_yan_mean[i]->SetBinContent(binE, mean_yan);
                    //h_yan_mean[i]->SetBinError(binE, sigma_yan);
                    h_lu_mean[i]->SetBinContent(binE, mean_lu);
                    //h_lu_mean[i]->SetBinError(binE, sigma_lu);
                    
                    h_yan_sigma[i]->SetBinContent(binE, sigma_yan);
                    h_lu_sigma[i]->SetBinContent(binE, sigma_lu);
                }
            }

            // Draw templates
            for (int i = 0; i < 6; ++i) {
                temp_hists[i]->SetLineColor(colors[i]);
                temp_hists[i]->SetLineWidth(2);
                if (i == 0) {
                    temp_hists[i]->SetTitle(Form("%s Templates Ek=%.2f-%.2f GeV/n", det.name, Ek, Ek2));
                    temp_hists[i]->GetXaxis()->SetTitle("1/Mass [1/GeV]");
                    temp_hists[i]->GetYaxis()->SetTitle("Normalized Entries");
                    temp_hists[i]->GetXaxis()->SetRangeUser(0.06, 0.22);
                    temp_hists[i]->GetYaxis()->SetRangeUser(0, maxY*1.2);
                    temp_hists[i]->Draw("HIST");
                } else {
                    temp_hists[i]->Draw("HIST SAME");
                }
                leg1->AddEntry(temp_hists[i], Form("%s %s", i<3?"Yan":"Lu", isotopes[i%3]), "L");
            }
            leg1->Draw();
            c1->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/MassCompare.pdf");
            delete c1;

            // Process ISS data
            TH1D* proj_yan_iss = h_yan_iss->ProjectionX(Form("yan_iss_%s_%.2f", det.name, Ek), binE, binE);
            TH1D* proj_lu_iss = h_lu_iss->ProjectionY(Form("lu_iss_%s_%.2f", det.name, Ek), luBinE, luBinE);

            double yan_integral = proj_yan_iss->Integral(proj_yan_iss->FindBin(0.06), proj_yan_iss->FindBin(0.22));
            double lu_integral = proj_lu_iss->Integral(proj_lu_iss->FindBin(0.06), proj_lu_iss->FindBin(0.22));

            h_yan_counts[idet]->SetBinContent(binE, yan_integral);
            h_yan_counts[idet]->SetBinError(binE, sqrt(yan_integral));
            h_lu_counts[idet]->SetBinContent(binE, lu_integral);
            h_lu_counts[idet]->SetBinError(binE, sqrt(lu_integral));

            if (yan_integral > 0 && lu_integral > 0) {
                TCanvas* c2 = new TCanvas("c2", "", 800, 600);
                TLegend* leg2 = new TLegend(0.65, 0.65, 0.85, 0.85);
                leg2->SetFillStyle(0);
                leg2->SetBorderSize(1);

                proj_yan_iss->Scale(1.0/yan_integral);
                proj_lu_iss->Scale(1.0/lu_integral);

                proj_yan_iss->SetLineColor(kBlack);
                proj_lu_iss->SetLineColor(kRed);
                proj_yan_iss->SetLineWidth(2);
                proj_lu_iss->SetLineWidth(2);

                proj_yan_iss->SetTitle(Form("%s ISS Ek=%.2f-%.2f GeV/n", det.name, Ek, Ek2));
                proj_yan_iss->GetXaxis()->SetTitle("1/Mass [1/GeV]");
                proj_yan_iss->GetYaxis()->SetTitle("Normalized Entries");
                proj_yan_iss->GetXaxis()->SetRangeUser(0.06, 0.22);
                double maxY_iss = std::max(proj_yan_iss->GetMaximum(), proj_lu_iss->GetMaximum());
                proj_yan_iss->GetYaxis()->SetRangeUser(0, maxY_iss*1.2);

                proj_yan_iss->Draw("HIST");
                proj_lu_iss->Draw("HIST SAME");

                leg2->AddEntry(proj_yan_iss, "Yan ISS Data", "L");
                leg2->AddEntry(proj_lu_iss, "Lu ISS Data", "L");
                leg2->Draw();

                TPaveText* pt = new TPaveText(0.2, 0.75, 0.35, 0.88, "NDC");
                pt->SetBorderSize(0);
                pt->SetFillColor(0);
                pt->SetFillStyle(0);
                pt->SetTextSize(0.03);
                pt->AddText(Form("Yan Counts: %.0f", yan_integral));
                pt->AddText(Form("Lu Counts: %.0f", lu_integral));
                pt->Draw();

                c2->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/MassCompare.pdf");
                delete c2;
            }

            // Cleanup
            for (auto h : temp_hists) delete h;
            delete proj_yan_iss;
            delete proj_lu_iss;
        }
        // Generate final plots for each detector
        TCanvas* c_final = new TCanvas("c_final", "", 800, 600);
        c_final->SetLogx();
        // Mean plot
        c_final->Clear();
        TLegend* leg_mean = new TLegend(0.65, 0.7, 0.85, 0.95);
        leg_mean->SetFillStyle(0);
        leg_mean->SetBorderSize(1);
        
        bool first = true;
        for(int i = 0; i < 3; i++) {
            if(first) {
                h_yan_mean[i]->SetTitle(Form("%s Mean vs Ek/n", det.name));
                h_yan_mean[i]->GetYaxis()->SetRangeUser(0.07, 0.19);
                h_yan_mean[i]->GetXaxis()->SetRangeUser(det.minE, det.maxE);
                h_yan_mean[i]->Draw("P");
                first = false;
            } else {
                h_yan_mean[i]->Draw("P SAME");
            }
            h_lu_mean[i]->Draw("P SAME");
            
            leg_mean->AddEntry(h_yan_mean[i], Form("Yan Be%d", int(masses[i])), "P");
            leg_mean->AddEntry(h_lu_mean[i], Form("Lu Be%d", int(masses[i])), "P");
        }
        leg_mean->Draw();
        c_final->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/%s_Mean.png", det.name));
        delete leg_mean;

        // Sigma plot
        c_final->Clear();
        TLegend* leg_sigma = new TLegend(0.2, 0.65, 0.5, 0.87);
        leg_sigma->SetFillStyle(0);
        leg_sigma->SetBorderSize(1);
        
        first = true;
        for(int i = 0; i < 3; i++) {
            if(first) {
                h_yan_sigma[i]->SetTitle(Form("%s Sigma vs Ek/n", det.name));
                h_yan_sigma[i]->GetYaxis()->SetRangeUser(0.009, 0.028);
                h_yan_sigma[i]->GetXaxis()->SetRangeUser(det.minE, det.maxE);
                h_yan_sigma[i]->Draw("P");
                first = false;
            } else {
                h_yan_sigma[i]->Draw("P SAME");
            }
            h_lu_sigma[i]->Draw("P SAME");
            
            leg_sigma->AddEntry(h_yan_sigma[i], Form("Yan Be%d", int(masses[i])), "P");
            leg_sigma->AddEntry(h_lu_sigma[i], Form("Lu Be%d", int(masses[i])), "P");
        }
        leg_sigma->Draw();
        c_final->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/%s_Sigma.png", det.name));
        delete leg_sigma;

        // Counts plot
        c_final->Clear();
        TLegend* leg_counts = new TLegend(0.74, 0.74, 0.88, 0.88);
        leg_counts->SetFillStyle(0);
        leg_counts->SetBorderSize(1);
        double ymax = h_lu_counts[idet]->GetMaximum(); 
        double ymin = h_yan_counts[idet]->GetMinimum(); 
        h_lu_counts[idet]->GetXaxis()->SetRangeUser(det.minE, det.maxE); 
        h_lu_counts[idet]->GetYaxis()->SetRangeUser(0.8*ymin, 1.2*ymax); 
        h_lu_counts[idet]->Draw("P");
        h_yan_counts[idet]->Draw("P SAME");
        leg_counts->AddEntry(h_yan_counts[idet], "Yan ISS", "P");
        leg_counts->AddEntry(h_lu_counts[idet], "Lu ISS", "P");
        leg_counts->Draw();
        
        c_final->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/%s_Counts.png", det.name));
        delete leg_counts;
        delete c_final;
    }
    c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/compare/MassCompare.pdf]");
    

    // Cleanup
    for(auto h : h_yan_mean) delete h;
    for(auto h : h_lu_mean) delete h;
    for(auto h : h_yan_sigma) delete h;
    for(auto h : h_lu_sigma) delete h;
    for(auto h : h_yan_counts) delete h;
    for(auto h : h_lu_counts) delete h;

    // Close files
    f_yan_7->Close();
    f_yan_9->Close();
    f_yan_10->Close();
    f_yan_iss->Close();
    f_lu_7->Close();
    f_lu_9->Close();
    f_lu_10->Close();
    f_lu_iss->Close();
}