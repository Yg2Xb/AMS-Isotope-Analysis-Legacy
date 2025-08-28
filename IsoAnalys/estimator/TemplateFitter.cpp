#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <vector>
#include <string>
#include <iostream>
#include "helper_func.cpp"

using namespace RooFit;

class TemplateFitter {
private:
    FitConfig config;
    TFile* issFile;
    TFile* be7File;
    TFile* be9File;
    TFile* be10File;
    TFile* outFile;
    
    std::vector<TH2D*> issHists;
    std::vector<TH2D*> be7Hists;
    std::vector<TH2D*> be9Hists;
    std::vector<TH2D*> be10Hists;
    
    // 用于存储结果的直方图
    TH2D* hBe7Frac;
    TH2D* hBe9Frac;
    TH2D* hBeSgf;
    TH2D* hChi2NDF;

public:
    TemplateFitter(const FitConfig& cfg) : config(cfg) {
        // 打开文件
        //issFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_iss.root");
        issFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber_temp.root");//use yan iss
        be7File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be7.root");
        be9File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be9.root");
        be10File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be10.root");
        std::string output_base = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/";
        outFile = new TFile((output_base + config.output_prefix + ".root").c_str(), "RECREATE");
        
        // 获取直方图
        const char* detname[] = {"ToF", "NaF", "AGL"};
        for (int i = 0; i < 3; ++i) {
            //TH2D* issHist = (TH2D*)issFile->Get(Form("hist0%d%d", i+1, config.hist_suffix));
            TH2D* issHist = (TH2D*)issFile->Get(Form("h_inv_%sMass_alpha_20", detname[i]));
            TH2D* be7Hist = (TH2D*)be7File->Get(Form("hist0%d%d", i+1, config.hist_suffix));
            TH2D* be9Hist = (TH2D*)be9File->Get(Form("hist0%d%d", i+1, config.hist_suffix));
            TH2D* be10Hist = (TH2D*)be10File->Get(Form("hist0%d%d", i+1, config.hist_suffix));
            
            issHists.push_back(issHist);
            be7Hists.push_back(be7Hist);
            be9Hists.push_back(be9Hist);
            be10Hists.push_back(be10Hist);
        }
        
        // 创建结果直方图
        hBe7Frac = new TH2D("hBe7Frac", "Be7 Fraction;E_{k}/n [GeV/n];Detector Type;Fraction", Be_bins.size()-1, &Be_bins[0], 3, 0, 3);
        hBe9Frac = new TH2D("hBe9Frac", "Be9 Fraction;E_{k}/n [GeV/n];Detector Type;Fraction", Be_bins.size()-1, &Be_bins[0], 3, 0, 3);
        hBeSgf = new TH2D("hBeSgf", "Significance;E_{k}/n [GeV/n];Detector Type;Significance", Be_bins.size()-1, &Be_bins[0], 3, 0, 3);
        hChi2NDF = new TH2D("hChi2NDF", "Chi2/NDF;E_{k}/n [GeV/n];Detector Type;Chi2/NDF", Be_bins.size()-1, &Be_bins[0], 3, 0, 3);
    }
    
    ~TemplateFitter() {

    }
    
    void FitBin(TH1D* hData, TH1D* hBe7, TH1D* hBe9, TH1D* hBe10, int EkBin, const char* detector, TPad* pad1, TPad* pad2, TString frametitle, std::string pdf_path) {
       
        TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
        hBe7->SetLineColor(kBlue);
        hBe9->SetLineColor(kGreen+2);
        hBe10->SetLineColor(kOrange+1);
        hBe7->Draw("hist");
        hBe9->Draw("histsame");
        hBe10->Draw("histsame");
        //c2->Print((pdf_path).c_str());

        // 创建变量
        RooRealVar x("x", config.variable_name.c_str(), config.fit_min, config.fit_max);
        
        // 创建数据集
        RooDataHist dataHist("data", "data", x, Import(*hData));
        RooDataHist be7Hist("be7", "be7", x, Import(*hBe7));
        RooDataHist be9Hist("be9", "be9", x, Import(*hBe9));
        RooDataHist be10Hist("be10", "be10", x, Import(*hBe10));
        
        // 创建模板PDF
        RooHistPdf pdfBe7("pdfBe7", "Be7 PDF", x, be7Hist);
        RooHistPdf pdfBe9("pdfBe9", "Be9 PDF", x, be9Hist);
        RooHistPdf pdfBe10("pdfBe10", "Be10 PDF", x, be10Hist);
        
        // 创建系数
        RooRealVar fracBe7("fracBe7", "Be7 Fraction", 0.6, 0., 1.);
        RooRealVar fracBe9("fracBe9", "Be9 Fraction", 0.3, 0., 1.);
        
        // 创建总PDF
        
        RooAddPdf totalPdf("totalPdf", "Total PDF",
                          RooArgList(pdfBe7, pdfBe9, pdfBe10),
                          RooArgList(fracBe7, fracBe9));
        
        /*
        // 执行拟合
        RooAddPdf totalPdf("totalPdf", "Total PDF",
                          RooArgList(pdfBe7, pdfBe9),
                          RooArgList(fracBe7));
        */

        cout<<"Detector:"<<detector<<" BeginRooFit"<<endl;
        RooFitResult* fitResult = totalPdf.fitTo(dataHist, Save());
        cout<<"Detector:"<<detector<<" EndRooFit"<<endl;
        if (fitResult->status() != 0) {
            std::cout << "NO!!!: Fit failed for " << frametitle.Data() << std::endl;
        }
        else{
            std::cout << "Yes!!!: Fit success for " << frametitle.Data() << std::endl;
        }
        
        // 计算显著度
        SignificanceResult result = findOptimalSignificance(x, &pdfBe10, dataHist, totalPdf, fracBe7, fracBe9, config.fit_min, config.fit_max);

        // 获取结果
        double bestEfficiency = result.efficiency;
        double bestUpperBound = result.upperBound;
        double N_sig = result.N_sig;
        double N_tot = result.N_tot;
        double significance = result.significance;
                
        // 绘制结果
        pad1->cd();
        RooPlot* frame = x.frame(Title(""));
        dataHist.plotOn(frame, Name("Data"), XErrorSize(0));
        totalPdf.plotOn(frame, Name("TotalModel"), LineColor(kRed));
        totalPdf.plotOn(frame, Components(pdfBe7), LineColor(kBlue), Name("Be7"));
        totalPdf.plotOn(frame, Components(pdfBe9), LineColor(kGreen+2), Name("Be9"));
        totalPdf.plotOn(frame, Components(pdfBe10), LineColor(kOrange+1), Name("Be10"));
        
        frame->SetYTitle("Events");
        frame->GetYaxis()->SetTitleOffset(1.2);
        frame->SetTitle(frametitle);
        frame->Draw();

        double chi2 = frame->chiSquare();
        int nBins = frame->GetNbinsX();
        int nParams = totalPdf.getParameters(dataHist)->getSize();
        int ndf = nBins - nParams;
        double chi2_ndf = chi2/ndf;
        double be7frac = fracBe7.getVal(), be7fracErr = fracBe7.getError();
        double be9frac = fracBe9.getVal(), be9fracErr = fracBe9.getError();
        if(fitResult->status() != 0)
        {
            be7frac = 0;
            be7fracErr = 0;
            be9frac = 0;
            be9fracErr = 0;
            N_sig = 0;
            N_tot = 0;
            chi2_ndf = 0;
            significance = 0;
        }
        
        // 添加图例
        TLegend* legend = createLegend(frame);
        legend->Draw();
        
        // 添加拟合信息
        TPaveText* pt = createFitInfo(config.l_xmin, config.l_xmax, 
                        be7frac, be7fracErr,
                        be9frac, be9fracErr,
                        chi2_ndf, fitResult->status(),
                        config.fit_min, bestUpperBound, bestEfficiency,
                        significance,
                        N_sig, N_tot);
        pt->Draw();
        
        // 绘制Pull plot
        pad2->cd();
        TGraphErrors* pullGraph = new TGraphErrors();
        setupPullPlot(pullGraph, config.variable_name, config.fit_min, config.fit_max);
        calculatePull(frame, pullGraph, config.fit_min, config.fit_max);
        pullGraph->GetXaxis()->SetRangeUser(config.fit_min, config.fit_max);
        pullGraph->Draw("AP");
        
        // 存储结果
        if(fitResult->status() == 0) {
            int detecType = detectorTypeMap.at(detector);
            hBe7Frac->SetBinContent(EkBin, detecType+1, fracBe7.getVal());
            hBe7Frac->SetBinError(EkBin, detecType+1, fracBe7.getError());
            hBe9Frac->SetBinContent(EkBin, detecType+1, fracBe9.getVal());
            hBe9Frac->SetBinError(EkBin, detecType+1, fracBe9.getError());
            hBeSgf->SetBinContent(EkBin, detecType+1, significance);
            hChi2NDF->SetBinContent(EkBin, detecType+1, chi2_ndf);
        }
    }
    
    void RunFitting() {
        // 创建PDF文件
        TCanvas* c = new TCanvas("c", "c", 800, 600);
        std::string pdf_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/" + 
                              config.output_prefix + ".pdf";
        c->Print((pdf_path + "[").c_str());
        
        for (int det = 0; det < detectors.size(); ++det) {
            TH2D* issHist = issHists[det];
            TH2D* be7Hist = be7Hists[det];
            TH2D* be9Hist = be9Hists[det];
            TH2D* be10Hist = be10Hists[det];
            
            // 对每个动能bin
            for (int i = 5; i < Be_bins.size()-1; ++i) {
                double ekn_low = Be_bins[i];
                double ekn_high = Be_bins[i+1];
                cout<<"bin:"<<i<<" ekn_low:"<<ekn_low<<" ekn_high:"<<ekn_high<<endl;
                
                // 检查是否在探测器范围内
                if (!isInDetectorRange(detectors[det], (ekn_low + ekn_high)/2)) 
                {
                    cout<<"inrange: "<<detectors[det]<<" "<<(ekn_low + ekn_high)/2<<endl;
                    continue;
                }
                
                // 获取对应的beta bin范围
                //auto beta_bins = getBetaBinRange(issHist, ekn_low, ekn_high);
                //int bin_low = beta_bins.first;
                //int bin_high = beta_bins.second;
                //if(bin_high<bin_low) continue;
                //cout<<"bin:"<<i<<" bin_low:"<<bin_low<<" bin_high:"<<bin_high<<endl;
                
                // 在beta bin范围内投影
                TH1D* hDataProj = issHist->ProjectionX(Form("data_%d_%d", det, i), i+1, i+1);
                TH1D* hBe7Proj = be7Hist->ProjectionY(Form("be7_%d_%d", det, i), i+1, i+1);
                TH1D* hBe9Proj = be9Hist->ProjectionY(Form("be9_%d_%d", det, i), i+1, i+1);
                TH1D* hBe10Proj = be10Hist->ProjectionY(Form("be10_%d_%d", det, i), i+1, i+1);
                
                if (hDataProj->Integral() < 10 || hDataProj->GetMaximum() < 20) {
                    delete hDataProj;
                    delete hBe7Proj;
                    delete hBe9Proj;
                    delete hBe10Proj;
                    continue;
                }
                /*
                hBe7Proj->Smooth(1); 
                hBe9Proj->Smooth(1); 
                hBe10Proj->Smooth(1); 
                */
                
                // 创建画布
                TCanvas* canvas = new TCanvas("canvas", "Template Fit", 800, 600);
                canvas->Divide(1, 2);
                TPad* pad1 = (TPad*)canvas->cd(1);
                pad1->SetPad(0, 0.25, 1, 1);
                TPad* pad2 = (TPad*)canvas->cd(2);
                pad2->SetPad(0, 0, 1, 0.27);
                pad2->SetBottomMargin(0.3);
                pad2->SetGridy();

                
                // 计算实际使用的动能范围
                //double actual_ekn_low = betaToEkn(issHist->GetXaxis()->GetBinLowEdge(bin_low));
                //double actual_ekn_high = betaToEkn(issHist->GetXaxis()->GetBinUpEdge(bin_high));
                
                // 设置标题
                TString frametitle = Form("%s, %.2f < E_{k}/n < %.2f GeV/n",
                                        detectors[det].c_str(),
                                        ekn_low, ekn_high);
                
                // 执行拟合
                //double beta_center = (issHist->GetXaxis()->GetBinCenter(bin_low) + issHist->GetXaxis()->GetBinCenter(bin_high)) / 2;
                FitBin(hDataProj, hBe7Proj, hBe9Proj, hBe10Proj, i+1, detectors[det].c_str(), pad1, pad2, frametitle, pdf_path);
                
                // 保存到PDF
                canvas->Print((pdf_path).c_str());
                delete canvas;
                
                // 清理
                delete hDataProj;
                delete hBe7Proj;
                delete hBe9Proj;
                delete hBe10Proj;
            }
        }
        
        // 关闭PDF
        c->Print((pdf_path + "]").c_str());
        delete c;
        
        // 保存结果
        outFile->cd();
         
        hBe7Frac->Write();
        hBe9Frac->Write();
        hBeSgf->Write();
        hChi2NDF->Write();
    }
};