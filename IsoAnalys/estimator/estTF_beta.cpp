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

using namespace RooFit;

// 计算每核能量（GeV/n）
double betaToEkn(double beta) {
    const double mass = 0.9315; // GeV, 质子质量
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    return (gamma - 1.0) * mass;
}

// 模板拟合类
class TemplateFitter {
private:
    TFile* issFile;
    TFile* be7File;
    TFile* be9File;
    TFile* be10File;
    TFile* outFile;
    
    std::vector<std::string> detectors = {"TOF", "NaF", "Agl"};
    std::vector<TH2D*> issHists;
    std::vector<TH2D*> be7Hists;
    std::vector<TH2D*> be9Hists;
    std::vector<TH2D*> be10Hists;
    
    // 用于存储结果的直方图
    TH1D* hBe7Frac;
    TH1D* hBe9Frac;
    TH1D* hChi2NDF;

public:
    TemplateFitter() {
        // 打开文件
        issFile = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_iss.root");
        be7File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be7.root");
        be9File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be9.root");
        be10File = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/mdbe/mdfil_Be10.root");
        outFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/estTF_beta.root", "RECREATE");
        
        // 获取直方图
        for (int i = 0; i < 3; ++i) {
            TH2D* issHist = (TH2D*)issFile->Get(Form("hist0%d3", i+1));
            TH2D* be7Hist = (TH2D*)be7File->Get(Form("hist0%d3", i+1));
            TH2D* be9Hist = (TH2D*)be9File->Get(Form("hist0%d3", i+1));
            TH2D* be10Hist = (TH2D*)be10File->Get(Form("hist0%d3", i+1));
            
            issHists.push_back(issHist);
            be7Hists.push_back(be7Hist);
            be9Hists.push_back(be9Hist);
            be10Hists.push_back(be10Hist);
        }
        
        // 创建结果直方图
        hBe7Frac = new TH1D("hBe7Frac", "Be7 Fraction;#beta;Fraction", 100, 0, 1);
        hBe9Frac = new TH1D("hBe9Frac", "Be9 Fraction;#beta;Fraction", 100, 0, 1);
        hChi2NDF = new TH1D("hChi2NDF", "Chi2/NDF;#beta;Chi2/NDF", 100, 0, 1);
    }
    
    ~TemplateFitter() {
        delete issFile;
        delete be7File;
        delete be9File;
        delete be10File;
        delete outFile;
    }
    
    void FitBin(TH1D* hData, TH1D* hBe7, TH1D* hBe9, TH1D* hBe10,
                double beta, const char* detector, TPad* pad1, TPad* pad2, TString frametitle) {
        // 创建变量
        RooRealVar x("x", "Estimator", 0.09, 0.71);
        
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
        
        // 执行拟合
        RooFitResult* fitResult = totalPdf.fitTo(dataHist, Save());
        
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
        
        // 添加图例
        TLegend* legend = new TLegend(0.2, 0.68, 0.45, 0.87);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->AddEntry(frame->findObject("Data"), "ISS Data", "p");
        legend->AddEntry(frame->findObject("Be7"), "Be7", "l");
        legend->AddEntry(frame->findObject("Be9"), "Be9", "l");
        legend->AddEntry(frame->findObject("Be10"), "Be10", "l");
        legend->Draw();
        
        // 添加拟合信息
        TPaveText* pt = new TPaveText(0.35, 0.68, 0.7, 0.87, "NDC");
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->AddText(Form("Chi2/NDF = %.2f", chi2_ndf));
        pt->AddText(Form("Be7 Frac = %.3f #pm %.3f", fracBe7.getVal(), fracBe7.getError()));
        pt->AddText(Form("Be9 Frac = %.3f #pm %.3f", fracBe9.getVal(), fracBe9.getError()));
        pt->Draw();
        
        // 绘制Pull plot
        pad2->cd();
        TGraphErrors* pullGraph = new TGraphErrors();
        pullGraph->SetTitle(";Estimator;Pull");
        
        RooHist* dataHist_frame = frame->getHist("Data");
        RooCurve* modelCurve = frame->getCurve("TotalModel");
        
        if (dataHist_frame && modelCurve) {
            double x_val, y_data;
            for (int i = 0; i < dataHist_frame->GetN(); ++i) {
                dataHist_frame->GetPoint(i, x_val, y_data);
                
                // 检查x值是否在合理范围内
                if (x_val >= 0.1 && x_val <= 0.65) {  // 根据实际范围调整
                    double y_model = modelCurve->Eval(x_val);
                    double data_err = dataHist_frame->GetErrorY(i);
                    double pull_err = 0;
                    
                    if (y_model != 0) {
                        double pull_val = y_data / y_model;
                        if (y_data != 0) {
                            pull_err = pull_val * (data_err / y_data);
                        }
                        pullGraph->SetPoint(i, x_val, pull_val);
                        pullGraph->SetPointError(i, 0, pull_err);
                    } else {
                        pullGraph->SetPoint(i, x_val, 0);
                        pullGraph->SetPointError(i, 0, 0);
                    }
                }
            }
        }
        
        pullGraph->SetMarkerStyle(20);
        //pullGraph->GetYaxis()->SetRangeUser(-5, 5);
        pullGraph->Draw("AP");
        
        // 存储结果
        int bin = hBe7Frac->FindBin(beta);
        hBe7Frac->SetBinContent(bin, fracBe7.getVal());
        hBe7Frac->SetBinError(bin, fracBe7.getError());
        hBe9Frac->SetBinContent(bin, fracBe9.getVal());
        hBe9Frac->SetBinError(bin, fracBe9.getError());
        hChi2NDF->SetBinContent(bin, chi2_ndf);
    }
    
    void RunFitting() {
        // 创建PDF文件
        TCanvas* c = new TCanvas("c", "c", 800, 800);
        c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/estTF_beta.pdf[");
        
        for (size_t det = 0; det < detectors.size(); ++det) {
            TH2D* issHist = issHists[det];
            TH2D* be7Hist = be7Hists[det];
            TH2D* be9Hist = be9Hists[det];
            TH2D* be10Hist = be10Hists[det];
            
            // 对每个beta bin进行拟合
            int nBinsX = issHist->GetNbinsX();
            for (int i = 1; i <= nBinsX; ++i) {
                double beta = issHist->GetXaxis()->GetBinCenter(i);
                double ekn = betaToEkn(beta);
                
                // 投影得到1D直方图
                TH1D* hDataProj = issHist->ProjectionY(Form("data_%d", i), i, i);
                TH1D* hBe7Proj = be7Hist->ProjectionY(Form("be7_%d", i), i, i);
                TH1D* hBe9Proj = be9Hist->ProjectionY(Form("be9_%d", i), i, i);
                TH1D* hBe10Proj = be10Hist->ProjectionY(Form("be10_%d", i), i, i);
                
                if (hDataProj->Integral() < 50 || hDataProj->GetMaximum() < 10) continue; // 跳过统计太少的bin
                
                // 创建画布
                TCanvas* canvas = new TCanvas("canvas", "Template Fit", 800, 800);
                canvas->Divide(1, 2);
                
                TPad* pad1 = (TPad*)canvas->cd(1);
                pad1->SetPad(0, 0.3, 1, 1);
                
                TPad* pad2 = (TPad*)canvas->cd(2);
                pad2->SetPad(0, 0, 1, 0.3);
                pad2->SetBottomMargin(0.2);
                pad2->SetGridy();
                
                // 设置标题
                pad1->cd();
                TString frametitle = Form("%s, %.2f < E_{k}/n < %.2f GeV/n",
                                   detectors[det].c_str(),
                                   betaToEkn(issHist->GetXaxis()->GetBinLowEdge(i)),
                                   betaToEkn(issHist->GetXaxis()->GetBinUpEdge(i)));
                
                // 执行拟合
                FitBin(hDataProj, hBe7Proj, hBe9Proj, hBe10Proj, beta, 
                       detectors[det].c_str(), pad1, pad2, frametitle);
                
                // 保存到PDF
                canvas->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/estTF_beta.pdf");
                delete canvas;
                
                // 清理
                delete hDataProj;
                delete hBe7Proj;
                delete hBe9Proj;
                delete hBe10Proj;
            }
        }
        
        // 关闭PDF
        c->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoEstimate/estTF_beta.pdf]");
        delete c;
        
        // 保存结果
        outFile->cd();
        hBe7Frac->Write();
        hBe9Frac->Write();
        hChi2NDF->Write();
    }
};

void estTF_beta() {
    TemplateFitter fitter;
    fitter.RunFitting();
}