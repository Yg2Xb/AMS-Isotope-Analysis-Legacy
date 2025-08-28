#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TColor.h>

void TofBorCP_Mean()
{
    // 打开ROOT文件
    TFile *file1 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/BerIsoResults/Beta/FitTOFinEkBin_Ber.root");
    TFile *file3 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/BerIsoResults/Beta/FitTOFinEkBin_Be10.root");

    // 从每个文件中提取直方图
    TH1D *hist1_NaF = (TH1D *)file1->Get("TOF_NaF_Mean_Ber")->Clone();
    TH1D *hist1_AGL = (TH1D *)file1->Get("TOF_AGL_Mean_Ber")->Clone();
    TH1D *hist3_NaF = (TH1D *)file3->Get("TOF_NaF_Mean_Be10")->Clone();
    TH1D *hist3_AGL = (TH1D *)file3->Get("TOF_AGL_Mean_Be10")->Clone();

    // 将NaF和AGL的bins合并到新直方图
    for (int i = 1; i <= hist1_NaF->GetNbinsX(); ++i)
    {
        if (i <= 100)
        {
            hist1_NaF->SetBinContent(i, hist1_NaF->GetBinContent(i));
            hist3_NaF->SetBinContent(i, hist3_NaF->GetBinContent(i));
            hist1_NaF->SetBinError(i, hist1_NaF->GetBinError(i));
            hist3_NaF->SetBinError(i, hist3_NaF->GetBinError(i));
        }
        else
        {
            hist1_NaF->SetBinContent(i, hist1_AGL->GetBinContent(i));
            hist3_NaF->SetBinContent(i, hist3_AGL->GetBinContent(i));
            hist1_NaF->SetBinError(i, hist1_AGL->GetBinError(i));
            hist3_NaF->SetBinError(i, hist3_AGL->GetBinError(i));
        }
    }
        // 创建画布
        TCanvas *c1 = new TCanvas("c1", "Combined Histograms", 800, 600);
        c1->cd();
        c1->SetLogx();
        //c1->SetGridx();
        //c1->SetGridy();

        hist3_NaF->SetMarkerColor(kBlue);
        hist3_NaF->SetMarkerSize(.8);
        hist1_NaF->SetMarkerSize(.8);
        hist3_NaF->SetLineColor(kBlue);
        // 绘制直方图
        c1->SetLeftMargin(.22);
        hist1_NaF->GetYaxis()->SetLabelSize(0.06);
        hist1_NaF->GetYaxis()->SetTitleSize(0.06);
        hist1_NaF->GetXaxis()->SetLabelSize(0.06);
        hist1_NaF->GetXaxis()->SetTitleSize(0.06);
        hist1_NaF->GetXaxis()->SetLabelOffset(-0.01);
        hist1_NaF->GetYaxis()->SetRangeUser(-0.02, 0.005);
        //hist1_NaF->GetYaxis()->SetRangeUser(0.002, 0.016);
        //hist1_NaF->GetXaxis()->SetRangeUser(0.5, 10);
        hist1_NaF->GetYaxis()->SetTitleOffset(1.6);
        hist1_NaF->GetYaxis()->SetNdivisions(508);
        hist1_NaF->SetTitle("");
        hist1_NaF->GetXaxis()->SetTitle("NaF Ek/n [GeV/n]");
        hist1_NaF->GetYaxis()->SetTitle("Mean of 1/Beta_{TOF}-1/Beta_{NaF}");
        hist1_NaF->Draw("e");
        hist3_NaF->Draw("eSAME");

        // 添加图s
/*
        TLegend *leg = new TLegend(0.26, 0.65, 0.66, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->Draw("same");
*/
        TLegend *leg2 = new TLegend(0.65, 0.2, 0.85, 0.35);
        leg2->SetBorderSize(0);
        leg2->SetFillColor(0);
        leg2->SetFillStyle(0);
        leg2->AddEntry(hist3_NaF, "MC", "P");
        leg2->AddEntry(hist1_NaF, "ISS", "P");
        leg2->Draw("same");
        c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/BerIsoResults/Beta/LitCPMean.png");
    }
