#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPad.h>
#include <iostream>

void DrawChi2() {
    // 设置风格
    gStyle->SetOptStat(0);

    // 文件路径
    TString path_AGL = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/result_AGLMass_Ber_WS_NoA.root";
    TString path_NaF = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/result_NaFMass_Ber_WS_NoA.root";
    TString path_ToF = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/result_ToFMass_Ber_WS_NoA.root";

    // 打开文件并读取直方图
    TFile *file_AGL = new TFile(path_AGL);
    TFile *file_NaF = new TFile(path_NaF);
    TFile *file_ToF = new TFile(path_ToF);

    TH2F *h_chi2_AGL = (TH2F*)file_AGL->Get("h_chi2_ndf");
    TH2F *h_chi2_NaF = (TH2F*)file_NaF->Get("h_chi2_ndf");
    TH2F *h_chi2_ToF = (TH2F*)file_ToF->Get("h_chi2_ndf");

    // 创建 PDF 文件以保存所有画布
    TString output_pdf = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/Chi2_Ber_WS_NoA.pdf";

    // 先创建 PDF 文件，并准备写入
    TCanvas *c = new TCanvas("c", "Chi2/NDF vs #alpha", 800, 600);
    c->SetLogy();
    c->Print(output_pdf + "[");

    // 数据文件对应的名字，用于画布标题
    TString mass_titles[3] = { "ToF Mass TempFit", "NaF Mass TempFit", "AGL Mass TempFit" };

    // 对每个文件进行处理
    for (int i = 0; i < 3; ++i) {
        // 选择对应的直方图
        TH2F *h_chi2 = nullptr;
        if (i == 0) h_chi2 = h_chi2_ToF;
        if (i == 1) h_chi2 = h_chi2_NaF;
        if (i == 2) h_chi2 = h_chi2_AGL;

        int n_xbins = h_chi2->GetNbinsX();

        // 遍历每个 Ek bin（x 轴）
        for (int ix = 1; ix <= n_xbins; ++ix) {
            // 获取该 bin 的边界 (GeV/n)
            double binlowedge = h_chi2->GetXaxis()->GetBinLowEdge(ix);
            double binupedge = h_chi2->GetXaxis()->GetBinUpEdge(ix);

            // 投影到 y 轴（alpha），生成一个 1D 直方图
            TH1D *hp = h_chi2->ProjectionY(Form("proj_%d_%d", i, ix), ix, ix);
            cout<<hp->GetXaxis()->GetBinLowEdge(1)<<endl;

            // 找到小于5000的最大值
            double max_chi2 = -1;
            for (int bin = 1; bin <= hp->GetNbinsX(); ++bin) {
                double bin_content = hp->GetBinContent(bin);
                if (bin_content < 5000 && bin_content > max_chi2) {
                    max_chi2 = bin_content;
                }
            }

            // 如果找到的最大值有效，设置y轴上限为该值的1.1倍
            if (max_chi2 > 0) {
                hp->GetYaxis()->SetRangeUser(0.6,max_chi2 * 1.01);
            }
                

            // 创建一个新的画布
            TCanvas *c_temp = new TCanvas(Form("c_temp_%d_%d", i, ix), "Chi2/NDF", 800, 600);
            c_temp->cd();
            c_temp->SetLogy(0);

            // 设置直方图的标题和轴标签
            hp->SetTitle(Form("%s: %.2f - %.2f GeV/n", mass_titles[i].Data(), binlowedge, binupedge));
            hp->GetXaxis()->SetLabelSize(0.05);
            hp->GetXaxis()->SetTitle("#alpha");
            hp->GetYaxis()->SetTitle("Chi2/NDF");
            hp->Draw("P");

            // 保存画布到 PDF 的一页
            c_temp->Print(output_pdf);

            delete c_temp;
            delete hp;
        }
    }

    // 关闭 PDF 文件
    c->Print(output_pdf + "]");

    // 关闭文件
    file_AGL->Close();
    file_NaF->Close();
    file_ToF->Close();

    delete file_AGL;
    delete file_NaF;
    delete file_ToF;
}