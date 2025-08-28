#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TLatex.h>
#include <vector>
#include "../Tool.h"

using namespace AMS_Iso;

struct DetectorConfig {
    const char* name;
    double ek_min;
    double ek_max;
};

void process(const char* input_path, int charge, int mass) {
    // 初始化配置
    const DetectorConfig detectors[] = {
        {"ToF", 0.22, 1.5},
        {"NaF", 0.8, 3.5},
        {"AGL", 2.3, 15.0}
    };

    // 获取输入文件前缀
    TString basename = gSystem->BaseName(input_path);
    TString prefix = basename(basename.Index("Mass"), basename.Length()-5);

    // 初始化输入输出文件
    TFile* fin = TFile::Open(input_path);
    TFile* fout = new TFile(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/%s_bestalpha.root", prefix.Data()), "RECREATE");
    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->Print(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/%s_Chi2.pdf[", prefix.Data()));

    // 获取能量分bin
    const auto& ek_bins = getKineticEnergyBins(charge, mass);

    for (const auto& det : detectors) {
        // 获取二维直方图
        TH2F* h2 = dynamic_cast<TH2F*>(fin->Get(Form("h_chi2_%s", det.name)));
        cout<<h2->GetMaximum()<<endl;
        if (!h2) continue;

        // 创建输出数据结构
        TH1D* h_best = new TH1D(Form("h_best_%s", det.name), Form("%s Best Alpha;E_{k}/n [GeV/n];#alpha", det.name), ek_bins.size()-1, ek_bins.data());
        TGraphAsymmErrors* g_err = new TGraphAsymmErrors();

        // 处理每个能量bin
        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double ek_center = h2->GetXaxis()->GetBinCenter(bx);
            if (ek_center < det.ek_min || ek_center > det.ek_max) continue;

            // 投影得到alpha分布
            TH1D* h_proj = h2->ProjectionY(Form("proj_%s_%d", det.name, bx), bx, bx);
            double max_chi2 = h_proj->GetMaximum(); 
            // 寻找最小值
            double min_chi2 = 1e9, best_alpha = 0;
            int min_bin = -1;
            for (int by = 1; by <= h_proj->GetNbinsX(); ++by) {
                double chi2 = h_proj->GetBinContent(by);
                if (chi2 > 0 && chi2 < min_chi2) {
                    min_chi2 = chi2;
                    best_alpha = h_proj->GetXaxis()->GetBinCenter(by);
                    min_bin = by; 
                }
            }

            // 寻找误差范围
            double alpha_low = best_alpha, alpha_high = best_alpha;
            for (int by = min_bin-1 ; by >= 1; --by) {
                if (h_proj->GetBinContent(by) > min_chi2+1) {
                    alpha_low = h_proj->GetXaxis()->GetBinCenter(by);
                    break;
                }
            }
            for (int by = min_bin+1; by <= h_proj->GetNbinsX(); ++by) {
                if (h_proj->GetBinContent(by) > min_chi2+1) {
                    alpha_high = h_proj->GetXaxis()->GetBinCenter(by);
                    break;
                }
            }
            //Printf("Ek = %.2f, best alpha = %.3f, alpha_low = %.3f, alpha_high = %.3f",h2->GetXaxis()->GetBinLowEdge(bx), best_alpha, alpha_low, alpha_high);

            // 存储结果
            h_best->Fill(ek_center, best_alpha);
            int n = g_err->GetN();
            g_err->SetPoint(n, ek_center, best_alpha);
            g_err->SetPointError(n, 0, 0, best_alpha-alpha_low, alpha_high-best_alpha);
            g_err->SetMarkerStyle(20);
            g_err->SetMarkerSize(0.8);

            // 绘图
            h_proj->SetTitle(Form("%s, %.2f - %.2f GeV/n", det.name, h2->GetXaxis()->GetBinLowEdge(bx), h2->GetXaxis()->GetBinUpEdge(bx)));
            h_proj->SetMinimum(0.85*min_chi2);
            h_proj->GetXaxis()->SetTitle("#alpha");
            h_proj->GetYaxis()->SetTitle("#chi^{2}");
            h_proj->Draw("P");
            // 绘制红色虚横线表示 min_chi2 + 1
            TLine* line = new TLine(h_proj->GetXaxis()->GetXmin(), min_chi2 + 1, h_proj->GetXaxis()->GetXmax(), min_chi2 + 1);
            line->SetLineColor(kRed); // 设置线条颜色为红色
            line->SetLineStyle(2);    // 设置线条样式为虚线
            line->SetLineWidth(2);    // 设置线条宽度
            line->Draw("same");       // 在同一画布上绘制
            c->Print(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/%s_Chi2.pdf", prefix.Data()));
            
            delete h_proj;
        }

        // 写入输出文件
        fout->cd();
        h_best->Write();
        g_err->Write(Form("g_err_%s", det.name));
        delete h_best;
        delete g_err;
    }

    // 清理资源
    c->Print(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/%s_Chi2.pdf]", prefix.Data()));
    delete c;
    fout->Close();
    fin->Close();
}

void Chi2AndAlpha(){
    process("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MassTF_Be_DiffBin_r1_Alpha_Use7.root",4,7);
    process("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MassTF_Be_DiffBin_r1_Alpha_Use9.root",4,9);
    process("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MassTF_Be_DiffBin_r1_Alpha_Use10.root",4,10);
}