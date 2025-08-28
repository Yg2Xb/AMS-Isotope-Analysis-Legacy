#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
// 拟合后添加置信区间
void DrawConfidenceBands(TGraphAsymmErrors* data_graph, TF1* fit_func) {
    // 获取协方差矩阵
    TFitResultPtr fit_result = data_graph->Fit(fit_func, "SQ");
    TMatrixD cov_matrix = fit_result->GetCovarianceMatrix();
    
    // 获取x轴范围
    double xmin = data_graph->GetXaxis()->GetXmin();
    double xmax = data_graph->GetXaxis()->GetXmax();
    
    // 创建用于误差带的图形
    const int npoints = 100;
    double dx = (xmax - xmin) / npoints;
    
    std::vector<double> x_points(npoints+1);
    std::vector<double> y_upper(npoints+1);
    std::vector<double> y_lower(npoints+1);
    
    for (int i = 0; i <= npoints; i++) {
        double x = xmin + i * dx;
        double y = fit_func->Eval(x);
        double error = 0;
        
        // 计算误差传播
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                error += cov_matrix(i1,i2) * pow(x,i1) * pow(x,i2);
            }
        }
        error = 2.5*sqrt(error);
        
        x_points[i] = x;
        y_upper[i] = y + error;
        y_lower[i] = y - error;
    }
    
    // 创建上下限曲线
    TGraph* upper_band = new TGraph(npoints+1, x_points.data(), y_upper.data());
    TGraph* lower_band = new TGraph(npoints+1, x_points.data(), y_lower.data());
    
    // 设置曲线样式
    upper_band->SetLineColor(kRed);
    upper_band->SetLineStyle(2);
    lower_band->SetLineColor(kRed);
    lower_band->SetLineStyle(2);
    
    // 绘制误差带
    upper_band->Draw("same L");
    lower_band->Draw("same L");
    
}

void FitBestAlpha() {
    // 输入文件路径
    const char* file_paths[] = {
        "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/MassTF_Be_DiffBin_r1_Alpha_Use7_bestalpha.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/MassTF_Be_DiffBin_r1_Alpha_Use9_bestalpha.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/MassTF_Be_DiffBin_r1_Alpha_Use10_bestalpha.root"
    };

    // 定义颜色和标签
    const int colors[] = {kBlue, kGreen + 1, kOrange + 2};
    const char* labels[] = {"Mass7 bin", "Mass9 bin", "Mass10 bin"};

    // 定义探测器名称
    const char* detectors[] = {"ToF", "NaF", "AGL"};

    // 遍历每个探测器
    for (const auto& det : detectors) {
        // 创建画布
        TCanvas* c = new TCanvas(Form("c_%s", det), Form("%s Best Alpha", det), 800, 600);
        c->SetGrid();

        // 创建图例
        TLegend* leg = new TLegend(0.65, 0.16, 0.88, 0.35);
        leg->SetBorderSize(0); // 无边框
        leg->SetFillStyle(0);  // 无色透明

        // 用于合并所有数据的 TGraphAsymmErrors
        TGraphAsymmErrors* merged_graph = new TGraphAsymmErrors();

        // 遍历每个文件
        for (int i = 0; i < 3; ++i) {
            TFile* file = TFile::Open(file_paths[i]);
            if (!file || file->IsZombie()) {
                std::cerr << "Error: Cannot open file " << file_paths[i] << std::endl;
                continue;
            }

            // 获取 TGraphAsymmErrors
            TGraphAsymmErrors* graph = dynamic_cast<TGraphAsymmErrors*>(file->Get(Form("g_err_%s", det)));
            if (!graph) {
                std::cerr << "Error: Cannot find graph g_err_" << det << " in file " << file_paths[i] << std::endl;
                file->Close();
                continue;
            }

            // 设置图形属性
            graph->SetMarkerColor(colors[i]);
            graph->SetLineColor(colors[i]);
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(1.2);

            // 绘制图形
            if (i == 0) {
                graph->Draw("AP");
                graph->GetYaxis()->SetRangeUser(0.98, 1.02); // 设置 Y 轴范围
                graph->SetTitle(Form("%s Best Alpha;E_{k}/n [GeV/n];#alpha", det));
            } else {
                graph->Draw("P SAME");
            }

            // 添加图例
            leg->AddEntry(graph, labels[i], "p");

            // 合并数据
            for (int j = 0; j < graph->GetN(); ++j) {
                double x, y, ex_low, ex_high, ey_low, ey_high;
                graph->GetPoint(j, x, y);
                ex_low = graph->GetErrorXlow(j);
                ex_high = graph->GetErrorXhigh(j);
                ey_low = graph->GetErrorYlow(j);
                ey_high = graph->GetErrorYhigh(j);
                int n = merged_graph->GetN();
                merged_graph->SetPoint(n, x, y);
                merged_graph->SetPointError(n, ex_low, ex_high, ey_low, ey_high);
            }

            file->Close();
        }

        // 绘制图例
        leg->Draw();

        // 拟合合并后的数据
        TF1* fit_func = new TF1("fit_func", "pol2", 0, 20); // 二次多项式拟合
        merged_graph->Fit(fit_func, "Q"); // 安静模式拟合
        fit_func->SetLineWidth(3);
        fit_func->SetLineColor(kRed);
        fit_func->Draw("same");
        
        DrawConfidenceBands(merged_graph, fit_func);

        double chi2 = fit_func->GetChisquare();
        int ndf = fit_func->GetNDF();
        double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;

        // 在画布上显示拟合结果
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.2, 0.8, Form("Fit: #alpha = %.4f + %.4f E_{k}/n + %.4f E_{k}^{2}/n",//, + %.4f E_{k}^{3}", 
                                         fit_func->GetParameter(0), 
                                         fit_func->GetParameter(1), 
                                         fit_func->GetParameter(2)));
                                         //fit_func->GetParameter(3))
        latex.DrawLatex(0.2, 0.75, Form("Errors: %.4f, %.4f, %.4f",//, %.4f", 
                                         fit_func->GetParError(0), 
                                         fit_func->GetParError(1), 
                                         fit_func->GetParError(2)));
                                         //fit_func->GetParError(3)
        latex.DrawLatex(0.2, 0.70, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf, chi2_ndf));

        // 保存画布为 PNG
        c->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/%s_BestAlpha_Fit.png", det));

        // 清理
        delete merged_graph;
        delete fit_func;
        delete leg;
        delete c;
    }
}