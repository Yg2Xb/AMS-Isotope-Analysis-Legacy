#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <string>
#include <iostream>

// --- 主分析函数 ---
void DrawTempR() {
    // --- 0. 全局样式和配置 ---
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111); 
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.2);

    // 定义输入输出文件路径
    const std::string basePath = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/";
    const std::string inputFileName = basePath + "chargeFit_B_EGE_detfuncT_unb.root";
    const std::string outputPdfName = basePath + "chargeFit_B_EGE_detfuncT_unb_ToyFitDetails.pdf";
    const std::string outputPngName = basePath + "chargeFit_B_EGE_detfuncT_unb_ToyComparison.png";
    const std::string outputRootName = basePath + "chargeFit_B_EGE_detfuncT_unb_ToyResults.root";

    // --- 1. 读取输入文件和直方图 ---
    auto inputFile = std::make_unique<TFile>(inputFileName.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "错误: 无法打开输入文件 " << inputFileName << std::endl;
        return;
    }

    auto h2_toy_be = (TH2D*)inputFile->Get("h2_toy_narrow_frac_be");
    auto h_be_fraction_orig = (TH1D*)inputFile->Get("h_be_fraction");

    if (!h2_toy_be || !h_be_fraction_orig) {
        std::cerr << "错误: 无法从文件中读取所需的直方图。" << std::endl;
        return;
    }
    std::cout << "成功读取输入文件和直方图。" << std::endl;

    // --- 2. 准备输出的直方图和多页PDF ---
    auto h_mean = (TH1D*)h_be_fraction_orig->Clone("h_mean_from_toy");
    h_mean->SetTitle("Mean of Be Fraction from Toy MC;E_{k} [GeV/n];Mean Fraction");
    h_mean->Reset();

    auto h_sigma = (TH1D*)h_be_fraction_orig->Clone("h_sigma_from_toy");
    h_sigma->SetTitle("Sigma of Be Fraction from Toy MC;E_{k} [GeV/n];Sigma (Uncertainty)");
    h_sigma->Reset();
    
    TCanvas *c_fits = new TCanvas("c_fits", "Toy Fit Details", 800, 600);
    c_fits->Print((outputPdfName + "[").c_str()); // 开启多页PDF

    // --- 3. 循环遍历每个能量仓，进行高斯拟合 ---
    int n_bins_y = h2_toy_be->GetYaxis()->GetNbins();
    for (int i = 1; i <= n_bins_y; ++i) {
        double ek_low = h2_toy_be->GetYaxis()->GetBinLowEdge(i);
        double ek_high = h2_toy_be->GetYaxis()->GetBinUpEdge(i);
        
        // 为每个能量仓的投影创建唯一的名称
        std::string proj_name = "proj_bin_" + std::to_string(i);
        TH1D* proj = h2_toy_be->ProjectionX(proj_name.c_str(), i, i);
        //proj->Rebin(2);
        //if(i >= 6 && i <= 13) proj->Rebin(4);
        //if(i > 13) proj->Rebin(2);
        if(i >= 1 && i < 6)proj->GetXaxis()->SetRangeUser(0,0.005);
        if(i >= 6)proj->GetXaxis()->SetRangeUser(0,0.005);

        if (proj->GetEntries() < 50) { // 如果统计量太少，则跳过拟合
            delete proj;
            continue;
        }
        

        // 创建高斯函数并设置初始参数以帮助拟合收敛
        TF1 *gaus_fit = new TF1("gaus_fit", "gaus", 0, 0.005);
        gaus_fit->SetParameter(1, h_be_fraction_orig->GetBinContent(i)); // 初始mean设为直方图的均值
        gaus_fit->SetParLimits(1, 0.8*h_be_fraction_orig->GetBinContent(i), 1.2*h_be_fraction_orig->GetBinContent(i)); // 初始mean设为直方图的均值

        gaus_fit->SetParameter(2, proj->GetStdDev()); // 初始sigma设为直方图的标准差

        // 执行拟合
        proj->Fit(gaus_fit, "RQ"); // "R"使用函数范围, "Q"安静模式

        // 提取并存储结果
        double mean = gaus_fit->GetParameter(1);
        double sigma = gaus_fit->GetParameter(2);
        double chi2 = gaus_fit->GetChisquare();
        int ndf = gaus_fit->GetNDF();
        double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;

        h_mean->SetBinContent(i, mean);
        h_mean->SetBinError(i, gaus_fit->GetParError(1));
        h_sigma->SetBinContent(i, sigma);
        h_sigma->SetBinError(i, gaus_fit->GetParError(2));

        // 绘制并保存到PDF
        proj->SetTitle(Form("Toy Be Fraction for E_{k} in [%.2f, %.2f] GeV/n", ek_low, ek_high));
        proj->GetXaxis()->SetTitle("Be Fraction");
        proj->GetYaxis()->SetTitle("Toy Experiments");
        proj->Draw("hist");
        gaus_fit->SetLineColor(kRed); 
        gaus_fit->Draw("same");
        
        TPaveText *pt = new TPaveText(0.6, 0.5, 0.9, 0.7, "NDC");
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(12);
        pt->AddText(Form("#mu = %.4f #pm %.4f", mean, gaus_fit->GetParError(1)));
        pt->AddText(Form("#sigma = %.4f #pm %.4f", sigma, gaus_fit->GetParError(2)));
        pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2_ndf));
        //pt->Draw();

        c_fits->Print(outputPdfName.c_str());

        delete proj;
        delete gaus_fit;
        delete pt;
    }
    c_fits->Print((outputPdfName + "]").c_str()); // 关闭多页PDF
    delete c_fits;
    std::cout << "所有能量仓拟合完成，详情已保存至 " << outputPdfName << std::endl;

    // --- 4. 绘制最终对比图 ---
    TCanvas *c_comp = new TCanvas("c_comp", "Comparison of Fit and Toy MC Results", 800, 600);
    c_comp->SetLogx();
    c_comp->SetGrid();

    // 从Toy MC结果创建TGraphErrors
    auto g_toy_results = new TGraphErrors();
    g_toy_results->SetTitle("Toy MC Results");
    int n_points = 0;
    for (int i = 1; i <= h_mean->GetNbinsX(); ++i) {
        //if (h_mean->GetBinContent(i) > 0) {
            double ek_center = h_mean->GetXaxis()->GetBinCenter(i);
            double ek_width = h_mean->GetXaxis()->GetBinWidth(i) / 2.0;
            double mean_val = h_mean->GetBinContent(i);
            cout<<h_mean->GetXaxis()->GetBinLowEdge(i)<<" "<<ek_center<<" "<<mean_val<<endl;
            double sigma_val = h_sigma->GetBinContent(i); // sigma是误差
            g_toy_results->SetPoint(n_points, ek_center, mean_val);
            g_toy_results->SetPointError(n_points, 0, sigma_val);
            n_points++;
        //}
    }
    
    // 设置原始拟合结果的样式
    h_be_fraction_orig->SetLineColor(kRed);
    h_be_fraction_orig->SetMarkerColor(kRed);
    h_be_fraction_orig->SetMarkerStyle(20);
    h_be_fraction_orig->GetYaxis()->SetTitleOffset(1.5);
    h_be_fraction_orig->GetXaxis()->SetNdivisions(505);
    h_be_fraction_orig->GetYaxis()->SetTitle("Be / (Be+B+C) Fraction");
    h_be_fraction_orig->SetTitle("Cut InnerQ 3.45-5.45, Be Fraction in unbL1Q 5.0-5.4");
    h_be_fraction_orig->Draw("p"); // 红点和线，带误差
    //h_be_fraction_orig->GetYaxis()->SetRangeUser(0,0.2);

    // 设置Toy MC结果的样式
    g_toy_results->SetLineColor(kBlue);
    g_toy_results->SetMarkerColor(kBlue);
    g_toy_results->SetMarkerStyle(21);
    g_toy_results->Draw("pz SAME"); // 蓝色点和线，带误差棒
    h_be_fraction_orig->Draw("pz same"); // 红点和线，带误差

    TLegend *leg = new TLegend(0.15, 0.67, 0.45, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_be_fraction_orig, "Fit Result", "pz");
    leg->AddEntry(g_toy_results, "Toy MC Result (#mu_{guas}#pm#sigma_{guas})", "pz");
    leg->Draw();

    c_comp->Print(outputPngName.c_str());
    std::cout << "最终对比图已保存至 " << outputPngName << std::endl;
    delete c_comp;

    // --- 5. 保存结果到新的ROOT文件 ---
    auto outputFile = std::make_unique<TFile>(outputRootName.c_str(), "RECREATE");
    h_mean->Write();
    h_sigma->Write();
    outputFile->Close();
    std::cout << "拟合结果 (mean, sigma) 已保存至 " << outputRootName << std::endl;
    
    inputFile->Close();
}