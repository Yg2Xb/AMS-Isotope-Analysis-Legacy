#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <TPad.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <map>
#include <memory> // For std::unique_ptr

// --- 辅助函数 ---

/**
 * @brief 为直方图设置统一的绘图样式
 * @param hist 指向TH1F对象的指针
 * @param color 颜色
 * @param markerStyle 点的样式
 * @param markerSize 点的大小
 */
void setHistStyle(TH1F* hist, int color, int markerStyle, double markerSize = 1.2) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetLineWidth(2);
    hist->SetStats(0);
}

/**
 * @brief 从文件中安全地获取直方图
 * @param file 指向TFile对象的指针
 * @param histName 直方图的名称
 * @return 返回直方图的指针，如果失败则返回nullptr
 */
TH1F* getHistogram(TFile* file, const TString& histName) {
    if (!file || file->IsZombie()) {
        std::cerr << "错误: 文件指针无效!" << std::endl;
        return nullptr;
    }
    TH1F* hist = nullptr;
    file->GetObject(histName, hist);
    if (!hist) {
        std::cerr << "错误: 在文件 " << file->GetName() << " 中未找到直方图: " << histName << std::endl;
        return nullptr;
    }
    // 创建一个克隆体，这样即使关闭文件，直方图依然可用
    return (TH1F*)hist->Clone(TString::Format("%s_clone", hist->GetName()));
}


// --- 主分析函数 ---
void compareB11Fragmentation() {
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    // --- 1. 配置输入和输出 ---
    // 更新为新的文件路径
    TString issFilePath = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/ISS_FragmentBer_MassFitAnalysis_unb_perfectQfit.root";
    TString mcB11FilePath = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_narrow_frag_unbL1_SDIATBeBcut.root";
    
    // 定义需要比较的直方图名称模板
    TString issHistName = "BeIso_B11Source_Ratio_Be10"; // 10Be from 11B / 11B
    TString mcHistNameTemplate = "SourceToIsoRatio_%s_Be10"; // %s will be replaced by detector name

    // 输出文件
    TString outputFileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Comparison_ISS_vs_B11MC.png";

    // --- 2. 加载数据 ---
    auto issFile = std::unique_ptr<TFile>(TFile::Open(issFilePath));
    auto mcB11File = std::unique_ptr<TFile>(TFile::Open(mcB11FilePath));

    // 加载ISS数据 (WideBin)
    TH1F* issHist = getHistogram(issFile.get(), issHistName);
    if (!issHist) return; // 如果加载失败则退出

    // 加载MC数据 (NarrowBin)，按探测器分开
    std::map<TString, TH1F*> mcDetHists;
    std::vector<TString> detectors = {"TOF", "NaF", "AGL"};
    for (const auto& det : detectors) {
        TString histName = TString::Format(mcHistNameTemplate.Data(), det.Data());
        mcDetHists[det] = getHistogram(mcB11File.get(), histName);
    }

    // --- 3. 处理MC Binning差异 ---
    // 创建一个新的MC直方图，其binning与ISS完全相同
    TH1F* mcFinalHist = (TH1F*)issHist->Clone("mcFinalHist");
    mcFinalHist->Reset(); // 清空内容，只保留结构

    // 对每个探测器的MC直方图进行Rebin
    for (auto const& [det, hist] : mcDetHists) {
        if (hist) {
            hist->Rebin(2);
            hist->Scale(0.5);

        }
    }

    // 填充最终的MC直方图
    // 循环遍历最终图的每一个bin (WideBins)
    for (int i = 1; i <= mcFinalHist->GetNbinsX(); ++i) {
        double energy = mcFinalHist->GetBinCenter(i);
        
        // 根据能量确定使用哪个探测器的数据
        TString currentDet;
        if (energy < 1.17) currentDet = "TOF";
        else if (energy < 3.23) currentDet = "NaF";
        else currentDet = "AGL";

        // 从对应的rebinned MC直方图中获取数据
        TH1F* sourceMCHist = mcDetHists[currentDet];
        if (sourceMCHist) {
            // **核心逻辑**: ISS的第i个bin对应rebinned MC的第(i+1)个bin
            int sourceBin = i + 1;
            
            if (sourceBin <= sourceMCHist->GetNbinsX()) {
                double content = sourceMCHist->GetBinContent(sourceBin);
                double error = sourceMCHist->GetBinError(sourceBin);
                mcFinalHist->SetBinContent(i, content);
                mcFinalHist->SetBinError(i, error);
                cout<<mcFinalHist->GetBinLowEdge(i)<<endl;
                cout<<mcFinalHist->GetBinContent(i)<<endl;
            }
        }
    }
    
    // --- 4. 绘图 ---
    TCanvas* canvas = new TCanvas("c_comparison", "Fragmentation Comparison: ISS vs B11 MC", 900, 700);
    TPad* pad = new TPad("pad", "pad", 0, 0, 1, 1);
    pad->SetGridx();
    pad->SetGridy();
    pad->Draw();
    pad->cd();

    // 设置坐标轴和标题
    issHist->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    issHist->GetYaxis()->SetTitle("Fragmented ^{10}Be / Source ^{11}B Ratio");
    issHist->GetXaxis()->SetRangeUser(0.6, 16.5);
    issHist->GetYaxis()->SetRangeUser(0.0, 0.02);
    issHist->SetTitle("Fragmentation Comparison");

    // 设置样式
    setHistStyle(issHist, kBlack, 20);
    setHistStyle(mcFinalHist, kRed, 20); // 使用不同样式以区分

    // 绘制
    issHist->Draw("PZ"); // "PZ" 绘制点和误差，避免0值被画出
    mcFinalHist->Draw("PZ SAME");

    // 创建图例
    TLegend* legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(issHist, "ISS Data", "ep");
    legend->AddEntry(mcFinalHist, "B11 MC Simulation", "ep");
    legend->Draw();
    
    // --- 5. 保存和清理 ---
    canvas->SaveAs(outputFileName);
    std::cout << "绘图已保存至: " << outputFileName << std::endl;

    // 释放克隆的直方图内存
    delete issHist;
    delete mcFinalHist;
    for (auto const& [det, hist] : mcDetHists) {
        delete hist;
    }
}

// 主函数入口点
void diffQcp() {
    try {
        compareB11Fragmentation();
    } catch (const std::exception& e) {
        std::cerr << "程序运行中发生异常: " << e.what() << std::endl;
    }
}