void B11FragBranchRcp() {
    // 设置样式
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    
    // 打开文件
    TFile* f_ISS = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/ISS_FragmentBer_MassFitAnalysis_unb_perfectQfit.root");
    TFile* f_MC = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_narrow_frag_unbL1_SDIATBeBcut.root");
    
    if (!f_ISS || !f_MC) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }
    
    // 获取ISS结果
    TH1F* h_ISS = (TH1F*)f_ISS->Get("h_Be10_from_B11_ratio");
    if (!h_ISS) {
        std::cerr << "Cannot find h_Be10_from_B11_ratio in ISS file!" << std::endl;
        return;
    }
    
    // 能量bin边界（与主程序一致）
    const double WideBins[] = {0.61, 0.86, 1.17, 1.55, 2.01, 2.57, 3.23, 4.00, 4.91, 5.99, 7.18, 8.60, 10.25, 12.13, 16.38};
    const int nBins = 14;
    
    // 创建MC结果的直方图
    TH1F* h_MC = new TH1F("h_MC_Be10_from_B11_ratio", "Be10 from B11 / Total Be from B11", nBins, WideBins);
    
    // 探测器名称数组
    const char* detectors[] = {"TOF", "NaF", "AGL"};
    
    // 预先读取并Rebin所有探测器的直方图
    std::map<std::string, TH1F*> h_Be7_map, h_Be9_map, h_Be10_map;
    
    for (const char* det : detectors) {
        h_Be7_map[det] = (TH1F*)f_MC->Get(Form("SourceToIsoCounts_%s_Be7", det));
        h_Be9_map[det] = (TH1F*)f_MC->Get(Form("SourceToIsoCounts_%s_Be9", det));
        h_Be10_map[det] = (TH1F*)f_MC->Get(Form("SourceToIsoCounts_%s_Be10", det));
        
        if (!h_Be7_map[det] || !h_Be9_map[det] || !h_Be10_map[det]) {
            std::cerr << "Cannot find MC histograms for " << det << std::endl;
            continue;
        }
        
        // Rebin to match wide bins (factor of 2) - 只做一次
        h_Be7_map[det]->Rebin(2);
        h_Be9_map[det]->Rebin(2);
        h_Be10_map[det]->Rebin(2);
    }
    
    // 计算每个能量bin的MC比值
    for (int bin = 1; bin <= nBins; bin++) {
        double binCenter = h_MC->GetBinCenter(bin);
        
        // 根据能量选择探测器
        const char* det;
        if (binCenter < 1.17) det = "TOF";
        else if (binCenter < 3.23) det = "NaF";
        else det = "AGL";
        
        // 获取已经Rebin过的直方图
        TH1F* h_Be7 = h_Be7_map[det];
        TH1F* h_Be9 = h_Be9_map[det];
        TH1F* h_Be10 = h_Be10_map[det];
        
        if (!h_Be7 || !h_Be9 || !h_Be10) {
            continue;
        }
        
        // 获取对应bin的值（注意bin索引偏移）
        double Be7_val = h_Be7->GetBinContent(bin + 1);
        double Be9_val = h_Be9->GetBinContent(bin + 1);
        double Be10_val = h_Be10->GetBinContent(bin + 1);
        
        double Be7_err = h_Be7->GetBinError(bin + 1);
        double Be9_err = h_Be9->GetBinError(bin + 1);
        double Be10_err = h_Be10->GetBinError(bin + 1);
        
        // 计算总和
        double total_Be = Be7_val + Be9_val + Be10_val;
        
        // 计算比值和误差
        double ratio = 0.0;
        double ratio_err = 0.0;
        
        if (total_Be > 0) {
            ratio = Be10_val / total_Be;
            // 误差传播 - 修正当Be10_val为0时的情况
            if (Be10_val > 0) {
                double Be_total_err = std::sqrt(Be7_err*Be7_err + Be9_err*Be9_err + Be10_err*Be10_err);
                ratio_err = ratio * std::sqrt(pow(Be10_err/Be10_val, 2) + pow(Be_total_err/total_Be, 2));
            } else {
                // 如果Be10为0，误差也设为0
                ratio_err = 0;
            }
        }
        
        h_MC->SetBinContent(bin, ratio);
        h_MC->SetBinError(bin, ratio_err);
        
        std::cout << "Bin " << bin << " (" << WideBins[bin-1] << "-" << WideBins[bin] 
                  << " GeV/n), Detector: " << det 
                  << ", MC Ratio: " << ratio << " +/- " << ratio_err << std::endl;
    }
    
    // 创建画布
    TCanvas* c1 = new TCanvas("c1", "Be10 from B11 Ratio Comparison", 800, 600);
    c1->SetLogx();
    c1->SetGrid();
    
    // 设置直方图样式
    h_ISS->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    h_ISS->GetYaxis()->SetTitle("B11->Be10 / B11->Be");
    h_ISS->GetYaxis()->SetRangeUser(0, 1);
    h_ISS->SetTitle("");
    
    h_MC->SetLineColor(kBlue);
    h_MC->SetMarkerColor(kBlue);
    //h_MC->SetMarkerStyle(21);
    //h_MC->SetMarkerSize(1.0);
    
    // 绘图
    h_ISS->Draw("E");
    h_MC->Draw("E SAME");
    
    // 保存
    c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be10fromB11_comparison.png");
    
    // 清理
    delete h_MC;
    f_ISS->Close();
    f_MC->Close();
    delete f_ISS;
    delete f_MC;
}