// 设置全局样式
void setStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTextSize(0.05);
}

// 创建并配置画布
TCanvas* createCanvas() {
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
    return canvas;
}

// 设置直方图样式
void setHistStyle(TH1* h, const char* xtitle, const char* ytitle, double xmin, double xmax, 
                 int color = kBlack, int markerStyle = 20) {
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(0.8);
    h->SetMarkerColor(color);
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    //h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleOffset(1.5);
}

// 添加拟合结果文本
void addFitText(double mean, double meanErr, double sigma, double sigmaErr, double x, double y, int color = kRed) {
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.045);
    latex->SetTextColor(color);
    
    TString meanText = TString::Format("Mean = %.6f #pm %.6f", mean, meanErr);
    TString sigmaText = TString::Format("Sigma = %.6f #pm %.6f", sigma, sigmaErr);
    
    latex->DrawLatex(x, y, meanText.Data());
    latex->DrawLatex(x, y-0.07, sigmaText.Data());
}

// 处理单个直方图
void processHistogram(const char* filename, const char* histname, const char* outputname, 
                     const char* ytitle, const char* xtitle,
                     double xmin, double xmax, double fitmin, double fitmax, int rebin = 1) {
    // 打开文件获取直方图
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    TH1* hist = (TH1*)file->Get(histname);
    if (!hist) {
        std::cerr << "Histogram not found: " << histname << std::endl;
        file->Close();
        return;
    }
    
    // 准备画布和直方图
    TCanvas* canvas = createCanvas();
    TH1* h = (TH1*)hist->Clone("h_clone");
    
    // 处理直方图
    if (rebin > 1) h->Rebin(rebin);
    h->Scale(1.0/h->Integral());
    setHistStyle(h, xtitle, ytitle, xmin, xmax);
    
    // 绘制直方图并拟合
    h->Draw("E");
    TF1* fitFunc = new TF1("fitFunc", "gaus", fitmin, fitmax);
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    h->Fit(fitFunc, "R");
    
    // 添加拟合文本
    addFitText(fitFunc->GetParameter(1), fitFunc->GetParError(1), 
               fitFunc->GetParameter(2), fitFunc->GetParError(2), 0.6, 0.85);
    
    // 保存结果并清理
    canvas->SaveAs(outputname);
    delete canvas;
    delete h;
    file->Close();
}

// 比较数据和MC直方图
void compareDataMC(const char* dataFile, const char* mcFile, const char* dataHist, const char* mcHist,
                  const char* outputname, const char* xtitle, 
                  double xmin, double xmax, double fitmin, double fitmax, int rebin = 1) {
    // 打开文件获取直方图
    TFile* fData = TFile::Open(dataFile);
    TFile* fMC = TFile::Open(mcFile);
    
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }
    
    TH1* hData = (TH1*)fData->Get(dataHist);
    TH1* hMC = (TH1*)fMC->Get(mcHist);
    
    if (!hData || !hMC) {
        std::cerr << "Histograms not found" << std::endl;
        fData->Close();
        fMC->Close();
        return;
    }
    
    // 准备画布和直方图
    TCanvas* canvas = createCanvas();
    TH1* hDataClone = (TH1*)hData->Clone("hData_clone");
    TH1* hMCClone = (TH1*)hMC->Clone("hMC_clone");
    
    // 处理直方图
    if (rebin > 1) {
        hDataClone->Rebin(rebin);
        hMCClone->Rebin(rebin);
    }
    hDataClone->Scale(1.0/hDataClone->Integral());
    hMCClone->Scale(1.0/hMCClone->Integral());
    
    setHistStyle(hDataClone, xtitle, "Normalized Events", xmin, xmax, kBlack, 20);
    setHistStyle(hMCClone, xtitle, "Normalized Events", xmin, xmax, kBlue, 24);
    
    // 绘制直方图并拟合
    hDataClone->Draw("E");
    hMCClone->Draw("E SAME");
    
    TF1* fitData = new TF1("fitData", "gaus", fitmin, fitmax);
    TF1* fitMC = new TF1("fitMC", "gaus", fitmin, fitmax);
    
    fitData->SetLineColor(kRed);
    fitData->SetLineWidth(2);
    fitMC->SetLineColor(kGreen+2);
    fitMC->SetLineWidth(2);
    
    hDataClone->Fit(fitData, "R");
    hMCClone->Fit(fitMC, "R");
    
    // 添加拟合文本
    addFitText(fitData->GetParameter(1), fitData->GetParError(1), 
              fitData->GetParameter(2), fitData->GetParError(2), 0.2, 0.85, kRed);
    
    addFitText(fitMC->GetParameter(1), fitMC->GetParError(1), 
              fitMC->GetParameter(2), fitMC->GetParError(2), 0.2, 0.62, kGreen+2);
    
    // 添加图例
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hDataClone, "DATA", "lep");
    legend->AddEntry(hMCClone, "MC", "lep");
    legend->Draw();
    
    // 保存结果并清理
    canvas->SaveAs(outputname);
    delete legend;
    delete canvas;
    delete hDataClone;
    delete hMCClone;
    fData->Close();
    fMC->Close();
}

// 处理能量相关的直方图
void processEkHistograms(const char* dataFile, const char* mcFile, const char* dataHist, const char* mcHist,
                        const char* outputPDF, const char* diffpng) {
    // 打开文件获取直方图
    TFile* fData = TFile::Open(dataFile);
    TFile* fMC = TFile::Open(mcFile);
    
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }
    
    TH2* hData2D = (TH2*)fData->Get(dataHist);
    TH2* hMC2D = (TH2*)fMC->Get(mcHist);
    
    if (!hData2D || !hMC2D) {
        std::cerr << "Histograms not found" << std::endl;
        fData->Close();
        fMC->Close();
        return;
    }
    
    // 创建点图和画布
    TGraph* diffGraph = new TGraph();
    TCanvas* canvas = createCanvas();
    canvas->Print(Form("%s[", outputPDF)); // 打开PDF
    
    // 处理每个Ek/n bin
    int nBinsY = hData2D->GetNbinsY();
    
    for (int i = 1; i <= nBinsY; i++) {
        double ek = hData2D->GetYaxis()->GetBinCenter(i);
        if (ek < 2 || ek > 25) continue;
        
        // 获取投影直方图
        TH1D* hDataProj = hData2D->ProjectionX(Form("hDataProj_%d", i), i, i);
        TH1D* hMCProj = hMC2D->ProjectionX(Form("hMCProj_%d", i), i, i);
        
        // 检查有效性
        if (hDataProj->Integral() < 10 || hMCProj->Integral() < 10) {
            delete hDataProj;
            delete hMCProj;
            continue;
        }
        
        // 处理直方图
        hDataProj->Rebin(4);
        hMCProj->Rebin(4);
        hDataProj->Scale(1.0/hDataProj->Integral());
        hMCProj->Scale(1.0/hMCProj->Integral());
        
        // 设置样式
        setHistStyle(hDataProj, "#Delta(R/#beta)/|R/#beta|", "Normalized Events", 
                    hDataProj->GetXaxis()->GetXmin(), hDataProj->GetXaxis()->GetXmax(), kBlack, 20);
        setHistStyle(hMCProj, "#Delta(R/#beta)/|R/#beta|", "Normalized Events", 
                    hMCProj->GetXaxis()->GetXmin(), hMCProj->GetXaxis()->GetXmax(), kBlue, 24);
        
        hDataProj->SetTitle(Form("E_{k}/n = %.2f GeV", ek));
        
        // 自动确定拟合范围
        double xmin = hDataProj->GetMean() - 2 * hDataProj->GetRMS();
        double xmax = hDataProj->GetMean() + 2 * hDataProj->GetRMS();
        
        // 拟合直方图
        TF1* fitData = new TF1(Form("fitData_%d", i), "gaus", xmin, xmax);
        TF1* fitMC = new TF1(Form("fitMC_%d", i), "gaus", xmin, xmax);
        
        fitData->SetLineColor(kRed);
        fitData->SetLineWidth(2);
        fitMC->SetLineColor(kGreen+2);
        fitMC->SetLineWidth(2);
        
        hDataProj->Fit(fitData, "RQ");
        hMCProj->Fit(fitMC, "RQ");
        
        // 获取拟合结果
        double meanData = fitData->GetParameter(1);
        double meanErrData = fitData->GetParError(1);
        double sigmaData = fitData->GetParameter(2);
        double sigmaErrData = fitData->GetParError(2);
        
        double meanMC = fitMC->GetParameter(1);
        double meanErrMC = fitMC->GetParError(1);
        double sigmaMC = fitMC->GetParameter(2);
        double sigmaErrMC = fitMC->GetParError(2);
        
        // 添加到差异图
        diffGraph->SetPoint(diffGraph->GetN(), ek, meanMC - meanData);
        
        // 绘制并保存
        canvas->Clear();
        hDataProj->Draw("E");
        hMCProj->Draw("E SAME");
        
        // 添加图例
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hDataProj, "DATA", "lep");
        legend->AddEntry(hMCProj, "MC", "lep");
        legend->Draw();
        
        // 添加拟合文本
        addFitText(meanData, meanErrData, sigmaData, sigmaErrData, 0.2, 0.85, kRed);
        addFitText(meanMC, meanErrMC, sigmaMC, sigmaErrMC, 0.2, 0.62, kGreen+2);
        
        canvas->Print(outputPDF);
        
        // 清理
        delete legend;
        delete hDataProj;
        delete hMCProj;
    }
    
    // 结束PDF
    canvas->Print(Form("%s]", outputPDF));
    
    // 创建差异图
    canvas->Clear();
    
    // 设置差异图样式
    diffGraph->SetTitle("Mean Difference vs E_{k}/n");
    //diffGraph->GetXaxis()->SetTitle("E_{k}/n [GeV]");
    //diffGraph->GetYaxis()->SetTitle("#mu_{MC} - #mu_{DATA}");
    diffGraph->GetXaxis()->SetRangeUser(2, 25);
    diffGraph->GetYaxis()->SetRangeUser(-0.07e-3, 0.07e-3);
    diffGraph->SetMarkerStyle(20);
    diffGraph->SetMarkerSize(1.2);
    diffGraph->SetMarkerColor(kBlue);
    
    diffGraph->Draw("AP");
    canvas->SaveAs(diffpng);
    cout<<11111<<endl;
    
    // 清理
    delete diffGraph;
    delete canvas;
    fData->Close();
    fMC->Close();
}

void richAGLCali() {
    // 基本设置
    const char* outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/richBeta";
    gSystem->mkdir(outputDir, kTRUE);
    setStyle();
    
    // 定义需要处理的直方图列表
    struct HistConfig {
        const char* filename;
        const char* histname;
        const char* outputname;
        const char* xtitle;
        double xmin, xmax, fitmin, fitmax;
    };
    
    // 单直方图处理配置
    vector<HistConfig> configs = {
        {"/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7_rich.root", "BetaCheck_Beryllium", 
         Form("%s/BetaCheck_Beryllium.png", outputDir), "#Delta#beta/|#beta|", 
         -0.0016, 0.0016, -0.0008, 0.0008},
        {"/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7_rich.root", "RigCheck_Beryllium", 
         Form("%s/RigCheck_Beryllium.png", outputDir), "#DeltaR/R", 
         -0.005, 0.008, -0.0006, 0.0016},
        {"/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7_rich.root", "RigBetaCheck_Beryllium", 
         Form("%s/RigBetaCheck_Beryllium.png", outputDir), "#Delta(R/#beta)/|R/#beta|", 
         -0.002, 0.002, -0.0008, 0.0008},
        {"/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Oxy16_rich.root", "RigBetaCheck_Oxygen", 
         Form("%s/RigBetaCheck_Oxygen.png", outputDir), "#Delta(R/#beta)/|R/#beta|", 
         -0.0013, 0.0013, -0.00055, 0.00062}
    };
    
    // 处理单直方图
    for (const auto& config : configs) {
        processHistogram(config.filename, config.histname, config.outputname, 
                         "Normalized Events", config.xtitle, 
                         config.xmin, config.xmax, config.fitmin, config.fitmax, 4);
    }
    
    // 处理DATA/MC比较
    compareDataMC(
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Oxy16_rich.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/O16_rich.root",
        "RigBetaCheck_Oxygen", "RigBetaCheck_Oxygen",
        Form("%s/RigBetaCheck_Oxygen_DataMC.png", outputDir),
        "#Delta(R/#beta)/|R/#beta|", -0.0013, 0.0013, -0.00055, 0.00062, 4
    );
    
    // 处理Ek相关的直方图
    processEkHistograms(
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Oxy16_rich.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/O16_rich.root",
        "RigBetaCheck_Ek_Oxygen", "RigBetaCheck_Ek_Oxygen",
        Form("%s/RigBetaCheck_Ek_Projections.pdf", outputDir),
        Form("%s/RigBetaCheck_Ek_Difference.png", outputDir)
    );
}