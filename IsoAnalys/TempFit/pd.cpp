void pd() {
    TFile* file = TFile::Open("root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber_temp_pure10.root");
    if (!file) {
        printf("Cannot open file\n");
        return;
    }

    // 设置画布样式
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    c->SetTickx();
    c->SetTicky();

    // Mass10的直方图名称
    const char* histNames[] = {
        "h_inv_ToFMass_BerUseMass10",
        "h_inv_NaFMass_BerUseMass10",
        "h_inv_AGLMass_BerUseMass10",
        "h_inv_ToFEst_BerUseMass10",
        "h_inv_NaFEst_BerUseMass10",
        "h_inv_AGLEst_BerUseMass10"
    };

    // 能量范围
    struct Range {
        double min;
        double max;
        const char* detector;
    };
    
    Range ranges[] = {
        {0.4, 1.1, "ToF"},
        {1.1, 3.2, "NaF"},
        {3.0, 7.0, "AGL"}
    };

    for (int i = 0; i < 6; i++) {
        TH2D* h2 = (TH2D*)file->Get(histNames[i]);
        if (!h2) continue;

        // 确定这是质量直方图还是估计值直方图
        bool isEst = strstr(histNames[i], "Est") != nullptr;
        
        // 确定探测器类型
        const char* detector = "";
        if (strstr(histNames[i], "ToF")) detector = "ToF";
        else if (strstr(histNames[i], "NaF")) detector = "NaF";
        else if (strstr(histNames[i], "AGL")) detector = "AGL";

        // 找到对应的能量范围
        Range* range = nullptr;
        for (auto& r : ranges) {
            if (strcmp(r.detector, detector) == 0) {
                range = &r;
                break;
            }
        }
        if (!range) continue;

        // 找到对应的bin范围
        int binLow = h2->GetYaxis()->FindBin(range->min);
        int binHigh = h2->GetYaxis()->FindBin(range->max);

        // 投影到x轴
        TH1D* proj = h2->ProjectionX(Form("%s_proj", histNames[i]), binLow, binHigh);
        proj->Rebin(i<=2?2:5);
        // 设置x轴范围
        if (isEst) {
            proj->GetXaxis()->SetRangeUser(0, 1);
        } else {
            proj->GetXaxis()->SetRangeUser(0.05, 0.225);
        }

        // 设置标题和标签
        proj->SetTitle(Form("%s Projection (%.1f-%.1f GeV)", detector, range->min, range->max));
        proj->GetXaxis()->SetTitle(isEst ? "Likelihood Ratio" : "1/Mass [1/GeV]");
        proj->GetYaxis()->SetTitle("Entries");

        // 绘制和保存
        proj->Draw("hist");
        c->SaveAs(Form("./%s_projection.png", histNames[i]));
    }

    file->Close();
    delete c;
}