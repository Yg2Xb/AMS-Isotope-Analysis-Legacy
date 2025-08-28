// 合并权重结构体
struct MCWeight { double wB11=1, wB10=1; };

// 权重函数
MCWeight getMCWeight(double B11Frac, const char* fB11, const char* fB10) {
    TFile *fileB11 = TFile::Open(fB11), *fileB10 = TFile::Open(fB10);
    if (!fileB11 || !fileB10) throw std::runtime_error("open MC root failed!");
    TH1D *hMCnumB11 = nullptr, *hMCnumB10 = nullptr;
    fileB11->GetObject("hMCnum", hMCnumB11);
    fileB10->GetObject("hMCnum", hMCnumB10);
    double nB11 = hMCnumB11->GetBinContent(1), nB10 = hMCnumB10->GetBinContent(1), s = nB11 + nB10;
    if(nB11<=0||nB10<=0) throw std::runtime_error("hMCnum bin content <=0");
    fileB11->Close(); fileB10->Close();

    return {(s/nB11)*B11Frac, (s/nB10)*(1-B11Frac)}; // B11权重=0.7
}

// 自动选择探测器
TString getDetectorByEnergy(double energy) {
    if (energy < 1.17) return "TOF";
    else if (energy < 3.23) return "NaF";
    else return "AGL";
}

void CpFragBeIso() {
    gStyle->SetOptStat(0);

    TFile *B11IsoFracFile = TFile::Open("/eos/ams/user/z/zetong/Boron_isotope/pre/B11f1.root"); 
    TH1D* h_B11_Frac = (TH1D*)B11IsoFracFile->Get("B11");

    // 文件路径
    TString issFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/ISS_FragmentBer_MassFitAnalysis_unb_perfectQfit.root";
    TString mcB11File = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_narrow_frag_unbL1_SDIATBeBcut.root";
    TString mcB10File = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_temp_narrow_frag_unbL1_SDIATBeBcut.root";
    
    std::vector<TString> beNames = {"Be7", "Be9", "Be10"};
    std::vector<int> beNumbers = {7, 9, 10};

    // --- 1. 读取ISS三种Be的Ratio直方图
    TFile* fISS = TFile::Open(issFile);
    if(!fISS || fISS->IsZombie()) { std::cerr<<"Failed to open ISS file"<<std::endl; return;}
    std::vector<TH1F*> hISS;
    for (const auto& be : beNames) {
        TString hname = "BeIso_BorSource_Ratio_" + be;
        TH1F* h = nullptr;
        fISS->GetObject(hname, h);
        if(!h) { std::cerr<<"No ISS hist: "<<hname<<std::endl; hISS.push_back(nullptr); continue; }
        h = (TH1F*)h->Clone("iss_"+be);
        h->SetDirectory(0);
        hISS.push_back(h);
    }

    TFile* fB11 = TFile::Open(mcB11File), *fB10 = TFile::Open(mcB10File);
    if (!fB11 || !fB10) { std::cerr<<"Failed to open MC files"<<std::endl; return;}

    // --- 3. 合并MC（按能量区间自动选探测器）
    std::vector<TH1F*> hMC_Be;
    for (size_t ibe=0; ibe<beNames.size(); ++ibe) {
        cout<<ibe<<endl;
        if(!hISS[ibe]) { hMC_Be.push_back(nullptr); continue; }

        // 为每个探测器缓存已处理的直方图
        std::map<TString, TH1F*> cachedHistsB11, cachedHistsB10;
        std::map<TString, TH1F*> cachedL1HistsB11, cachedL1HistsB10;

        TH1F* hCombined = (TH1F*)hISS[ibe]->Clone("hCombined");
        hCombined->Reset();

        // 合并L1 B
        TH1F* hL1BCombined = (TH1F*)hISS[ibe]->Clone("hL1BCombined");
        hL1BCombined->Reset();

        for (int bin=1; bin<=hCombined->GetNbinsX(); ++bin) {
            cout<<bin<<endl;
            double B11Frac = h_B11_Frac->GetBinContent(bin);
            //if(B11Frac <= 0) B11Frac = 0.7;
            B11Frac = 0.7;
            MCWeight w = getMCWeight(B11Frac,
                "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_bkg.root",
                "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_bkg.root"
            );
            double energy = hCombined->GetBinCenter(bin);
            TString det = getDetectorByEnergy(energy);

            TString hnameB11 = TString::Format("SourceToIsoCounts_%s_%s", det.Data(), beNames[ibe].Data());
            TString hnameB10 = TString::Format("SourceToIsoCounts_%s_%s", det.Data(), beNames[ibe].Data());
            
            // 检查是否已经处理过这个探测器的直方图
            TH1F *h_BeFromB11 = nullptr, *h_BeFromB10 = nullptr;
            
            if(cachedHistsB11.find(det) == cachedHistsB11.end()) {
                // 第一次遇到这个探测器，获取并处理直方图
                fB11->GetObject(hnameB11, h_BeFromB11);
                fB10->GetObject(hnameB10, h_BeFromB10);
                if(h_BeFromB11 && h_BeFromB10) {
                    h_BeFromB11 = (TH1F*)h_BeFromB11->Clone(Form("cloned_B11_%s_%s", det.Data(), beNames[ibe].Data()));
                    h_BeFromB10 = (TH1F*)h_BeFromB10->Clone(Form("cloned_B10_%s_%s", det.Data(), beNames[ibe].Data()));
                    h_BeFromB11->Rebin(2);
                    h_BeFromB10->Rebin(2);
                    cout<<"Iso rebin for detector: "<<det<<endl;
                    
                    // 缓存处理后的直方图
                    cachedHistsB11[det] = h_BeFromB11;
                    cachedHistsB10[det] = h_BeFromB10;
                }
            } else {
                // 使用已经处理过的直方图
                h_BeFromB11 = cachedHistsB11[det];
                h_BeFromB10 = cachedHistsB10[det];
            }

            if(h_BeFromB11 && h_BeFromB10) {
                double valB11 = h_BeFromB11->GetBinContent(h_BeFromB11->FindBin(energy));
                double valB10 = h_BeFromB10->GetBinContent(h_BeFromB10->FindBin(energy));
                double errB11 = h_BeFromB11->GetBinError(h_BeFromB11->FindBin(energy));
                double errB10 = h_BeFromB10->GetBinError(h_BeFromB10->FindBin(energy));
                double valBe = w.wB11*valB11 + w.wB10*valB10;
                double errBe = sqrt(pow(w.wB11*errB11,2) + pow(w.wB10*errB10,2));
                hCombined->SetBinContent(bin, valBe);
                hCombined->SetBinError(bin, errBe);
            }

            TString hL1B11name = TString::Format("TruthL1SourceCounts_%s", det.Data());
            TString hL1B10name = TString::Format("TruthL1SourceCounts_%s", det.Data());
            
            TH1F *hL1B11 = nullptr, *hL1B10 = nullptr;
            
            if(cachedL1HistsB11.find(det) == cachedL1HistsB11.end()) {
                // 第一次遇到这个探测器，获取并处理直方图
                fB11->GetObject(hL1B11name, hL1B11);
                fB10->GetObject(hL1B10name, hL1B10);
                if(hL1B11 && hL1B10) {
                    hL1B11 = (TH1F*)hL1B11->Clone(Form("cloned_L1B11_%s", det.Data()));
                    hL1B10 = (TH1F*)hL1B10->Clone(Form("cloned_L1B10_%s", det.Data()));
                    hL1B11->Rebin(2);
                    hL1B10->Rebin(2);
                    cout<<"L1 rebin for detector: "<<det<<endl;
                    
                    // 缓存处理后的直方图
                    cachedL1HistsB11[det] = hL1B11;
                    cachedL1HistsB10[det] = hL1B10;
                }
            } else {
                // 使用已经处理过的直方图
                hL1B11 = cachedL1HistsB11[det];
                hL1B10 = cachedL1HistsB10[det];
            }

            cout<<energy<<" "<<det<<" ";
            if(hL1B11 && hL1B10) {
                double vB11 = hL1B11->GetBinContent(hL1B11->FindBin(energy));
                double vB10 = hL1B10->GetBinContent(hL1B10->FindBin(energy));
                double eB11 = hL1B11->GetBinError(hL1B11->FindBin(energy));
                double eB10 = hL1B10->GetBinError(hL1B10->FindBin(energy));
                double vL1B = w.wB11*vB11 + w.wB10*vB10;
                double eL1B = sqrt(pow(w.wB11*eB11,2) + pow(w.wB10*eB10,2));
                hL1BCombined->SetBinContent(bin, vL1B);
                hL1BCombined->SetBinError(bin, eL1B);
                cout<<hL1B11->FindBin(energy)<<" "<<hL1B11->GetBinLowEdge(hL1B11->FindBin(energy))<<" "<<vB11<<endl;
            }
        }

        hCombined->Divide(hL1BCombined);
        hCombined->SetLineColor(kRed);
        hCombined->SetMarkerColor(kRed);
        hCombined->SetMarkerStyle(20);
        hCombined->SetLineWidth(2);
        hMC_Be.push_back(hCombined);

        delete hL1BCombined;

        // 清理缓存的直方图
        for(auto& pair : cachedHistsB11) delete pair.second;
        for(auto& pair : cachedHistsB10) delete pair.second;
        for(auto& pair : cachedL1HistsB11) delete pair.second;
        for(auto& pair : cachedL1HistsB10) delete pair.second;
    }

    // --- 4. 画图
    for(int i=0; i<3; ++i) {
        if(!hISS[i]) continue;
        // 自动y轴上下限
        double minY=1e9, maxY=-1e9;
        for(int bin=1; bin<=hISS[i]->GetNbinsX(); ++bin) {
            double y1 = hISS[i]->GetBinContent(bin), y2 = hMC_Be[i] ? hMC_Be[i]->GetBinContent(bin) : 0;
            if(y1>0) { minY=std::min(minY,y1); maxY=std::max(maxY,y1);}
            if(y2>0) { minY=std::min(minY,y2); maxY=std::max(maxY,y2);}
        }
        double dy = (maxY-minY)*0.2;
        double ylow = std::max(0., minY-dy), yhigh = maxY+dy;

        TCanvas* c = new TCanvas(TString::Format("c_Be%d", beNumbers[i]), TString::Format("L2 Fragment ^{%d}Be / L1 ^{11}B", beNumbers[i]), 900, 700);
        hISS[i]->SetTitle(TString::Format("L2 Fragment ^{%d}Be / L1 Boron", beNumbers[i]));
        hISS[i]->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
        hISS[i]->GetYaxis()->SetTitle(TString::Format("^{%d}Be / Boron Ratio", beNumbers[i]));
        hISS[i]->GetYaxis()->SetTitleOffset(1.62);
        hISS[i]->GetXaxis()->SetRangeUser(0.0, 16.3);
        hISS[i]->GetYaxis()->SetRangeUser(ylow, yhigh);
        hISS[i]->SetTitleSize(0.035,"t"); // title字小一点
        hISS[i]->Draw("pz");

        TLegend* legend = new TLegend(0.65,0.68,0.88,0.88); // legend往下
        legend->SetBorderSize(0);
        legend->SetFillColorAlpha(0, 0.6);
        legend->SetTextSize(0.035);
        legend->AddEntry(hISS[i], "ISS Data", "ep");

        if(hMC_Be[i]) {
            hMC_Be[i]->Draw("pz SAME");
            legend->AddEntry(hMC_Be[i], "B10 B11 MC Mix", "ep");
        }

        legend->Draw();
        c->SetGrid();
        c->SaveAs(TString::Format("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be%d_ISS_MC_unb_perfectQfit.png", beNumbers[i]));
        if(i==2)
        {
            hISS[0]->Draw("pz");
            hISS[0]->SetTitle("L2 Fragment ^{7}Be ^{9}Be / L1 Boron;E_{k}/n [GeV/n];Ratio");
            hISS[0]->SetMarkerColor(kOrange+1); 
            hISS[0]->SetLineColor(kOrange+1); 
            hISS[1]->SetMarkerColor(kGreen+2); 
            hISS[1]->SetLineColor(kGreen+2); 
            hISS[0]->GetYaxis()->SetRangeUser(0.0, 0.006);
            hISS[1]->Draw("pz same");
            legend->Clear();
            legend->AddEntry(hISS[0], "ISS L2 Be7 / L1 B", "ep");
            legend->AddEntry(hISS[1], "ISS L2 Be9 / L1 B", "ep");
            legend->Draw();
            c->SetGrid();
            c->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/Be79_ISS_unb_perfectQfit.png");
        }
        delete c;
    }
    fISS->Close(); fB11->Close(); fB10->Close();
}