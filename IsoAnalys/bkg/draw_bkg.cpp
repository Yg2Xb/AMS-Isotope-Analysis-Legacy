void draw_bkg(const string& source_particle = "Boron11") {
    // 解析粒子名称和质量数
    string particle_name;
    int mass_number;
    for(size_t i = 0; i < source_particle.length(); ++i) {
        if(isdigit(source_particle[i])) {
            particle_name = source_particle.substr(0, i);
            mass_number = stoi(source_particle.substr(i));
            break;
        }
    }

    // 构建文件路径
    string file_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/" + 
                      particle_name + std::to_string(mass_number) + "_bkg_analysis.root";
    string pdf_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/" + 
                     particle_name + std::to_string(mass_number) + "_bkg_analysis.pdf";

    // 打开输入文件
    TFile *f = TFile::Open(file_path.c_str());
    if (!f) {
        cout << "Cannot open file: " << file_path << endl;
        return;
    }

    // 创建画布和PDF
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Print((pdf_path + "[").c_str()); // 打开PDF文件

    // 定义常量数组
    const char* rwTags[] = {"orig", "reweight"};
    const char* selTypes[] = {"Complete", "Loose"};
    const char* bkgTypes[] = {"Bkg", "NoBkg", "StrictBkg"};
    const char* detNames[] = {"TOF", "NaF", "AGL"};
    const char* beIsotopes[] = {"Be7", "Be9", "Be10"};

    // ------------- 绘制生成数直方图 (h_num_) -------------
    {
        TCanvas *numCanvas = new TCanvas("num_canvas", "Generated Numbers", 1200, 800);
        numCanvas->Divide(3, 2);
        
        int pad = 1;
        for(int iso = 0; iso < 3; ++iso) {
            for(int rw = 0; rw < 2; ++rw) {
                numCanvas->cd(pad++);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                
                // 构建直方图名称
                TString histName = Form("h_num_%s_%s", rwTags[rw], beIsotopes[iso]);
                TH1F* h = (TH1F*)f->Get(histName);
                if(!h) {
                    cout << "Histogram not found: " << histName << endl;
                    continue;
                }
                
                h->SetLineColor(kBlue);
                h->SetLineWidth(2);
                h->SetTitle(Form("%s %s Generated Events", rwTags[rw], beIsotopes[iso]));
                h->Draw("hist");
            }
        }
        numCanvas->Print(pdf_path.c_str());
        delete numCanvas;
    }

    // ------------- 绘制背景研究直方图 (h_bkg_study_) -------------
    // 对每个同位素单独绘制
    for(int iso = 0; iso < 3; ++iso) {
        const char* isotopeName = beIsotopes[iso];
        
        // 对每个探测器类型
        for(int det = 0; det < 3; ++det) {
            const char* detName = detNames[det];
            
            // 对原始/重新加权分别绘制
            for(int rw = 0; rw < 2; ++rw) {
                const char* rwTag = rwTags[rw];
                
                // 创建一个画布，显示所有sel和bkg组合
                TCanvas *bkgCanvas = new TCanvas(Form("bkg_%s_%s_%s", isotopeName, detName, rwTag), 
                                                "Background Study", 1200, 800);
                bkgCanvas->Divide(3, 2);
                
                int pad = 1;
                for(int sel = 0; sel < 2; ++sel) {
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        bkgCanvas->cd(pad++);
                        gPad->SetLogx(1);
                        gPad->SetLogy(1);
                        
                        // 构建直方图名称
                        TString histName = Form("h_bkg_%s_%s_%s_%s_%s", 
                                               rwTag, selTypes[sel], bkgTypes[bkg], detName, isotopeName);
                        TH1F* h = (TH1F*)f->Get(histName);
                        if(!h) {
                            cout << "Histogram not found: " << histName << endl;
                            continue;
                        }
                        
                        h->SetLineColor(kBlue);
                        h->SetLineWidth(2);
                        h->SetTitle(Form("%s %s %s %s %s", 
                                        rwTag, selTypes[sel], bkgTypes[bkg], detName, isotopeName));
                        h->Draw("hist");
                    }
                }
                bkgCanvas->Print(pdf_path.c_str());
                delete bkgCanvas;
            }
        }
    }

    // ------------- 绘制接收度直方图 (acceptance_hists_) -------------
    // 对每个同位素单独绘制
    for(int iso = 0; iso < 3; ++iso) {
        const char* isotopeName = beIsotopes[iso];
        
        // 对每个探测器类型
        for(int det = 0; det < 3; ++det) {
            const char* detName = detNames[det];
            
            // 对原始/重新加权分别绘制
            for(int rw = 0; rw < 2; ++rw) {
                const char* rwTag = rwTags[rw];
                
                // 创建一个画布，显示所有sel和bkg组合
                TCanvas *accCanvas = new TCanvas(Form("acc_%s_%s_%s", isotopeName, detName, rwTag), 
                                               "Acceptance", 1200, 800);
                accCanvas->Divide(3, 2);
                
                int pad = 1;
                for(int sel = 0; sel < 2; ++sel) {
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        accCanvas->cd(pad++);
                        gPad->SetLogx(1);
                        // 不对接收度使用对数Y轴
                        gPad->SetLogy(0);
                        
                        // 构建直方图名称
                        TString histName = Form("h_acc_%s_%s_%s_%s_%s", 
                                              rwTag, selTypes[sel], bkgTypes[bkg], detName, isotopeName);
                        TH1F* h = (TH1F*)f->Get(histName);
                        if(!h) {
                            cout << "Histogram not found: " << histName << endl;
                            continue;
                        }
                        
                        h->SetLineColor(kRed);
                        h->SetLineWidth(2);
                        h->SetMarkerStyle(20);
                        h->SetMarkerColor(kRed);
                        h->SetMarkerSize(0.8);
                        h->Draw("hist p");
                        
                        // 设置Y轴范围，避免太小的值显示不清
                        double maxY = h->GetMaximum();
                        h->GetYaxis()->SetRangeUser(0, maxY * 1.2);
                    }
                }
                accCanvas->Print(pdf_path.c_str());
                delete accCanvas;
            }
        }
    }

    // ------------- 比较不同背景选择对接收度的影响 -------------
    // 对每个同位素单独绘制
    for(int iso = 0; iso < 3; ++iso) {
        const char* isotopeName = beIsotopes[iso];
        
        // 对每个探测器类型
        for(int det = 0; det < 3; ++det) {
            const char* detName = detNames[det];
            
            // 创建一个画布，比较不同bkg类型
            TCanvas *compCanvas = new TCanvas(Form("comp_%s_%s", isotopeName, detName), 
                                            "Background Comparison", 1200, 800);
            compCanvas->Divide(2, 2);
            
            // 颜色数组，用于区分不同曲线
            int colors[] = {kRed, kBlue, kGreen+2};
            
            // 对不同sel和rw组合
            for(int sel = 0; sel < 2; ++sel) {
                for(int rw = 0; rw < 2; ++rw) {
                    compCanvas->cd(sel*2 + rw + 1);
                    gPad->SetLogx(1);
                    gPad->SetLogy(0);
                    
                    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                    double maxY = 0;
                    
                    // 叠加绘制不同bkg类型
                    for(int bkg = 0; bkg < 3; ++bkg) {
                        TString histName = Form("h_acc_%s_%s_%s_%s_%s", 
                                              rwTags[rw], selTypes[sel], bkgTypes[bkg], detName, isotopeName);
                        TH1F* h = (TH1F*)f->Get(histName);
                        if(!h) continue;
                        
                        h->SetLineColor(colors[bkg]);
                        h->SetMarkerColor(colors[bkg]);
                        h->SetLineWidth(2);
                        h->SetMarkerStyle(20);
                        h->SetMarkerSize(0.8);
                        
                        maxY = std::max(maxY, h->GetMaximum());
                        
                        if(bkg == 0) {
                            h->SetTitle(Form("%s %s %s Acceptance Comparison", 
                                           rwTags[rw], selTypes[sel], detName));
                            h->Draw("hist p");
                        } else {
                            h->Draw("hist p same");
                        }
                        
                        leg->AddEntry(h, bkgTypes[bkg], "lp");
                    }
                    
                    // 调整Y轴范围
                    TH1F* firstHist = (TH1F*)f->Get(Form("h_acc_%s_%s_%s_%s_%s", 
                                                       rwTags[rw], selTypes[sel], bkgTypes[0], detName, isotopeName));
                    if(firstHist) {
                        firstHist->GetYaxis()->SetRangeUser(0, maxY * 1.2);
                    }
                    
                    leg->Draw();
                }
            }
            compCanvas->Print(pdf_path.c_str());
            delete compCanvas;
        }
    }

    // ------------- 绘制源粒子通量拟合 -------------
    {
        TF1* flux = (TF1*)f->Get("f_SourceFlux");
        if(flux) {
            TCanvas *fc = new TCanvas("fc", "Flux Fit", 800, 600);
            fc->SetLogx(1);
            fc->SetLogy(1);
            
            flux->SetTitle(Form("%s Flux Spline Fit Func;Rigidity;Flux", particle_name.c_str()));
            flux->SetLineColor(kBlue);
            flux->SetLineWidth(2);
            
            TFile *fluxFile = TFile::Open(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_%sFlux.root",
                particle_name.c_str()));
            if(fluxFile) {
                TGraphAsymmErrors* graph = (TGraphAsymmErrors*)fluxFile->Get("graph1");
                if(graph) {
                    graph->SetMarkerStyle(20);
                    graph->SetMarkerColor(kRed);
                    graph->SetLineColor(kRed);
                    
                    flux->Draw();
                    graph->Draw("P SAME");
                    
                    // 添加图例
                    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                    leg->AddEntry(flux, "Fit", "l");
                    leg->AddEntry(graph, "Data", "p");
                    leg->Draw();
                }
                fluxFile->Close();
            } else {
                flux->Draw();
            }
            fc->Print(pdf_path.c_str());
            delete fc;
        }
    }

    // 关闭PDF文件
    c->Print((pdf_path + "]").c_str());

    delete c;
    f->Close();
}