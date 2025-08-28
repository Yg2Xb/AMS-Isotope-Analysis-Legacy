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

    // 设置画布属性
    c->SetLogx(1);
    c->SetLogy(1);
    
    // 定义数组
    const char* dets[] = {"TOF", "NaF", "AGL"};
    const char* bes[] = {"Be7", "Be9", "Be10"};
    const char* types[] = {"NumEk", "genEk", "recEk", "genEk_reweighted", "recEk_reweighted", "NumEk_reweighted"};
    const char* acc_types[] = {"gen_", "rec_", "gen_reweighted", "rec_reweighted"};

    // 绘制事例数直方图
    for (const char* det : dets) {
        for (const char* be : bes) {
            // 创建多重图
            TCanvas *mc = new TCanvas(Form("mc_%s_%s", det, be), "", 1600, 1200);
            mc->Divide(3, 2);
            
            for (int i = 0; i < 6; ++i) {
                mc->cd(i + 1);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                
                TString histName = Form("h_%s_%s_%s", types[i], det, be);
                TH1F* h = (TH1F*)f->Get(histName);
                if (!h) continue;

                h->SetLineColor(kBlue);
                h->SetLineWidth(2);
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.8);
                h->Draw("hist");
            }
            mc->Print(pdf_path.c_str());
            delete mc;
        }
    }

// 绘制acceptance直方图
for (const char* det : dets) {
    for (const char* be : bes) {
        TCanvas *ac = new TCanvas(Form("ac_%s_%s", det, be), "", 1600, 1200);
        ac->Divide(2, 2);
        
        for (int i = 0; i < 4; ++i) {
            ac->cd(i + 1);
            gPad->SetLogx(1);
            gPad->SetLogy(1);
            
            TString histName = Form("h_acc_%s_%s_%s", acc_types[i], det, be);
            TH1F* h = (TH1F*)f->Get(histName);
            if (!h) continue;

            // 修改标题，添加探测器和同位素信息
            TString newTitle = Form("%s Acceptance (%s detector, %s);%s;Acceptance (m^2 sr)", 
                acc_types[i], det, be, h->GetXaxis()->GetTitle());
            //h->SetTitle(newTitle);

            h->SetLineColor(kRed);
            h->SetLineWidth(2);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.8);
            h->Draw("hist");
        }
        ac->Print(pdf_path.c_str());
        delete ac;
    }
}

    // 绘制源粒子通量拟合
    TF1* flux = (TF1*)f->Get("f_SourceFlux");
    if (flux) {
        TCanvas *fc = new TCanvas("fc", "Flux Fit", 800, 600);
        fc->SetLogx(1);
        fc->SetLogy(1);
        
        flux->SetTitle(Form("%s Flux Spline Fit Func;Rigidity;Flux", particle_name.c_str()));
        flux->SetLineColor(kBlue);
        flux->SetLineWidth(2);
        
        TFile *fluxFile = TFile::Open(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_%sFlux.root",
            particle_name.c_str()));
        if (fluxFile) {
            TGraphAsymmErrors* graph = (TGraphAsymmErrors*)fluxFile->Get("graph1");
            if (graph) {
                flux->Draw();
                graph->Draw("P SAME");
                
                // 添加图例
                TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                leg->AddEntry(flux, "Fit", "l");
                leg->AddEntry(graph, "Data", "p");
                leg->Draw();
            }
            fluxFile->Close();
        }
        fc->Print(pdf_path.c_str());
        delete fc;
    }

    // 关闭PDF文件
    c->Print((pdf_path + "]").c_str());

    delete c;
    f->Close();
}