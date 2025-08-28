void draw_RTIcut() {
    // 输出目录
    TString outdir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/";
    // 打开root文件
    TFile *f = TFile::Open("/cvmfs/ams.cern.ch/Offline/AMSDataDir/v6.00/RTIcut.root");
    if (!f || f->IsZombie()) {
        printf("File not found!\n");
        return;
    }

    // hist名字和输出文件名对应表
    const char* histnames[] = {"hRcut", "hRcutMin", "hRcutAcc", "hRcutAccE"};
    const int nhist = 4;

    // 画图风格
    gStyle->SetOptStat(0);

    for (int i = 0; i < nhist; ++i) {
        // 取出hist
        TH2F* hist = (TH2F*) f->Get(histnames[i]);
        if (!hist) {
            printf("Histogram %s not found!\n", histnames[i]);
            continue;
        }

        // -- 1. 绘制二维分布 --
        TCanvas* c2d = new TCanvas(Form("c2d_%s",histnames[i]), histnames[i], 800, 600);
        c2d->SetRightMargin(0.15);  // 右边距大一点
        hist->GetXaxis()->SetRangeUser(0.1, 30);
        hist->GetYaxis()->SetRangeUser(0.1, 30);
        hist->GetXaxis()->SetTitleSize(0.06);
        hist->GetYaxis()->SetTitleSize(0.06);
        hist->GetXaxis()->SetLabelSize(0.06);
        hist->GetYaxis()->SetLabelSize(0.06);
        hist->Draw("colz");
        TString outname2d = outdir + Form("%s.png", histnames[i]);
        c2d->SaveAs(outname2d);
        delete c2d;
        printf("Saved %s\n", outname2d.Data());

        // -- 2. X/Y投影并叠加 --
        TString projtitle = Form("%s Projections", histnames[i]);
        TCanvas* cproj = new TCanvas(Form("cproj_%s",histnames[i]), projtitle, 800, 600);
        cproj->SetRightMargin(0.08);

        TH1D* hx = hist->ProjectionX(Form("%s_projX", histnames[i]));
        TH1D* hy = hist->ProjectionY(Form("%s_projY", histnames[i]));

        // 设置样式
        hx->SetLineColor(kRed);
        hx->SetLineWidth(2);
        hx->SetTitle(projtitle);
        hx->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
        hx->GetYaxis()->SetTitle("Counts");
        hx->GetXaxis()->SetRangeUser(0.1, 30);
        hx->GetXaxis()->SetTitleSize(0.06);
        hx->GetYaxis()->SetTitleSize(0.06);
        hx->GetXaxis()->SetLabelSize(0.06);
        hx->GetYaxis()->SetLabelSize(0.06);

        hy->SetLineColor(kBlack);
        hy->SetLineWidth(2);

        // 归一化到最大值为1方便对比（如不需要可注释掉）
        //hx->Scale(1.0/hx->GetMaximum());
        //hy->Scale(1.0/hy->GetMaximum());

        // 设置最大值以便两条曲线都显示完整
        double maxy = std::max(hx->GetMaximum(), hy->GetMaximum());
        hx->SetMaximum(maxy*1.15);

        hx->Draw("hist");
        hy->Draw("hist same");

        // 图例
        TLegend *leg = new TLegend(0.65,0.78,0.93,0.93);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hx, "Projection X", "l");
        leg->AddEntry(hy, "Projection Y", "l");
        leg->Draw();

        TString outnameproj = outdir + Form("%s_proj.png", histnames[i]);
        cproj->SaveAs(outnameproj);
        delete cproj; delete hx; delete hy; delete leg;
        printf("Saved %s\n", outnameproj.Data());
    }
    f->Close();
}