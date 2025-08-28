#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void Drawtofcp() {
    // 打开 ROOT 文件
    TFile *fileBe10 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7.root");
    TFile *fileBer10 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root");

    // 检查文件是否成功打开
    if (!fileBe10 || fileBe10->IsZombie() || !fileBer10 || fileBer10->IsZombie()) {
        printf("Error opening one or both of the files!\n");
        return;
    }

    // 从文件中获取 TH2F 直方图
    TH2F *hBe10_2D = (TH2F*)fileBe10->Get("RBeta_Ek_TOF_NaF_Berlium");
    TH2F *hBer10_2D = (TH2F*)fileBer10->Get("RBeta_Ek_TOF_NaF_Berlium");

    // 检查是否成功获取直方图
    if (!hBe10_2D || !hBer10_2D) {
        printf("Error retrieving histograms from files!\n");
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "TOF Beta Comparison", 800, 600);
    c1->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/compareBeta7.pdf[");
    for(int i = 5; i < 20; i++){
        // 将 TH2F 直方图在 y 轴的第 2 个 bin 上进行投影
        TH1 *hBe10_proj = hBe10_2D->ProjectionX("hBe10_proj", i, i);
        TH1 *hBer10_proj = hBer10_2D->ProjectionX("hBer10_proj", i, i);

        hBer10_proj->Rebin(100);
        hBe10_proj->Rebin(100);
        hBer10_proj->GetXaxis()->SetRangeUser(-0.06,.06);
        hBer10_proj->SetLineColor(kBlack);
        hBer10_proj->SetMarkerColor(kBlack);
        hBer10_proj->SetMarkerStyle(20);
        hBer10_proj->SetMarkerSize(1.0);

        hBe10_proj->SetLineColor(kRed);
        hBe10_proj->SetMarkerColor(kRed);
        hBe10_proj->SetMarkerStyle(20);
        hBe10_proj->SetMarkerSize(1.0);

        // 设置 X 轴和 Y 轴的标签和标题样式
        hBer10_proj->SetTitle(";1/beta_{TOF}-1/beta_{NaF};Normalized Events");
        hBer10_proj->SetTitle(Form("[%.2f-%.2f] GeV/n", hBer10_2D->GetYaxis()->GetBinLowEdge(i), hBer10_2D->GetYaxis()->GetBinLowEdge(i+1)));
        hBer10_proj->GetXaxis()->SetTitleOffset(1.15);
        hBer10_proj->GetYaxis()->SetTitleOffset(1.15);
        hBer10_proj->GetXaxis()->SetLabelSize(0.06);  // 设置 x 轴标签大小
        hBer10_proj->GetXaxis()->SetTitleSize(0.06);  // 设置 x 轴标题大小
        hBer10_proj->GetYaxis()->SetTitleSize(0.06);  // 设置 y 轴标题大小
        hBer10_proj->GetYaxis()->SetLabelSize(0.06);  // 设置 y 轴标签大小

        // 绘制直方图
        hBer10_proj->DrawNormalized("P");  // "E P" 表示绘制误差条和点
        hBe10_proj->DrawNormalized("P SAME");  // "SAME HIST" 表示在同一张图上绘制线条直方图

        // 创建无边框、透明背景的图例
        TLegend *legend = new TLegend(0.65, 0.65, 0.88, 0.88);
        legend->SetBorderSize(0);  // 无边框
        legend->SetFillStyle(0);   // 透明背景
        legend->SetTextSize(0.04); // 文本大小

        // 在图例中添加条目
        legend->AddEntry(hBer10_proj, "Ber ISS", "p");
        legend->AddEntry(hBe10_proj, "Be7 MC", "p");
        legend->Draw();

        c1->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/compareBeta7.pdf");
        // 保存为 PNG 文件
    }
    c1->Print("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/tofBeta/compareBeta7.pdf]");

    // 关闭 ROOT 文件
    fileBe10->Close();
    fileBer10->Close();
}