#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TLine.h>
#include <TLegend.h>
#include "../Tool.h"
using namespace AMS_Iso;

void Beta_Cutoff() {
    // 设置文件路径
    const char* inputFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root";
    const char* outputDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/";

    // 打开文件
    TFile *f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) {
        printf("Error opening file: %s\n", inputFile);
        return;
    }

    // 获取树
    TTree *tree = (TTree*)f->Get("saveTree");
    if (!tree) {
        printf("Error getting tree from file\n");
        f->Close();
        return;
    }

    // 设置画布和直方图样式
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    // 创建画布
    TCanvas *c = new TCanvas("c", "Beta vs Cutoff", 800, 600);
    c->SetRightMargin(0.16);
    c->SetLogz();

    // NaF部分
    c->cd();
    tree->Draw("richBeta:cutOffRig>>hNaF(250,0,30,250,0.75,1.005)",
               "(cutStatus & 0x04) == 0x04 && rich_NaF", "colz");
    TH2D *hNaF = (TH2D*)gDirectory->Get("hNaF");
    if (hNaF) {
        hNaF->SetTitle("Beta vs Cutoff (NaF);Cutoff Rigidity [GV];#beta");
        hNaF->GetXaxis()->SetTitleOffset(1.2);
        hNaF->GetYaxis()->SetTitleOffset(1.3);
    }
    DrawTheoryLines(c, 0.75, 1.005);
    c->SaveAs(Form("%sBeta_Cutoff_NaF.png", outputDir));

    // AGL部分
    c->Clear();
    tree->Draw("richBeta:cutOffRig>>hAGL(250,0,30,250,0.95,1.005)",
               "(cutStatus & 0x04) == 0x04 && !rich_NaF", "colz");
    TH2D *hAGL = (TH2D*)gDirectory->Get("hAGL");
    if (hAGL) {
        hAGL->SetTitle("Beta vs Cutoff (AGL);Cutoff Rigidity [GV];#beta");
        hAGL->GetXaxis()->SetTitleOffset(1.2);
        hAGL->GetYaxis()->SetTitleOffset(1.3);
    }
    DrawTheoryLines(c, 0.95, 1.005);
    c->SaveAs(Form("%sBeta_Cutoff_AGL.png", outputDir));

    // 清理
    f->Close();
    delete c;
}

void DrawBetaCutoff() {
    Beta_Cutoff();
}