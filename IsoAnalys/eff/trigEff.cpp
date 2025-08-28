#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

void setStyle() {
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(1.);
    gStyle->SetLineWidth(1);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
}

void trigEff() {
    setStyle();
    
    // 只处理Ber文件
    const char* inputFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Li6.root";
    
    TFile* f = TFile::Open(inputFile);
    if(!f) {
        cout << "Cannot open file " << inputFile << endl;
        return;
    }
    
    TH1D* hPhysTrig = (TH1D*)f->Get("Event_PhysTrig_Lithium");
    TH1D* hUnPhys = (TH1D*)f->Get("Event_UnPhys_Lithium");
    //hPhysTrig->Rebin(3);
    //hUnPhys->Rebin(3);

    
    //hUnPhys->Scale(100.);
    
    // 制作总事例数直方图
    TH1D* hTotal = (TH1D*)hUnPhys->Clone("Total");
    hTotal->Add(hPhysTrig);
    
    // 计算效率
    TH1D* hEff = (TH1D*)hPhysTrig->Clone("Eff");
    hEff->Divide(hPhysTrig, hTotal, 1., 1., "B");
    hEff->GetXaxis()->SetRangeUser(0.8,1000);
    
    TCanvas* c = new TCanvas("c", "Trigger Efficiency", 800, 600);
    c->SetLogx();
    
    hEff->SetMarkerStyle(20);
    hEff->SetMarkerColor(kRed);
    hEff->SetLineColor(kRed);
    hEff->GetYaxis()->SetRangeUser(0.98, 1.001);
    hEff->GetYaxis()->SetTitle("Trigger Efficiency");
    hEff->GetXaxis()->SetTitle("Rigidity [GV]");
    hEff->SetTitle("Li6 Trigger Efficiency");
    
    hEff->Draw("P");
    
    c->SaveAs("TriggerEfficiency_Li6.png");
    
    // 保存到root文件
    TFile* outFile = new TFile("trigEff_Li6.root", "RECREATE");
    hEff->Write();
    outFile->Close();
    
    f->Close();
    delete c;
}