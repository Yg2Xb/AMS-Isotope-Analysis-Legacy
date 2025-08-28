#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>

void RichBetaAngle() {
    // 输入文件路径
    TString inputFiles[] = {
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Li6.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Li7.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Lit.root"
    };
    
    // 输出文件路径
    TString outputFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/richBeta/RichBetaAngle.pdf";
    
    // 直方图名称
    const char* histNames[] = {
        "RBeta_Phi_AGL_Lithium", "RBeta_Phi_NaF_Lithium",
        "RBeta_Theta_AGL_Lithium", "RBeta_Theta_NaF_Lithium"
    };
    const int nHists = 4;
    
    // 创建多页 PDF 文件
    TCanvas *c = new TCanvas("c", "RichBetaAngle", 800, 600);
    c->SetLogz();  // 设置 logz 轴
    c->SetRightMargin(0.2);  // 设置左边距为 2.0（20%）
    c->Print((outputFile + "[").Data());  // 打开 PDF 文件
    
    // 遍历文件
    for (int f = 0; f < 3; ++f) {
        TFile *file = TFile::Open(inputFiles[f]);
        if (!file || file->IsZombie()) {
            std::cerr << "无法打开文件: " << inputFiles[f] << std::endl;
            continue;
        }
        
        // 文件名称前缀
        TString prefix;
        if (f == 0) prefix = "Li6";
        else if (f == 1) prefix = "Li7";
        else prefix = "Lit";

        // 遍历直方图
        for (int i = 0; i < nHists; ++i) {
            // 获取直方图
            TH2D *hist = (TH2D*)file->Get(histNames[i]);
            if (!hist) {
                std::cerr << "无法找到直方图: " << histNames[i] << " 在文件: " << inputFiles[f] << std::endl;
                continue;
            }
            if(f==2)
            {
                hist->RebinX(1);
                hist->RebinY(1);
            }
            
            // 设置直方图标题
            if (TString(histNames[i]).Contains("AGL")) {
                hist->SetTitle(prefix + " InnerRig > 100GV");  // 设置 AGL 的标题
                if (TString(histNames[i]).Contains("Phi")) {
                    hist->GetYaxis()->SetRangeUser(0.997, 1.003);     // Y 轴范围 for phi
                }
                if (TString(histNames[i]).Contains("Theta")) {
                    hist->GetXaxis()->SetRangeUser(2.79, 3.15);       // X 轴范围 for theta
                    hist->GetYaxis()->SetRangeUser(0.99, 1.02);      // Y 轴范围 for theta
                }
            } else if (TString(histNames[i]).Contains("NaF")) {
                hist->SetTitle(prefix + " InnerRig > 50GV");  // 设置 NaF 的标题
                if (TString(histNames[i]).Contains("Phi")) {
                    hist->GetYaxis()->SetRangeUser(0.993, 1.007);     // Y 轴范围 for phi
                }
                if (TString(histNames[i]).Contains("Theta")) {
                    //hist->RebinY(2);
                    hist->GetXaxis()->SetRangeUser(2.79, 3.15);       // X 轴范围 for theta
                    hist->GetYaxis()->SetRangeUser(0.99, 1.02);      // Y 轴范围 for theta
                }
            }

            // 绘图
            c->cd();
            hist->Draw("COLZ");  // 绘制直方图
            c->Print(outputFile);  // 输出到 PDF 文件
        }
        
        file->Close();  // 关闭文件
        delete file;
    }
    
    c->Print((outputFile + "]").Data());  // 关闭 PDF 文件
    delete c;
    
    std::cout << "绘图完成并保存到: " << outputFile << std::endl;
}
