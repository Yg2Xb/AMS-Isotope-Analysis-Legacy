#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <vector>
#include <iostream>

void DrawMCBkgStudy() {
    std::vector<std::string> fileNames = {
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag_4.8.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag_4.8to5.4.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag_SDIATQ.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag_fitstudy.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag.root"
    };
    std::vector<std::string> labels = {
        "L1Q: 4.8-5.5", "L1Q: 4.8-5.4", "SDIAT Q Cut", "L1Q: 4.5-5.5", "old"
    };
    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kBlack};
    std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

    for (const auto& det : detectors) {
        TCanvas* c = new TCanvas(Form("c_%s", det.c_str()), Form("B11 to Be10 Ratio (%s)", det.c_str()), 900, 700);
        TLegend* legend = new TLegend(0.65, 0.65, 0.89, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        std::vector<TFile*> files;
        std::vector<TH1F*> hists;

        // 先全部打开文件、取出hist
        for (size_t i = 0; i < fileNames.size(); i++) {
            TFile* file = TFile::Open(fileNames[i].c_str());
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << fileNames[i] << std::endl;
                files.push_back(nullptr);
                hists.push_back(nullptr);
                continue;
            }
            files.push_back(file);

            TH1F* hist = (TH1F*)file->Get(Form("SourceToIsoRatio_%s_Be10", det.c_str()));
            if (!hist) {
                std::cerr << "Cannot find histogram SourceToIsoRatio_" << det << "_Be10 in " << fileNames[i] << std::endl;
                hists.push_back(nullptr);
                continue;
            }
            // 让hist留在内存（不close file之前不要让file析构）
            hist->SetDirectory(0);
            hists.push_back(hist);
        }

        // 绘图
        bool first = true;
        for (size_t i = 0; i < hists.size(); i++) {
            if (!hists[i]) continue;
            hists[i]->SetMarkerColor(colors[i]);
            hists[i]->SetLineColor(colors[i]);
            if (first) {
                hists[i]->SetTitle(Form("B11 to Be10 Fragmentation Ratio (%s);E_{k}/n [GeV/n];L2 Fragmented Be10 / L1 B11", det.c_str()));
                hists[i]->Draw("PE");
                first = false;
            } else {
                hists[i]->Draw("PE SAME");
            }
            legend->AddEntry(hists[i], labels[i].c_str(), "EP");
        }
        legend->Draw();
        c->SaveAs(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MCBe10_B11_Ratio_%s.png", det.c_str()));

        // 清理
        for (auto file : files) if (file) file->Close();
        delete c;
    }
    std::cout << "比较完成，图像已保存。" << std::endl;
}