#include "../basic_var.h"
#include "../Tool.h"


using namespace AMS_Iso;

struct IsotopeConfig {
    std::string name;
    int charge;
    std::vector<int> masses;
    std::vector<std::string> legends;
    std::vector<int> colors;

    IsotopeConfig(const std::string& elementName) {
        // 从IsotopeData查找对应元素
        for (const auto& iso : IsotopeData) {
            if (ConvertElementName(iso.name_, false) == elementName) {
                name = elementName;
                charge = iso.charge_;
                // 添加非零质量数
                for (int i = 0; i < iso.isotope_count_; ++i) {
                    if (iso.mass_[i] != 0) {
                        masses.push_back(iso.mass_[i]);
                        legends.push_back(elementName + std::to_string(iso.mass_[i]));
                        colors.push_back(iso.color_[i]);
                    }
                }
                break;
            }
        }
    }
};

void plotIsoET(const std::string& element, const std::vector<std::string>& inputFiles) {
    // 获取同位素配置
    IsotopeConfig config(element);
    if (config.masses.empty()) {
        std::cerr << "Error: Unknown element " << element << std::endl;
        return;
    }

    const std::vector<std::string> detectors = {"TOF", "NaF", "Aero"};
    const int nMasses = config.masses.size();
    const int nDetectors = detectors.size();

    // 创建存储直方图的二维数组
    std::vector<std::vector<TH1D*>> histograms(nMasses);
    for (int i = 0; i < nMasses; ++i) {
        // 打开对应质量数的文件
        TFile* file = TFile::Open(inputFiles[i].c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Unable to open file " << inputFiles[i] << std::endl;
            return;
        }

        histograms[i].resize(nDetectors);
        for (int j = 0; j < nDetectors; ++j) {
            std::string histName = Form("HExpBeta_M%d_%s", config.masses[i], detectors[j].c_str());
            histograms[i][j] = (TH1D*)file->Get(histName.c_str());
            histograms[i][j]->SetTitle(""); 
            if (!histograms[i][j]) {
                std::cerr << "Error: Cannot find histogram " << histName << std::endl;
                return;
            }
            // 克隆直方图以避免文件关闭后的问题
            histograms[i][j] = (TH1D*)histograms[i][j]->Clone();
        }
    }

    // 创建画布
    TCanvas *c1 = new TCanvas("c1", "Exposure Time Comparison", 1200, 800);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.1);
    c1->SetTopMargin(0.1);
    c1->SetBottomMargin(0.15);
    c1->Divide(3, 1);

    // 从BetaTypes获取能量范围
    std::vector<std::pair<double, double>> xRanges;
    for (const auto& beta : Detector::BetaTypes) {
        xRanges.push_back({std::max(beta.Ekn_range_[0],0.5), beta.Ekn_range_[1]});
    }

    // 绘制部分保持不变...
    for (int det = 0; det < nDetectors; ++det) {
        c1->cd(det + 1);
        
        gPad->SetLeftMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetLogx();
        
        TLegend *legend = new TLegend(0.21, 0.65, 0.52, 0.87);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.055);

        for (int mass = 0; mass < nMasses; ++mass) {
            TH1D* hist = histograms[mass][det];
            hist->SetLineColor(config.colors[mass]);
            hist->SetLineWidth(3);
            hist->SetStats(0);
            
            hist->GetXaxis()->SetRangeUser(det == 0? 0.45 : xRanges[det].first, xRanges[det].second);
            hist->GetYaxis()->SetRangeUser(0, 320e6);
            
            hist->GetXaxis()->SetTitleOffset(0.5);
            hist->GetYaxis()->SetTitleOffset(1.25);
            hist->GetXaxis()->SetLabelOffset(-0.05);
            hist->GetXaxis()->SetLabelSize(0.08);
            hist->GetYaxis()->SetLabelSize(0.08);
            hist->GetXaxis()->SetTitleSize(0.08);
            hist->GetYaxis()->SetTitleSize(0.08);
            
            hist->GetXaxis()->SetTitle("Ek/n[GeV/n]");
            hist->GetYaxis()->SetTitle("");
            hist->SetTitle("");

            if (mass == 0) {
                hist->Draw("hist");
            } else {
                hist->Draw("hist same");
            }

            legend->AddEntry(hist, Form("%s %s", config.legends[mass].c_str(), 
                                      detectors[det].c_str()), "l");
        }
        legend->Draw();
    }

    // 保存图像
    std::string outputFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/ExpoTime/" + 
                            element + "_Ek.pdf";
    c1->Print(outputFile.c_str());

    // 清理内存
    for (auto& histVec : histograms) {
        for (auto* hist : histVec) {
            delete hist;
        }
    }
}

void plotET()
{
    // 使用示例：处理Be的三个同位素
    std::vector<std::string> beFiles = {
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root"
    };
    std::vector<std::string> bFiles = {
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root",
        "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root"
    };
    plotIsoET("Be", beFiles);
    //plotIsoET("B", bFiles);
}