#include "../Tool.h"
#include <memory>
#include <vector>

using namespace AMS_Iso;

// 简单的同位素信息结构体
struct IsotopeInfo {
    const char* name;
    int Z{4}, A;
    std::unique_ptr<TH1D> flux;
    Color_t color;
    double fracYmin;
    double fracYmax;

    // 构造函数
    IsotopeInfo(const char* n, int z, int a, Color_t c, double ymin, double ymax)
        : name(n), Z(z), A(a), color(c), fracYmin(ymin), fracYmax(ymax) {}

    // 移动构造函数
    IsotopeInfo(IsotopeInfo&& other) noexcept
        : name(other.name), Z(other.Z), A(other.A),
          flux(std::move(other.flux)),
          color(other.color),
          fracYmin(other.fracYmin),
          fracYmax(other.fracYmax) {}
};

void DrawSun_new() {
    // 通用画图样式设置
    auto SetStyle = [](auto* h, const char* ytitle, double ymin = -99, double ymax = -99) {
        h->SetTitle("");
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1);
        h->GetXaxis()->SetTitle("Ek/n [GeV/n]");
        for (auto axis : {h->GetXaxis(), h->GetYaxis()}) {
            axis->SetTitleSize(0.065);
            axis->SetLabelSize(0.055);
            axis->SetTitleFont(62);
            axis->SetLabelFont(62);
            axis->SetTickLength(0.03);
        }
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitle(ytitle);
        h->GetYaxis()->SetNdivisions(508);
        if (ymin >= -90) h->GetYaxis()->SetRangeUser(ymin, ymax);
    };

    // 初始化
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    const std::string outDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/SunResult/";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    // 初始化同位素信息
    std::vector<IsotopeInfo> isotopes;
    isotopes.emplace_back("Be7", 4, 7, kBlue, 0.4, 0.8);
    isotopes.emplace_back("Be9", 4, 9, kGreen+2, 0.1, 0.6);
    isotopes.emplace_back("Be10", 4, 10, kOrange+1, 0.0, 0.16);

    // 打开文件并获取直方图
    std::unique_ptr<TFile> f(TFile::Open("/eos/ams/user/z/zetong/ISOTOPE/Results/Sun_results.root"));
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // 获取直方图
    for (auto& iso : isotopes) {
        iso.flux.reset((TH1D*)f->Get(Form("szt_fluxbe%d", iso.A)));
        if (!iso.flux) {
            std::cerr << "Error: isotope flux is null" << std::endl;
            return;
        }
    }

    // 创建画布并绘制flux
    std::unique_ptr<TCanvas> c(new TCanvas("c", "c", 800, 600));
    c->SetLogx();
    c->SetLogy(0);
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.165);

    // 画原始flux
    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;

        c->Clear();
        SetStyle(iso.flux.get(), Form("#Phi^{%d}Be [m^{-2}sr^{-1}s^{-1}GeV/n^{-1}]", iso.A));
        iso.flux->SetMarkerColor(iso.color);
        iso.flux->SetLineColor(iso.color);
        iso.flux->SetMarkerSize(1.);
        iso.flux->SetLineWidth(2);
        iso.flux->Draw("PE");
        c->SaveAs(Form("%sflux_%s.png", outDir.c_str(), iso.name));
    }

    // 画flux × (bin center)^3
    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;

        // 构造新TH1D
        auto* h3 = (TH1D*)iso.flux->Clone(Form("flux3_%s", iso.name));
        h3->Reset();
        for (int i = 1; i <= iso.flux->GetNbinsX(); ++i) {
            double val = iso.flux->GetBinContent(i);
            double err = iso.flux->GetBinError(i);
            double cen = iso.flux->GetBinCenter(i);
            h3->SetBinContent(i, val * pow(cen,1));
            h3->SetBinError(i, err * pow(cen,1));
        }

        c->Clear();
        SetStyle(h3, Form("#Phi^{%d}Be #times (E_{k}/n) [m^{-2}sr^{-1}s^{-1}]", iso.A));
        h3->SetMarkerColor(iso.color);
        h3->SetLineColor(iso.color);
        h3->SetMarkerSize(1.);
        h3->SetLineWidth(2);
        h3->GetYaxis()->SetRangeUser(h3->GetMinimum()*0.01, h3->GetMaximum()*1.45);
        h3->Draw("PE");
        c->SaveAs(Form("%sflux3_%s.png", outDir.c_str(), iso.name));

        delete h3;
    }

    // 计算总flux（保持不变）
    auto totalFlux = std::unique_ptr<TH1D>((TH1D*)isotopes[0].flux->Clone("totalFlux"));
    totalFlux->Reset();

    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;
        totalFlux->Add(iso.flux.get());
    }

    // 分别绘制每个同位素的fraction
    c->SetLogy(0);
    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;

        c->Clear();
        auto fraction = std::unique_ptr<TH1D>((TH1D*)iso.flux->Clone(Form("frac_%s", iso.name)));
        fraction->GetListOfFunctions()->Clear();  
        fraction->Divide(totalFlux.get());
        SetStyle(fraction.get(), Form("#Phi^{%d}Be Fraction", iso.A), iso.fracYmin, iso.fracYmax);
        fraction->SetMarkerColor(iso.color);
        fraction->SetLineColor(iso.color);
        fraction->Draw("PE");
        c->SaveAs(Form("%sfraction_%s.png", outDir.c_str(), iso.name));
    }
}