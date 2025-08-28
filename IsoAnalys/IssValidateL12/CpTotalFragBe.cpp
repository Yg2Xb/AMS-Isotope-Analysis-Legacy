#include <TFile.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>

// --- 带误差结构体 ---
struct ValueWithError {
    double v=0, e=0;
    ValueWithError() = default;
    ValueWithError(double v_, double e_) : v(v_), e(e_) {}
    ValueWithError operator-(const ValueWithError &o) const { return {v-o.v, std::hypot(e,o.e)}; }
    ValueWithError operator/(const ValueWithError &o) const {
        if(o.v==0) return {0,0};
        double q=v/o.v;
        return {q, std::hypot(e/o.v, (v*o.e)/(o.v*o.v))};
    }
    ValueWithError operator*(const ValueWithError& o) const {
        double p_val = v * o.v;
        double p_err = std::hypot(e * o.v, v * o.e);
        return {p_val, p_err};
    }
    ValueWithError operator*(double scalar) const { return {v*scalar, std::abs(e*scalar)}; }
    friend ValueWithError operator*(double scalar, const ValueWithError& ve) { return ve*scalar; }
};

// --- 文件操作 ---
std::unique_ptr<TFile> openFile(const char* fn, const char* opt="READ") {
    auto f=std::unique_ptr<TFile>(TFile::Open(fn,opt));
    if(!f||f->IsZombie()) throw std::runtime_error("OpenFileFail: "+std::string(fn));
    return f;
}
TH1* getHist(TFile* f, const TString& n) {
    TObject* obj = f->Get(n);
    if (!obj) throw std::runtime_error("HistNotFound: "+std::string(n.Data()));
    auto* h = dynamic_cast<TH1*>(obj);
    if (!h) throw std::runtime_error("Object is not TH1: "+std::string(n.Data()));
    return h;
}
// --- 类型安全hist获取 ---
template<typename H>
H* getHistChecked(TFile* f, const TString& n) {
    auto* h = dynamic_cast<H*>(getHist(f, n));
    if(!h) throw std::runtime_error("Hist "+std::string(n.Data())+" not matching type "+typeid(H).name());
    return h;
}

// --- 合并三探测器 ---
TH1F* mergeDetectors(TH1F* h_tof, TH1F* h_naf, TH1F* h_agl, const char* name) {
    int n = h_tof->GetNbinsX(); float* bins = new float[n+1];
    for(int i=0;i<=n;++i) bins[i]=h_tof->GetBinLowEdge(i+1);
    auto h = new TH1F(name, name, n, bins); delete[] bins;
    for(int i=1;i<=n;++i) {
        double e=h->GetBinCenter(i);
        cout<<e<<endl;
        if(e<1.17)      {h->SetBinContent(i, h_tof->GetBinContent(i)); h->SetBinError(i,h_tof->GetBinError(i));}
        else if(e<3.23) {h->SetBinContent(i, h_naf->GetBinContent(i)); h->SetBinError(i,h_naf->GetBinError(i));}
        else            {h->SetBinContent(i, h_agl->GetBinContent(i)); h->SetBinError(i,h_agl->GetBinError(i));}
    } return h;
}

// --- 合并权重 ---
struct MCWeight { double w11=1, w10=1; };
// 修改权重函数，接受能量bin作为参数
MCWeight getMCWeight(double B11Frac, const char* f11, const char* f10) {
    auto b11 = openFile(f11), b10 = openFile(f10);
    auto h11 = getHistChecked<TH1D>(b11.get(), "hMCnum");
    auto h10 = getHistChecked<TH1D>(b10.get(), "hMCnum");
    double n11=h11->GetBinContent(1), n10=h10->GetBinContent(1), s=n11+n10;
    if(n11<=0||n10<=0) throw std::runtime_error("hMCnum bin content <=0");
    return {s/n11*B11Frac, s/n10*(1.0-B11Frac)};
}

// 修改加权合并直方图函数，接受B11分数直方图作为参数
template <typename H>
H* mixHistWithFraction(const char* name, H* h11, H* h10, TH1D* fracHist) {
    if(!h11 || !h10) throw std::runtime_error("mixHistWithFraction: input hist is nullptr!");
    if(!fracHist) throw std::runtime_error("mixHistWithFraction: fracHist is nullptr!");
    
    gROOT->cd();
    H* h = (H*)h11->Clone(name);
    h->SetDirectory(0);
    h->Reset();
    
    for(int bin=1; bin<=h11->GetNbinsX(); bin++) {
        double energy = h11->GetBinCenter(bin);
        int fracBin = fracHist->FindBin(energy);
        double B11Frac = fracHist->GetBinContent(fracBin);
        if(B11Frac <= 0) B11Frac = 0.7; // 默认值
        
        // 获取B11和B10的值和误差
        double val11 = h11->GetBinContent(bin);
        double err11 = h11->GetBinError(bin);
        double val10 = h10->GetBinContent(bin);
        double err10 = h10->GetBinError(bin);
        
        // 计算MC生成数量以获取权重
        auto b11 = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_bkg.root");
        auto b10 = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_bkg.root");
        auto mcNum11 = getHistChecked<TH1D>(b11.get(), "hMCnum");
        auto mcNum10 = getHistChecked<TH1D>(b10.get(), "hMCnum");
        double n11 = mcNum11->GetBinContent(1);
        double n10 = mcNum10->GetBinContent(1);
        double s = n11 + n10;
        
        // 计算加权值和误差
        double weight11 = s/n11 * B11Frac;
        double weight10 = s/n10 * (1.0-B11Frac);
        double mixedVal = weight11 * val11 + weight10 * val10;
        double mixedErr = sqrt(pow(weight11 * err11, 2) + pow(weight10 * err10, 2));
        
        // 设置bin内容
        h->SetBinContent(bin, mixedVal);
        h->SetBinError(bin, mixedErr);
    }
    
    return h;
}

// 修改获取MC比率函数，使用B11分数直方图
void getMCratios(const char* f11, const char* f10, TH1D* fracHist, const char* det, TH1F*& h_cut, TH1F*& h_truth) {
    auto b11 = openFile(f11), b10 = openFile(f10);
    auto h_cut11   = getHistChecked<TH1F>(b11.get(), Form("FragNucCounts_%s", det));
    auto h_cut10   = getHistChecked<TH1F>(b10.get(), Form("FragNucCounts_%s", det));
    auto h_b11     = getHistChecked<TH1F>(b11.get(), Form("TruthL1SourceCounts_%s", det));
    auto h_b10     = getHistChecked<TH1F>(b10.get(), Form("TruthL1SourceCounts_%s", det));
    auto h_be711   = getHistChecked<TH1F>(b11.get(), Form("SourceToIsoCounts_%s_Be7", det));
    auto h_be911   = getHistChecked<TH1F>(b11.get(), Form("SourceToIsoCounts_%s_Be9", det));
    auto h_be1011  = getHistChecked<TH1F>(b11.get(), Form("SourceToIsoCounts_%s_Be10", det));
    auto h_be710   = getHistChecked<TH1F>(b10.get(), Form("SourceToIsoCounts_%s_Be7", det));
    auto h_be910   = getHistChecked<TH1F>(b10.get(), Form("SourceToIsoCounts_%s_Be9", det));
    auto h_be1010  = getHistChecked<TH1F>(b10.get(), Form("SourceToIsoCounts_%s_Be10", det));
    
    // 使用分数直方图合并碎裂Be
    h_cut = mixHistWithFraction(Form("MC_CutBe_%s",det), h_cut11, h_cut10, fracHist);
    TH1F* h_L1B = mixHistWithFraction(Form("MC_L1B_%s",det), h_b11, h_b10, fracHist);
    
    // 合并真Be
    TH1F* h_truth11 = (TH1F*)h_be711->Clone(); h_truth11->Add(h_be911); h_truth11->Add(h_be1011);
    TH1F* h_truth10 = (TH1F*)h_be710->Clone(); h_truth10->Add(h_be910); h_truth10->Add(h_be1010);
    h_truth = mixHistWithFraction(Form("MC_TruthBe_%s",det), h_truth11, h_truth10, fracHist);
    
    // 输出一些调试信息
    if(true) {
        cout << "Check MC mixing for " << det << " with B11 fraction histogram" << endl;
        for(int bin=1; bin<=h_truth->GetNbinsX(); bin++) {
            double energy = h_truth->GetBinCenter(bin);
            int fracBin = fracHist->FindBin(energy);
            double B11Frac = fracHist->GetBinContent(fracBin);
            if(B11Frac <= 0) B11Frac = 0.7;
            cout << "Bin " << bin << ", Energy " << energy << ", B11Frac " << B11Frac;
            cout << ", Cut " << h_cut->GetBinContent(bin) << ", Truth " << h_truth->GetBinContent(bin) << endl;
        }
    }
    
    // 除以L1上B
    h_cut->Divide(h_L1B); 
    h_truth->Divide(h_L1B);
}

void CpTotalFragBe() {
    try {
        // ========== 加载B11分数直方图 ==========
        auto B11FracFile = openFile("/eos/ams/user/z/zetong/Boron_isotope/pre/B11f1.root");
        auto h_B11_Frac = getHistChecked<TH1D>(B11FracFile.get(), "B11");
        
        // ========== ISS ==========
        auto iss_counts = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_temp_wide_bkg_L1Inner510cutoff.root");
        auto iss_ratiof = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_combined.root");
        std::array<const char*,3> dets{"TOF","NaF","AGL"};
        std::vector<TH1F*> iss_num(3), iss_den(3), iss_ratio(3);
        
        for(int k=0;k<3;++k) {
            auto h_frag = getHistChecked<TH1F>(iss_counts.get(), Form("FragNucCounts_%s",dets[k]));
            auto h_boron = getHistChecked<TH1F>(iss_counts.get(), Form("CutL1SourceCounts_%s",dets[k]));
            auto h_be_ratio = getHistChecked<TH1D>(iss_ratiof.get(), "h_be_fraction");
            auto h_boron_ratio = getHistChecked<TH1D>(iss_ratiof.get(), "h_b_fraction");
            
            for(int i=1;i<=h_frag->GetNbinsX();++i) {
                double rawB = h_boron->GetBinContent(i), rawB_err = std::sqrt(std::max(0.0,rawB));
                double rawBe = h_frag->GetBinContent(i), rawBe_err = std::sqrt(std::max(0.0,rawBe));
                double beRat = h_be_ratio->GetBinContent(i), beRat_err = h_be_ratio->GetBinError(i);
                double bRat = h_boron_ratio->GetBinContent(i), bRat_err = h_boron_ratio->GetBinError(i);
                ValueWithError contamBe = ValueWithError(rawB, rawB_err) * ValueWithError(beRat, beRat_err);
                ValueWithError boron_final = ValueWithError(rawB, rawB_err) * ValueWithError(bRat, bRat_err);
                ValueWithError be_num = ValueWithError(rawBe,rawBe_err) - contamBe;
                
                if(i==1) {
                    cout<<"check iss"<<endl;
                    cout<<rawB<<endl;
                    cout<<rawBe<<endl;
                    cout<<bRat<<endl;
                    cout<<beRat<<endl;
                    cout<<be_num.v<<endl;
                }
                
                if(!iss_num[k]) iss_num[k] = (TH1F*)h_frag->Clone(Form("ISS_NumBe_%s",dets[k]));
                if(!iss_den[k]) iss_den[k] = (TH1F*)h_boron->Clone(Form("ISS_DenB_%s",dets[k]));
                iss_num[k]->SetBinContent(i, be_num.v); iss_num[k]->SetBinError(i,be_num.e);
                iss_den[k]->SetBinContent(i, boron_final.v); iss_den[k]->SetBinError(i, boron_final.e);
            }
            
            iss_ratio[k] = (TH1F*)iss_num[k]->Clone(Form("ISS_Be_B_Ratio_%s",dets[k]));
            iss_ratio[k]->Divide(iss_den[k]);
        }
        
        auto iss_merge = mergeDetectors(iss_ratio[0], iss_ratio[1], iss_ratio[2], "ISS_Be_B_Ratio");
        iss_merge->SetLineColor(kBlack); iss_merge->SetMarkerStyle(20);

        // ========== MC ==========
        std::vector<TH1F*> mc_cut(3, nullptr), mc_truth(3, nullptr);
        
        for(int k=0;k<3;++k){
            getMCratios(
                "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_temp_wide_bkg_frag_L1Inner.root",
                "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_temp_wide_bkg_frag_L1Inner.root",
                h_B11_Frac, dets[k], mc_cut[k], mc_truth[k]
            );
            
            if(!mc_cut[k]) { std::cerr<<"mc_cut["<<k<<"] is nullptr!"<<std::endl; }
            else { std::cerr<<"mc_cut["<<k<<"] not nullptr"<<std::endl; }
        }
        
        auto mc_cut_merge   = mergeDetectors(mc_cut[0], mc_cut[1], mc_cut[2], "MC_CutBe_B_Ratio");
        auto mc_truth_merge = mergeDetectors(mc_truth[0], mc_truth[1], mc_truth[2], "MC_TruthBe_B_Ratio");
        mc_cut_merge->SetLineColor(kBlue); mc_cut_merge->SetMarkerColor(kBlue);
        mc_truth_merge->SetLineColor(kRed); mc_truth_merge->SetMarkerColor(kRed);

        // ========== 输出画图 ==========
        gStyle->SetOptStat(0); gStyle->SetErrorX(0);
        auto c = new TCanvas("c","Be/B Comparison",900,700);
        mc_truth_merge->Draw("pz");
        mc_truth_merge->SetTitle("L2 Fragmented Be / L1 Selected Boron;E_{k}/n[GeV/n]; Ratio");
        mc_truth_merge->GetXaxis()->SetRangeUser(0,16.3);
        mc_truth_merge->GetYaxis()->SetRangeUser(0,0.02);
        //mc_cut_merge->Draw("E1X0 SAME");
        iss_merge->Draw("pz SAME");
        
        auto leg = new TLegend(0.55,0.65,0.87,0.87);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(0, 0.6);
        leg->SetTextSize(0.033);
        leg->AddEntry(iss_merge,"ISS Data","ep");
        //leg->AddEntry(mc_cut_merge,"MC mix(cut Be/B)","ep");
        leg->AddEntry(mc_truth_merge,"B10 B11 MC mix","ep");
        leg->Draw();
        
        c->SetGrid();
        c->SetTopMargin(0.1);
        c->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/TotalBeB_comp_yzxchargefit.png");
        std::cout << "Saved: .../TotalBeB_comp.png/root" << std::endl;
    } catch(const std::exception& e) {
        std::cerr << "\nException in CpTotalFragBe: " << e.what() << std::endl;
        throw;
    }
}