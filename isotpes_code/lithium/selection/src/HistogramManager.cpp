/***********************************************************
 *  File: HistogramManager.cpp
 *
 *  Implementation file for AMS Isotopes histogram management.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include "HistogramManager.h"
#include <stdexcept>
#include <iostream>
#include "TOFCut.h"  
#include "RICHCut.h"

namespace AMS_Iso {
// Forward declarations
class IsotopeVar;  // 只需要这一个前向声明

void HistogramManager::createHistograms(const IsotopeVar* isotope, int useMass) {
    if (!isotope) {
        throw std::runtime_error("Null isotope pointer in createHistograms");
    }
    isotope_ = isotope;
    useMass_ = useMass; 

    createMCHistograms();
    createExposureHistograms();
    createEventAndKineticHistograms();
    createResolutionHistograms();
    createEfficiencyHistograms();
    //createChargeTempFitHistograms();
}

void HistogramManager::createMCHistograms() {
    auto& collection = histograms[HistType::MC];
    collection.hist1D["MCNum"] = std::make_unique<TH1D>(
        "hMCnum", "MC_TotNumber;;", 1, 1, 2);
}

void HistogramManager::createExposureHistograms() {
    auto& collection = histograms[HistType::Exposure];
    
    // 创建Rigidity曝光时间直方图
    collection.hist1D["ExpRig"] = std::make_unique<TH1D>(
        "HExposureRig",
        "Histogram for Exposure Time in Rigidity Bins;Rig[GV];ExpoTime[s]",
        Constants::RIGIDITY_BINS,
        Binning::RigidityBins.data()  // 使用bin边界数组
    );

    // 创建Beta曝光时间直方图
    auto& betaCollection = histograms[HistType::BetaExposure];
    betaCollection.arrays.resize(isotope_->getIsotopeCount());
    for (int is = 0; is < isotope_->getIsotopeCount(); ++is) {
        //2025.3.4, only expotime in usemass bin
        //int charge = isotope_->getCharge();
        //int IsotopeIndex = findIsotopeIndex(useMass_, charge);
        //for (int is = IsotopeIndex; is <= IsotopeIndex; ++is) {
        betaCollection.arrays[is].resize(Constants::BETA_TYPES);
        for (int ib = 0; ib < Constants::BETA_TYPES; ++ib) {
            const auto& kinEnergyBins = Binning::KineticEnergyBins[isotope_->getCharge()-1][is];
            betaCollection.arrays[is][ib] = std::make_unique<TH1D>(
                Form("HExpBeta_M%d_%s",
                     isotope_->getMass(is),
                     Detector::BetaTypes[ib].getName().c_str()),
                Form("Histogram for Exposure Time in %s Beta for Mass%d;Ek[GeV/n];ExpoTime[s]",
                     Detector::BetaTypes[ib].getName().c_str(),
                     isotope_->getMass(is)),
                Constants::RIGIDITY_BINS,
                kinEnergyBins.data()  // 使用bin边界数组
            );
        }
    }
}

//acceptance
void HistogramManager::createEventAndKineticHistograms() {
    auto& collection = histograms[HistType::Event];
    const std::string& istp_name = isotope_->getName();
    int charge = isotope_->getCharge();
    std::cout<<"Using Mass = "<<useMass_<<" FOR BINNING"<<std::endl;
    int IsotopeIndex = findIsotopeIndex(useMass_, charge);

    // 创建基本事例计数直方图 (Generated Ek)
    collection.hist1D["EventEk"] = std::make_unique<TH1D>(
        Form("Event_Ek_%s", istp_name.c_str()),
        Form("Hist for Events Count in Ek %s;Generated Ek[GeV/n];Events", istp_name.c_str()),
        Constants::RIGIDITY_BINS,
        Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
    );
    
    // TOF事例计数直方图 (Generated Ek)
    TString tfname;
    for (int cn1 = 0; cn1 <= 10; ++cn1) {
        tfname += Detector::TOFCuts[cn1];
        collection.hist1D[Form("EventEkTOF_%d", cn1)] = std::make_unique<TH1D>(
            Form("Event_Ek_TOF_%s_%d", istp_name.c_str(), cn1),
            Form("cut:%s;Generated Ek[GeV/n];Events", tfname.Data()),
            Constants::RIGIDITY_BINS,
            Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
        );
    }

    // RICH事例计数直方图 (Generated Ek)
    TString rhname;
    for (int cn1 = 0; cn1 <= 10; ++cn1) {
        if (cn1 < 10) rhname += Detector::RICHCuts[cn1];
        // NaF
        collection.hist1D[Form("EventEkNaF_%d", cn1)] = std::make_unique<TH1D>(
            Form("Event_Ek_NaF_%s_%d", istp_name.c_str(), cn1),
            Form("cut:%s;Generated Ek[GeV/n];Events", rhname.Data()),
            Constants::RIGIDITY_BINS,
            Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
        );
        
        // AGL
        collection.hist1D[Form("EventEkAGL_%d", cn1)] = std::make_unique<TH1D>(
            Form("Event_Ek_AGL_%s_%d", istp_name.c_str(), cn1),
            Form("cut:%s;Generated Ek[GeV/n];Events", rhname.Data()),
            Constants::RIGIDITY_BINS,
            Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
        );
    }

    // TOF特殊区域动能直方图 (Generated Ek)
    const std::vector<std::string> special_regions = {
        "Exclude4", "NaF", "isNaF", "NaFGeo", "NaFGeo2"
    };

    for (const auto& region : special_regions) {
        collection.hist1D[Form("EventEkTOF_%s", region.c_str())] = std::make_unique<TH1D>(
            Form("Event_Ek_TOF_%s_%s", region.c_str(), istp_name.c_str()),
            Form("TOF Kinetic Energy %s;Generated Ek[GeV/n];Events", region.c_str()),
            Constants::RIGIDITY_BINS,
            Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
        );
    }

    // 实际测量的Ek直方图
    collection.hist1D["MeasuredEkTOF"] = std::make_unique<TH1D>(
        Form("Event_Measured_Ek_TOF_%s", istp_name.c_str()),
        Form("TOF Measured Kinetic Energy;Ek[GeV/n];Events"),
        Constants::RIGIDITY_BINS,
        Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
    );

    collection.hist1D["MeasuredEkNaF"] = std::make_unique<TH1D>(
        Form("Event_Measured_Ek_NaF_%s", istp_name.c_str()),
        Form("NaF Measured Kinetic Energy;Ek[GeV/n];Events"),
        Constants::RIGIDITY_BINS,
        Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
    );

    collection.hist1D["MeasuredEkAGL"] = std::make_unique<TH1D>(
        Form("Event_Measured_Ek_AGL_%s", istp_name.c_str()),
        Form("AGL Measured Kinetic Energy;Ek[GeV/n];Events"),
        Constants::RIGIDITY_BINS,
        Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
    );
}

void HistogramManager::createResolutionHistograms() {
    auto& collection = histograms[HistType::Resolution];
    const std::string& istp_name = isotope_->getName();
    int charge = isotope_->getCharge();
    int IsotopeIndex = findIsotopeIndex(useMass_, charge);

    // Beta分辨直方图
    for (int i = 0; i < 2; ++i) {
        // RICH Beta分辨率
        collection.hist1D[Form("RSL_RichBeta_%s", Detector::RichDetectorNames[i].c_str())] =
            std::make_unique<TH1D>(
                Form("R_%sBeta_%s", Detector::RichDetectorNames[i].c_str(), istp_name.c_str()),
                Form("1/%sBeta_%s;1/Beta_{%s};Events",
                     Detector::RichDetectorNames[i].c_str(),
                     istp_name.c_str(),
                     Detector::RichDetectorNames[i].c_str()),
                Detector::RichBins[i][charge - 1]*100,
                1 - Detector::RichAxis[i],
                1 + Detector::RichAxis[i]
            );
        // RICH Beta to theta and phi
        collection.hist2D[Form("RBeta_Theta_%s", Detector::RichDetectorNames[i].c_str())] =
            std::make_unique<TH2D>(
                Form("RBeta_Theta_%s_%s", Detector::RichDetectorNames[i].c_str(), istp_name.c_str()),
                Form("RBeta_Theta_%s_%s;rich_theta;1/beta_{%s};Events",
                     Detector::RichDetectorNames[i].c_str(),
                     istp_name.c_str(),
                     Detector::RichDetectorNames[i].c_str()),
                1000, 2.7, 3.2,
                1000, 0.95, 1.05
            );
        collection.hist2D[Form("RBeta_Phi_%s", Detector::RichDetectorNames[i].c_str())] =
            std::make_unique<TH2D>(
                Form("RBeta_Phi_%s_%s", Detector::RichDetectorNames[i].c_str(), istp_name.c_str()),
                Form("RBeta_Phi_%s_%s;rich_phi;1/beta_{%s};Events",
                     Detector::RichDetectorNames[i].c_str(),
                     istp_name.c_str(),
                     Detector::RichDetectorNames[i].c_str()),
                1000, -3.2, 3.2,
                1000, 0.95, 1.05
            );

        // TOF-RICH Beta for TOF rsl
        collection.hist2D[Form("RSL_TOFBeta_Ek_%s", Detector::RichDetectorNames[i].c_str())] =
            std::make_unique<TH2D>(
                Form("RBeta_Ek_TOF_%s_%s", Detector::RichDetectorNames[i].c_str(), istp_name.c_str()),
                //Form("RBeta_Ek_TOF_%s_%s;1/beta_{TOF}-1/beta_{%s};Ek [GeV/n];Events",
                Form("RBeta_Ek_TOF_%s_%s;1/beta_{TOF}-1/beta_{%s};Beta*Rig [GV];Events",
                     Detector::RichDetectorNames[i].c_str(),
                     istp_name.c_str(),
                     Detector::RichDetectorNames[i].c_str()),
                Detector::TOFRichBins[i][charge - 1]*100,
                -Detector::TOFAxis[i],
                Detector::TOFAxis[i],
                500, 0.1, 50.1
                //Constants::RIGIDITY_BINS,
                //Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
            );

        // TOF other GeoCut 分辨直方图
        const std::vector<std::string> resolution_types = {
            "AllEdgeTOFBeta_Ek", "NaFGeoTOFBeta_Ek"
        };

        for (const auto& type : resolution_types) {
            collection.hist2D[Form("RSL_%s_%s", type.c_str(), Detector::RichDetectorNames[i].c_str())] =
                std::make_unique<TH2D>(
                    Form("%s_%s_%s", type.c_str(), Detector::RichDetectorNames[i].c_str(), istp_name.c_str()),
                    //Form("%s;1/beta_{TOF}-1/beta_{%s};Ek[GeV/n];Events",type.c_str(),
                    Form("%s;1/beta_{TOF}-1/beta_{%s};Beta*Rig [GV];Events",type.c_str(),
                         Detector::RichDetectorNames[i].c_str()),
                    Detector::TOFRichBins[i][charge - 1]*100,
                    -Detector::TOFAxis[i],
                    Detector::TOFAxis[i],
                    500, 0.1, 50.1 
                    //Constants::RIGIDITY_BINS,
                    //Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
                );
        }
    }
            // AGL Beta分辨率
            collection.hist1D["BetaCheck"] =
            std::make_unique<TH1D>(
                Form("BetaCheck_%s", istp_name.c_str()),
                Form("1/#beta_{truth}-1/#beta_{AGL};1/#beta_{truth}-1/#beta_{AGL};Events"),
                1000, -0.002, 0.002
            );
    
            // Rigidity分辨率  
            collection.hist1D["RigCheck"] =
            std::make_unique<TH1D>(
                Form("RigCheck_%s", istp_name.c_str()),
                Form("1/Rig_{truth}-1/GBLInnerRig;1/Rig_{truth}-1/GBLInnerRig;Events"),
                1000, -0.008, 0.016
            );
    
            // Tracker-AGL Beta分辨率
            collection.hist1D["RigBetaCheck"] =
            std::make_unique<TH1D>(
                Form("RigBetaCheck_%s", istp_name.c_str()),
                Form("1/#beta_{tracker}-1/#beta_{AGL};1/#beta_{tracker}-1/#beta_{AGL};Events"),
                1000, -0.005, 0.005
            );
    
            // Tracker-AGL Beta分辨率 Ek
            collection.hist2D["RigBetaCheck_Ek"] =
            std::make_unique<TH2D>(
                Form("RigBetaCheck_Ek_%s", istp_name.c_str()),
                Form("1/#beta_{tracker}-1/#beta_{AGL} in Ek;1/#beta_{tracker}-1/#beta_{AGL};Ek/n[GeV/n];Events"),
                1000, -0.005, 0.005, 15 - 1, Binning::BkgEkWideBin.data()  
            );
}

void HistogramManager::createEfficiencyHistograms() {
    auto& collection = histograms[HistType::Efficiency];
    const std::string& istp_name = isotope_->getName();
    int charge = isotope_->getCharge();
    int IsotopeIndex = findIsotopeIndex(useMass_, charge);

    // 物理触发直方图
    const std::vector<std::pair<std::string, std::string>> trigger_types = {
        {"PhysTrig", "NucleiCut_PhysTrigger"},
        {"UnPhys", "NucleiCut_UnPhysTrigger"},
        {"TOFPhysTrig", "Nuclei+TOFCut_PhysTrigger"},
        {"TOFUnPhys", "Nuclei+TOFCut_UnPhysTrigger"},
        {"RICHPhysTrig", "Nuclei+RICHCut_PhysTrigger"},
        {"RICHUnPhys", "Nuclei+RICHCut_UnPhysTrigger"}
    };

    for (const auto& [type, title] : trigger_types) {
        collection.hist1D[type] = std::make_unique<TH1D>(
            Form("Event_%s_%s", type.c_str(), "Carbon"),
            Form("%s; InnerRig [GV]; Events", title.c_str()),
            Constants::RIGIDITY_BINS, Binning::RigidityBins.data()
        );
    }

    // 效率直方图
    const std::vector<std::string> eff_types = {
        "L1PickUp", "RICHRec"
    };

    for (const auto& type : eff_types) {
        for (int p = 0; p < 4; ++p) {
            const char* effname = (p % 2 == 0 ? "den" : "num");
            const char* data_type = (p < 2 ? "ISS" : "MC");

            collection.hist1D[Form("%sEff_%s_Ek_%s_%d", 
                                 type.c_str(), effname, data_type, p)] =
                std::make_unique<TH1D>(
                    Form("%sEff_%s_Ek_%s_%s", 
                         type.c_str(), effname, data_type, istp_name.c_str()),
                    Form("%sEff_%s_Ek_%s_%s;Ek[GeV/n];Events",
                         type.c_str(), effname, data_type, istp_name.c_str()),
                    Constants::RIGIDITY_BINS,
                    Binning::KineticEnergyBins[charge-1][IsotopeIndex].data()
                );
        }
    }
}

void HistogramManager::createChargeTempFitHistograms() {
    auto& collection = histograms[HistType::ChargeTempFit];
    const std::string& istp_name = isotope_->getName();
    int charge = isotope_->getCharge();
    int IsotopeIndex = findIsotopeIndex(useMass_, charge);

    const std::vector<std::string> select_types = {"L1Temp", "L2Temp"};
    const std::vector<std::string> l1_types = {"L1Normal", "L1Unbiased"};
    const std::vector<double> coe_values = {0.2, 0.4, 1.0};
    const std::vector<std::string> elements = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    const std::vector<std::string> bkg_types = {"Bkg"};

    auto getXAxisTitle = [](const std::string& select, const std::string& l1_type) {
        if (select == "L2Temp") return "TrackerL2 Q";
        return l1_type == "L1Normal" ? "TrackerL1 Q" : "TrackerL1 UnbiasedQ";
    };

    for (const auto& select : select_types) {
        for (const auto& l1_type : l1_types) {
            for (const auto& coe : coe_values) {
                if(select == "L1Temp" && coe < 1.0) continue;
                for (const auto& elem : elements) {
                    for (const auto& bkg : bkg_types) {
                        std::string coe_str = Form("%.1f", coe);
                        std::replace(coe_str.begin(), coe_str.end(), '.', 'p');
                        
                        std::string name = Form("%s_%s_Coe%s_%s_%s", 
                            select.c_str(), l1_type.c_str(), coe_str.c_str(), elem.c_str(), bkg.c_str());

                        collection.hist2D[name] = std::make_unique<TH2D>(
                            name.c_str(),
                            Form("%s_%s_Coe%s_%s_%s; %s; Ek/n[GeV/n]; Events", 
                                select.c_str(), l1_type.c_str(), coe_str.c_str(), elem.c_str(), bkg.c_str(),
                                getXAxisTitle(select, l1_type)),
                            600, 3, 9, 15 - 1, Binning::BkgEkWideBin.data()
                        );
                    }
                }
            }
        }
    }
    
    //bkg check
    const std::vector<std::string> det_types = {"TOF", "NaF", "AGL"}; 
    for (const auto& det_type : det_types) {
        std::string fragname = Form("FragNucCounts_%s", det_type.c_str());
        collection.hist1D[fragname] = std::make_unique<TH1D>(
            fragname.c_str(), Form("L2 Frag; Ek[GeV/n]; Events"),
            15 - 1, Binning::BkgEkWideBin.data()
        );
        std::string sourcename = Form("CutL1SourceCounts_%s", det_type.c_str());
        collection.hist1D[sourcename] = std::make_unique<TH1D>(
            sourcename.c_str(), Form("L1 Source; Ek[GeV/n]; Events"),
            15 - 1, Binning::BkgEkWideBin.data()
        );
    }

}

TH1* HistogramManager::getHistogram(HistType type, const std::string& name) const {
    auto it = histograms.find(type);
    if (it != histograms.end()) {
        // 先查找1D直方图
        auto it1D = it->second.hist1D.find(name);
        if (it1D != it->second.hist1D.end()) {
            return it1D->second.get();
        }
        // 再查找2D直方图
        auto it2D = it->second.hist2D.find(name);
        if (it2D != it->second.hist2D.end()) {
            return it2D->second.get();
        }
    }
    return nullptr;
}

TH1D* HistogramManager::getHist1D(HistType type, const std::string& name) const {
    auto it = histograms.find(type);
    if (it != histograms.end()) {
        auto hist_it = it->second.hist1D.find(name);
        if (hist_it != it->second.hist1D.end()) {
            return hist_it->second.get();
        }
    }
    return nullptr;
}

TH2D* HistogramManager::getHist2D(HistType type, const std::string& name) const {
    auto it = histograms.find(type);
    if (it != histograms.end()) {
        auto hist_it = it->second.hist2D.find(name);
        if (hist_it != it->second.hist2D.end()) {
            return hist_it->second.get();
        }
    }
    return nullptr;
}

TH1D* HistogramManager::getBetaExposureHist(int isotopeIdx, int betaTypeIdx) const {
    auto it = histograms.find(HistType::BetaExposure);
    if (it != histograms.end() && 
        isotopeIdx < static_cast<int>(it->second.arrays.size()) &&
        betaTypeIdx < static_cast<int>(it->second.arrays[isotopeIdx].size())) {
        return it->second.arrays[isotopeIdx][betaTypeIdx].get();
    }
    return nullptr;
}

TH1D* HistogramManager::getEventHist(DetType det, int cutIdx) const {
    std::string name;
    switch (det) {
        case DetType::TOF: name = Form("EventEkTOF_%d", cutIdx); break;
        case DetType::NaF: name = Form("EventEkNaF_%d", cutIdx); break;
        case DetType::AGL: name = Form("EventEkAGL_%d", cutIdx); break;
    }
    return getHist1D(HistType::Event, name);
}

void HistogramManager::write(TDirectory* dir) {
    if (!dir) {
        throw std::runtime_error("Invalid directory for writing histograms");
    }

    dir->cd();

    // 写入所有直方图
    for (const auto& [type, collection] : histograms) {
        // 写入1D直方图
        for (const auto& [name, hist] : collection.hist1D) {
            if (hist) hist->Write();
        }
        // 写入2D直方图
        for (const auto& [name, hist] : collection.hist2D) {
            if (hist) hist->Write();
        }
        // 写入数组直方图
        for (const auto& array : collection.arrays) {
            for (const auto& hist : array) {
                if (hist) hist->Write();
            }
        }
    }
}

} // namespace AMS_Iso