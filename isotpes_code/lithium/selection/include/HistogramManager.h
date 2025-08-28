/***********************************************************
 *  File: HistogramManager.h
 *
 *  Header file for AMS Isotopes histogram management.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <string_view>

#include <string>
#include <vector>
#include <TDirectory.h>
#include <memory>
#include <map>

#include <TH1D.h>
#include <TH2D.h>
#include "basic_var.h"

namespace AMS_Iso {

// 直方图类型枚举
enum class HistType {
    MC,                 // MC总数
    Exposure,           // 曝光时间
    BetaExposure,      // Beta曝光时间
    Event,             // 事例计数
    Resolution,    // RICH分辨率
    Efficiency,         // 效率
    ChargeTempFit 
};

class HistogramManager {
public:
    HistogramManager() = default;
    ~HistogramManager() = default;

    // 创建直方图
    void createHistograms(const IsotopeVar* isotope, int useMass);

    // 访问直方图
    TH1* getHistogram(HistType type, const std::string& name) const;
    TH1D* getHist1D(HistType type, const std::string& name) const;
    TH2D* getHist2D(HistType type, const std::string& name) const;
    TH1D* getBetaExposureHist(int isotopeIdx, int betaTypeIdx) const;
    TH1D* getEventHist(DetType det, int cutIdx) const;

    // 写入直方图
    void write(TDirectory* dir);

private:
    // 直方图集合结构
    struct HistCollections {
        std::map<std::string, std::unique_ptr<TH1D>> hist1D;
        std::map<std::string, std::unique_ptr<TH2D>> hist2D;
        std::vector<std::vector<std::unique_ptr<TH1D>>> arrays;
    };

    // 创建不同类型的直方图
    void createMCHistograms();
    void createExposureHistograms();
    void createEventAndKineticHistograms();
    void createResolutionHistograms();
    void createEfficiencyHistograms();
    void createChargeTempFitHistograms();

    // 成员变量
    std::map<HistType, HistCollections> histograms;
    const IsotopeVar* isotope_ = nullptr;
    int useMass_ = -1; 
};

} // namespace AMS_Iso