//this is estTempFit.cpp
#include "TemplateFitter.cpp"

void estTempFit() {
    // 配置1：1/mass
    FitConfig config1 = {
        "1/Mass",    // 变量名
        0.065,        // 最小值
        0.21,        // 最大值
        0.08, 0.12,
        0.64, 0.88,
        1,           // 直方图后缀
        "RMassTempFit"  // 输出文件前缀
    };

    // 配置2：Estimator
    FitConfig config2 = {
        "Estimator",  // 变量名
        0.1,         // 最小值
        0.64,          // 最大值
        0.1, 0.28,
        0.3, 0.7,
        3,            // 直方图后缀
        "EstTempFit"  // 输出文件前缀
    };

    // 选择使用哪个配置
    TemplateFitter fitter(config1);  
    fitter.RunFitting();
}