#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double MASS_UNIT = 0.931;  // GeV/c^2

double kineticEnergyToBeta(double kineticEnergy) {
    if (kineticEnergy < 0.0) return -9.0;
    
    double gamma = kineticEnergy / MASS_UNIT + 1.0;
    return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

void printBetaBin() {
    const int RIGIDITY_BINS = 73;
    const double KineticEnergyBins[] = {0.08, 0.13, 0.17, 0.21, 0.27, 0.33, 0.41, 0.49, 0.59, 0.70, 0.82, 
        0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16, 2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 
        5.64, 6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51, 11.45, 12.45, 13.50, 14.65, 15.84, 17.14, 
        18.54, 20.04, 21.64, 23.34, 25.19, 27.13, 29.23, 31.48, 33.93, 36.53, 39.33, 42.33, 45.58, 
        49.08, 53.08, 57.08, 61.58, 66.58, 72.57, 79.07, 86.57, 95.07, 104.57, 115.57, 128.57, 
        144.57, 164.07, 188.57, 219.57, 261.57, 329.07, 439.07, 649.07, 1649.07};

    std::cout << "Rbins_beta[" << RIGIDITY_BINS + 1 << "] = {";
    for (int i = 0; i <= RIGIDITY_BINS; i++) {
        if (i % 10 == 0) std::cout << "\n    ";
        std::cout << std::fixed << std::setprecision(8) << kineticEnergyToBeta(KineticEnergyBins[i]);
        if (i != RIGIDITY_BINS) std::cout << ", ";
    }
    std::cout << "\n};\n";

}