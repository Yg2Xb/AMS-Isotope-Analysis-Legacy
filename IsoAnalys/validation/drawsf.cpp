#include <TCanvas.h>
#include <TGraph.h>
#include <cmath>
#include <iostream>

constexpr double Mass_Unit = 0.931;
constexpr double Charge = 4;  // charge of Beryllium (Be4)
constexpr double Mass_9 = 9;
constexpr double Mass_10 = 10;
constexpr double SafetyFactor = 1.0; // Default Safety factor if needed

// Function to calculate the safety factor (sf)
double calculateSF(double rigidity) {
    double denom = std::sqrt(std::pow(9 * Mass_Unit, 2) + rigidity * rigidity * 16);  // numerator for sf
    double num = std::sqrt(std::pow(10 * Mass_Unit, 2) + rigidity * rigidity * 16);  // denominator for sf
    return std::sqrt(num / denom);
}

void drawsf() {
    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Safety Factor vs Rigidity", 800, 600);
    
    // Define the range of rigidity
    int nPoints = 100;
    double rigidityMin = 1.0;  // Minimum rigidity
    double rigidityMax = 30.0;  // Maximum rigidity
    double rigidityStep = (rigidityMax - rigidityMin) / (nPoints - 1);
    
    // Create arrays to store the rigidity values and the corresponding sf values
    double rigidity[nPoints], sf[nPoints];
    
    // Fill the arrays
    for (int i = 0; i < nPoints; ++i) {
        rigidity[i] = rigidityMin + i * rigidityStep;
        sf[i] = calculateSF(rigidity[i]);
    }

    // Create a graph from the rigidity and sf data
    TGraph* graph = new TGraph(nPoints, rigidity, sf);
    graph->SetTitle("Safety Factor vs Rigidity");
    graph->GetXaxis()->SetTitle("Rigidity (GV)");
    graph->GetYaxis()->SetTitle("Safety Factor (sf)");

    // Set graph style
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetMarkerSize(1.2);
    
    // Draw the graph
    graph->Draw("ALP");

    // Save the plot to a file (optional)
    canvas->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/Validation/safety_factor_vs_rigidity.png");

    // Keep the canvas open
    canvas->Update();
}

