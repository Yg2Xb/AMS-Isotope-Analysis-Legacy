#include <iostream>
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "../Tool.h"
using namespace AMS_Iso;
void plot_newbeta() {
    // 定义参数范围
    std::vector<double> alpha_values;
    for (double alpha = 1.; alpha <= 1.005; alpha += 0.001) {
        alpha_values.push_back(alpha);
    }

    int nBeta = 1000;
    double beta_min = 0.65, beta_max = 0.9;

    TCanvas* c1 = new TCanvas("c1", "NewBeta vs Beta", 1600, 1200);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.65, 0.15, 0.88, 0.45);
    int c = 0;
    for (double alpha : alpha_values) {

        std::vector<double> beta_values, newbeta_values;
        
        for (int i = 0; i <= nBeta; ++i) {
            double beta = beta_min + i * (beta_max - beta_min) / nBeta;
            double newbeta = alpha * beta / std::sqrt(1 - beta * beta + std::pow(alpha * beta, 2));
            //double newbeta = 1./(alpha+1./beta);
            double gamma = 1./sqrt(1-beta*beta);
            double newgamma = 1./sqrt(1-newbeta*newbeta);
            double rig = 7*0.931*beta*gamma / 4;
            double Rmass = newbeta*newgamma / (rig*4);
            cout<<1/newbeta<<", "<<1/beta<<endl;

            double diff = abs(Rmass - 1/(7*0.931))/abs(1/(7*0.931));
            cout<<1/newbeta - 1/beta<<endl;

            double ek  = betaToKineticEnergy(beta);

            beta_values.push_back(ek);
            newbeta_values.push_back(1/newbeta - 1/beta);
        }
        c1->SetLeftMargin(0.2);
        TGraph* graph = new TGraph(nBeta + 1, beta_values.data(), newbeta_values.data());
        graph->GetYaxis()->SetTitleOffset(1.5);
        graph->SetLineWidth(2);
        graph->SetLineColor(colors[c]);
        graph->SetTitle(Form("#alpha = %.3f", alpha));

        mg->Add(graph);
        leg->AddEntry(graph, Form("#alpha = %.3f", alpha), "l");
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        c++;
    }

    //mg->GetYaxis()->SetTitleOffset(1.25);
    mg->Draw("AL");
    mg->SetTitle("; Ek/n; 1/tuned_beta - 1/beta");
    leg->Draw();

    // 保存 PNG
    c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/AlphaTuning/betadiff.png");

    delete c1;
}
