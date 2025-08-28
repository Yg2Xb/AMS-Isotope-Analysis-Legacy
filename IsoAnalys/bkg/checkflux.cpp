#include "../Tool.h"  
using namespace AMS_Iso;

void checkflux(){
    double rigidity, flux;
    
    TFile *inputFile = new TFile("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BoronFlux.root", "READ");

    TGraphAsymmErrors *boronFlux = (TGraphAsymmErrors *)inputFile->Get("graph1");
    
    boronFlux->GetPoint(17, rigidity, flux);
    
    double kineticEnergy = rigidityToKineticEnergy(rigidity, 5, 11);
    double dRdEk = dR_dEk(kineticEnergy, 5, 11);
    double convertedFlux = flux * dRdEk;
    double fraction = 0.7;
    double C12Ekflux = convertedFlux*fraction; 
    
    cout << "Rigidity: " << rigidity << " GV" << endl;
    cout << "Boron Rig Flux(Total): " << flux << " m^-2 sr^-1 s^-1 GV^-1" << endl;
    cout << "Kinetic Energy: " << kineticEnergy << " GeV/n" << endl;
    cout << "dR/dEk: " << dRdEk << endl;
    cout << "Converted Boron Ek Flux: " << convertedFlux << " m^-2 sr^-1 s^-1 (GeV/n)^-1" << endl;
    cout << "Converted B11 Ek Flux: " << convertedFlux*fraction << " m^-2 sr^-1 s^-1 (GeV/n)^-1" << endl;

    inputFile->Close();
    delete inputFile;
}