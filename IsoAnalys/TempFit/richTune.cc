#include "richTune.hh"

// Function to get the smeared RICH beta value
double GetSmearRichBeta(int iz, double beta, bool isNaF) {
    if (beta <= 0) return beta;
    
    // Inverse the beta as per the given logic
    beta = 1 / beta;
    
    // Get the width based on charge and radiator type
    double wid = GetRichWidth(iz, isNaF);
    
    // Initialize random number generator with the provided seed
    std::random_device rd;  // 使用随机设备生成种子
    TRandom3 rand(rd()); 
    
    // Smear the beta using a Gaussian distribution
    double smear = rand.Gaus();
    double newbeta = beta + wid * smear;
    newbeta = 1 / newbeta;  // Invere again to get the smeared beta
    
    return newbeta;
}

// Function to get the RICH width based on particle charge and radiator type
double GetRichWidth(int iz, bool isNaF) {
    const int nch = 7;
    int zch[nch] = {2, 3, 4, 5, 6, 7, 8};  // Z values for different charge states
    
    const int nrad = 2;
    double datasig[nrad][nch]={//%core sigma after fixing core fraction to 80%
        0.0021171,  0.00168615, 0.00129161, 0.00108304, 0.00100818,0.000897614,0.000859017,
        0.000679357,0.000509869,0.000421296,0.000379986,0.00034206,0.000328994,0.000317768,
    };
    double datamcdiff[nrad][nch]={
        1.12,1.14,1.14,1.14,1.21,1.21,1.21,
        1.08,1.13,1.13,1.13,1.13,1.13,1.13,
    };
    
    // Determine if we are using NaF or AGL
    int irad = (isNaF) ? 0 : 1;

    // Find the index corresponding to the given charge Z
    int ich = 0;
    if (iz < zch[0]) {
        ich = 0;
    } else if (iz > zch[nch - 1]) {
        ich = nch - 1;
    } else {
        for (int i = 0; i < nch; i++) {
            if (zch[i] == iz) {
                ich = i;
                break;
            }
        }
    }
    
    // Retrieve the data sigma and MC sigma
    double dsig = datasig[irad][ich];
    double msig = datasig[irad][ich] / datamcdiff[irad][ich];

    double datasigBer[nrad]={//%core sigma after fixing core fraction to 90%
        0.001425, 0.000455
    }; 
    double datamcdiffBer[nrad][3]={
        1.1709120789,1.1632653061, 1.1623164763,
        1.1375,      1.1346633416, 1.131840796
    };
    /*
    if(iz == 4)
    {
        dsig = datasigBer[irad];
        msig = datasigBer[irad] / datamcdiffBer[irad][ia];
    }
    */
    // Calculate the width (wch) based on the difference between data and MC sigmas
    double wch = 0;
    if (dsig > msig) {
        wch = std::sqrt(std::fabs(std::pow(dsig, 2) - std::pow(msig, 2)));
    }
    
    return wch;
}