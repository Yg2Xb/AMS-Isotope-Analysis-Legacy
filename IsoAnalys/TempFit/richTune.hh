#ifndef RICHTUNE_HH
#define RICHTUNE_HH

#include <cmath>       // For mathematical functions like sqrt and fabs
#include <TRandom3.h>  // For TRandom3 random number generator
#include <iostream>    // For debugging or any output if needed

// Function to get the smeared RICH beta value.
// iz: particle charge Z
// beta: CIEMAT beta after applying corrections
// seed: a random number for each event (e.g., Run + Event number)
// isNaF: true for NaF, false for AGL
double GetSmearRichBeta(int iz, double beta, int seed, bool isNaF);

// Function to get the RICH width based on particle charge and radiator type.
// iz: particle charge Z
// isNaF: true for NaF, false for AGL
double GetRichWidth(int iz, bool isNaF);

#endif // RICHTUNE_HH