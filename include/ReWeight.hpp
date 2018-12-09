#ifndef REWEIGHT_HPP
#define REWEIGHT_HPP


#include "TH1.h"
#include "TH2.h"


#include <cmath>
#include <iostream>


// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug);

// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight2(const Double_t T1, const Double_t T2, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug);

// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon_baseline
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight3(const Double_t T1, const Double_t T2,
                  const Double_t epsilon_baseline, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug);

#endif
