#ifndef AUX_HPP
#define AUX_HPP

#include "TF1.h"
//#include <TF2.h>
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"


void print_histo_text(std::ostream& dst, const std::string& title, const TH1* const histo);

/*
Double_t fit_function_2d(Double_t *_x, Double_t *par);
*/
// TODO: this method for 1d fit


// should the parameter go in the "function" or the "data"?
// currently I have it in the data, not the function, which is weird?

// x[0] stores the current x value, par[0] stores the "parameter" (amplitude)
// par[1 ...] stores the data for the other histogram, however these parameters
// are constants
// NOTE: par also includes the "position of the end of the final bin"
Double_t fit_function(Double_t *x, Double_t *par);


// chi square test histogram to TGraph
Double_t chi_square_test(const TGraph* const graph_fit, const TH1* const histo_data);


// histo_fit is the fit function
// histo_data is the data
// TODO: assumes ranges are the same
Double_t chi_square_test(const TH1* const histo_fit, const TH1* const histo_data, const Double_t min, const Double_t max);


// histo_fit is the fit function
// histo_data is the data
Double_t chi_square_test(const TH1* const histo_fit, const TH1* const histo_data);


// driver for chisquare test
void driver(const TH1* const histo_fit, const TH1* const histo_data, Double_t &param);


// histo_fit is the fit function
// histo_data is the data
// 2d version
Double_t chi_square_test(const TH2* const histo_fit, const TH2* const histo_data);
Double_t chi_square_test(const TH2* const histo_fit, const TH2* const histo_data, const Double_t min, const Double_t max, Int_t &non_empty_bins_2d);


#endif
