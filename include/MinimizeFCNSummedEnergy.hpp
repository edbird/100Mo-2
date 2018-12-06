#ifndef HISTOGRAM_FCN_HPP
#define HISTOGRAM_FCN_HPP


// root headers
#include "TH1.h"


// C++ headers
#include <cassert>


////////////////////////////////////////////////////////////////////////////////
//
// HistogramFCN
//
// This class fits h1 * amplitude to h2
//
////////////////////////////////////////////////////////////////////////////////
class HistogramFCN
{

    public:

    // notes on usage for summed energy histograms:
    // the reweighted histogram is passed into h1
    // the original histogram is passed into h2
    // h1 is fit to h2 using the parameter amplitude,
    // (h1 * amplitude - h2) / error is minimized
    // the errors are obtained from original, h2
    // therefore h1 is the "function" and h2 is the "data"
    HistogramFCN(const TH1* const h1, const TH1* const h2)
        : h1(h1)
        , h2(h2)
        , amplitude(1.0)
    {
        assert(h1 != nullptr);
        assert(h2 != nullptr);
        assert(h1->GetNbinsX() == h2->GetNbinsX());
    }

    ~HistogramFCN()
    {
    }

    virtual double up()
    {
        return error_def;
    }

    virtual double operator()(const std::vector<double>& const param)
    {
        assert(par.size() == 1);
        assert(h1 != nullptr);
        assert(h2 != nullptr);
        assert(h1->GetNbinsX() == h2->GetNbinsX());

        double chi2{0.0};
 
        for(int ix{0}; ix <= h1->GetNbinsX(); ++ ix)
        {
            double c1{h1->GetBinContent(ix)};
            double c2{h2->GetBinContent(ix)};

            double error{c2->GetBinError(ix)};

            double chi{(c1 * amplitude - c2) / error};
            chi2 += (chi * chi);
        }

        return chi2;
    }

    void SetErrorDef(const double error_def)
    {
        this->error_def = error_def;
    }






    private:

    const TH1* h1;
    const TH1* h2;
    double amplitude;
    double error_def;

};



#endif // HISTOGRAM_FCN_HPP
