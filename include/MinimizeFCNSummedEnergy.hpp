#ifndef MINIMIZE_FCN_SUMMED_ENERGY_HPP
#define MINIMIZE_FCN_SUMMED_ENERGY_HPP


// Root headers, minuit
#include "Minuit2/FCNBase.h"


// Root headers, other
#include "TH1.h"


// C++ headers
#include <cassert>
#include <iostream>



////////////////////////////////////////////////////////////////////////////////
//
// MinimizeFCNSummedEnergy
//
// This class fits h1 * amplitude to h2
//
////////////////////////////////////////////////////////////////////////////////
class MinimizeFCNSummedEnergy : public ROOT::Minuit2::FCNBase
{

    public:

    // notes on usage for summed energy histograms:
    // the reweighted histogram is passed into h1
    // the original histogram is passed into h2
    // h1 is fit to h2 using the parameter amplitude,
    // (h1 * amplitude - h2) / error is minimized
    // the errors are obtained from original, h2
    // therefore h1 is the "function" and h2 is the "data"
    MinimizeFCNSummedEnergy(const TH1* const h1, const TH1* const h2)
        : h1(h1)
        , h2(h2)
        , error_def{1.0}
    {
        assert(h1 != nullptr);
        assert(h2 != nullptr);
        assert(h1->GetNbinsX() == h2->GetNbinsX());
    }

    virtual
    ~MinimizeFCNSummedEnergy()
    {
    }

    virtual double Up() const
    {
        return error_def;
    }

    virtual double operator()(const std::vector<double>& param) const
    {
        assert(param.size() == 1);
        assert(h1 != nullptr);
        assert(h2 != nullptr);
        assert(h1->GetNbinsX() == h2->GetNbinsX());

        double amplitude{param.at(0)};

        double chi2{0.0};
 
        for(int ix{1}; ix <= h1->GetNbinsX(); ++ ix)
        {
            double c1{h1->GetBinContent(ix)};
            double c2{h2->GetBinContent(ix)};

            double error{h2->GetBinError(ix)};

            if(error == 0.0)
            {
                if(c2 == 0.0)
                {
                    continue;
                }
                else
                {
                    std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
                }
            }

            double chi{(c1 * amplitude - c2) / error};
            chi2 += (chi * chi);
        }

        //std::cout << "chi2=" << chi2 << std::endl;
        return chi2;
    }

    void SetErrorDef(const double error_def)
    {
        this->error_def = error_def;
    }






    private:

    const TH1* h1;
    const TH1* h2;
    double error_def;

};



#endif // MINIMIZE_FCN_SUMMED_ENERGY_HPP
