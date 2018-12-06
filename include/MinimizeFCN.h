#ifndef MINIMIZE_FCN
#define MINIMIZE_FCN


#include "Minuit/FCNBase.h"


#include <vector>


class MinimizeFCN : public FCNBase
{

    public:

    MinimizeFCN()
        : _error_def(1.0)
    {
    }


    ~MinimizeFCN()
    {
    }


    virtual double up() const
    {
        return _error_def;
    };

    virtual double operator()(const std::vector<double>& par) const
    {
        assert(par.size() == 1);

        // get parameter
        double epsilon_31 = par.at(0);

        MinimizeFunction(epsilon_31);

        double chi2 = 0.0;

    }

    void setErrorDef(const double def)
    {
        _error_def = def;
    }

    private:

    double _error_def;

}


#endif // MINIMIZE_FCN
