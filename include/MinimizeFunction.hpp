#ifndef MINIMIZE_FUNCTION
#define MINIMIZE_FUNCTION




#include <vector>
#include <cassert>


class MinimizeFunction
{

    public:

    MinimizeFunction(const double epsilon_31)
        : epsilon_31(epsilon_31)
    {

        // initialization here
        this->ReadData();


    }

    ~MinimizeFunction()
    {
    }

    // get

    double operator()(const std::vector<double>& parameter)
    {

        assert(par.size() == 1);


    }

    double up()
    {

    }
    

    // others

    void ReadData();


    private:

    double epsilon_31;


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // data histograms
    TH2D *h_nEqNull;
    TH2D *h_nEqTwo;

};



#endif // MINIMIZE_FUNCTION
