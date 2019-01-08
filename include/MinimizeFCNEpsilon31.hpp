#ifndef MINIMIZE_FCN_EPSILON_31_HPP
#define MINIMIZE_FCN_EPSILON_31_HPP


// Local headers
#include "StaticsGroup.hpp"
#include "ReWeight.hpp"
#include "MinimizeFCNSummedEnergy.hpp"
#include "aux.hpp"


// Root headers, minuit
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"


// Root headers, other
//#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"


// C++ headers
#include <cassert>
#include <string>
#include <sstream>





////////////////////////////////////////////////////////////////////////////////
//
// MinimizeFCNEpsilon31
//
// This class fits epsilon31 parameter (xi)
//
////////////////////////////////////////////////////////////////////////////////
class MinimizeFCNEpsilon31 : public ROOT::Minuit2::FCNBase
{


    public:

        MinimizeFCNEpsilon31(const StaticsGroup& staticsgroup)
            : error_def{1.0}
            , staticsgroup{staticsgroup}
            , print_iteration_enable{true}
        {
            // if statistics robustness test is enabled, create a histogram
            // with gaussian distributed numbers to use as parameters for
            // applying randomness to reweight histograms
            // (parameters for creating pseudodata)
            if(staticsgroup.GetStatisticsRobustnessTestEnable() == true)
            {
                CreateRandomizedStatistics();
            }

        }


        virtual
        ~MinimizeFCNEpsilon31()
        {
            std::cout << "~MinimizeFCNEpsilon31" << std::endl;
        }


        virtual double Up() const
        {
            return error_def;
        }

        virtual double operator()(const std::vector<double>& param) const
        {

            assert(param.size() == 1);

            const double epsilon_31{param.at(0)};

            //std::cout << "operator() : epsilon_31=" << epsilon_31 << std::endl;
            //std::cin.get();


            ////////////////////////////////////////////////////////////////////
            // reweight histograms have scope limited to this function
            ////////////////////////////////////////////////////////////////////

            // reweighted electron energy after reweighting by epsilon (model)
            TH1D *h_el_energy_reweight{nullptr};
            // reweigted summed electron energy (model)
            TH1D *h_el_energy_sum_reweight{nullptr};

            ////////////////////////////////////////////////////////////////////
            // reweight pseudodata histograms, used for statistics robustness
            // testing
            ////////////////////////////////////////////////////////////////////

            // reweighted electron energy after reweighting by epsilon (model) - pseudodata version
            TH1D *h_el_energy_reweight_pseudo{nullptr};
            // reweigted summed electron energy (model) - pseudodata version
            TH1D *h_el_energy_sum_reweight_pseudo{nullptr};

            // allocate
            Int_t num_bins{staticsgroup.Get_num_bins()};
            std::string h_name_append;

            h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
            h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
            
            // get histograms stored in StaticsGroup
            const TH1D *h_el_energy_sum_original{staticsgroup.Get_h_el_energy_sum_original()};
            /*const*/ TH1D *h_el_energy_original{staticsgroup.Get_h_el_energy_original()};
            
            // fill histograms
            EventLoop(epsilon_31, h_el_energy_reweight, h_el_energy_sum_reweight, h_el_energy_sum_original);

            // if statistics robustness test is enabled, create fluctuation data
            // apply statistics robustness test if enabled
            if(staticsgroup.GetStatisticsRobustnessTestEnable() == true)
            {
                
                h_el_energy_reweight_pseudo = new TH1D((std::string("h_el_energy_reweight_pseduo") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
                h_el_energy_sum_reweight_pseudo = new TH1D((std::string("h_el_energy_sum_reweight_pseudo") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);

                ApplyRandomizedStatistics(h_el_energy_reweight, h_el_energy_sum_reweight, h_el_energy_reweight_pseudo, h_el_energy_sum_reweight_pseudo);
            }
            else
            {
                // if statistics robustness test not enabled, then "copy"
                // randomized reweighted histograms to regular reweighted
                // histograms
                // output = input
                h_el_energy_reweight_pseudo = h_el_energy_reweight;
                h_el_energy_sum_reweight_pseudo = h_el_energy_sum_reweight;
            }
            

            // fit histograms (summed energy)
            double amplitude{1.0};
            Fit(h_el_energy_sum_reweight_pseudo, h_el_energy_sum_original, amplitude);

            // scale
            // scale the red histogram using the amplitude parameter
            h_el_energy_sum_reweight_pseudo->Scale(amplitude);
            h_el_energy_reweight_pseudo->Scale(amplitude);
            //h_el_energy_2d_reweight->Scale(amplitude);

            // calculate chi2
            double chi2{0.0};

            // Note: need to change to use systematics chisquare measurement
            // implemented in 100-Mo (original code)
            chi2 = chi_square_test(h_el_energy_reweight_pseudo, h_el_energy_original);
            // TODO: write chi2 to output and params [done]
            // TODO: subrange fit and el_energy cut
            

            // save iteration output, also saves params
            if(print_iteration_enable)
            {
                PrintIteration(h_el_energy_original, h_el_energy_reweight_pseudo, epsilon_31, chi2, amplitude);
            }

            // read only object, must manage all memory here
            delete h_el_energy_reweight;
            h_el_energy_reweight = nullptr;
            delete h_el_energy_sum_reweight;
            h_el_energy_sum_reweight = nullptr;


            if(staticsgroup.GetStatisticsRobustnessTestEnable() == true)
            {
                // does not do anything if pointer set to null if statistics
                // robustness test disabled
                delete h_el_energy_reweight_pseudo;
                h_el_energy_reweight_pseudo = nullptr;
                delete h_el_energy_sum_reweight_pseudo;
                h_el_energy_sum_reweight_pseudo = nullptr;
            }


            // increment the call counter
            ++ static_iteration_counter;

            // set last chi2
            last_chi2 = chi2;

            return chi2;

        }


        // fill reweight histograms
        void
        EventLoop(const double epsilon_31,
                  TH1D *h_el_energy_reweight,
                  TH1D *h_el_energy_sum_reweight,
                  const TH1D *h_el_energy_sum_original) const;


        // fill vector with random gaussian parameters
        void
        CreateRandomizedStatistics();

        // print random gaussian parameters in histogram
        void
        PrintRandomizedStatistics() const;


        // create pseudodata histograms from reweighted histograms and
        // random gaussian parameters
        void
        ApplyRandomizedStatistics(TH1D * const h_el_energy_reweight,
                                  TH1D * const h_el_energy_sum_reweight,
                                  TH1D * const h_el_energy_reweight_pseduo,
                                  TH1D * const h_el_energy_sum_reweight_pseudo) const;

        
        // fit amplitude of summed electron energy histogram (param A)
        void
        Fit(/*const double epsilon_31,*/
            TH1D *h_el_energy_sum_reweight,
            const TH1D *h_el_energy_sum_original,
            double &amplitude) const;


        void
        PrintIteration(TH1* const h_el_energy_original,
                       TH1* const h_el_energy_reweight,
                       const double epsilon_31,
                       const double chi2,
                       const double amplitude) const;

        
        void
        SetPrintIterationEnable(const bool enable);


        double
        GetLastChi2() const
        {
            return last_chi2;
        }

        /*
        double
        GetLastEpsilon31() const
        {
            return last_epsilon_31;
        }
        */

        static void
        ResetStaticIterationCounter()
        {
            static_iteration_counter = 0;
        }


    protected:

        static unsigned long long static_iteration_counter;
        bool print_iteration_enable;

        double error_def;

        const StaticsGroup& staticsgroup;

        // store Gaussian distribution generated parameters for creating
        // fluctuated pseudodata
        std::vector<double> gen_gaussian_params;

        // last obtained chi2 result
        mutable double last_chi2;

};


#endif // MINIMIZE_FCN_EPSILON_31_HPP
