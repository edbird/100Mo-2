#ifndef STATICS_GROUP
#define STATICS_GROUP


// Local headers
#include "read_data.hpp"


// Root headers
#include "TFile.h"
#include "TTree.h"





////////////////////////////////////////////////////////////////////////////////
//
// StaticsGroup
//
// Put all "static" data in this class (data which does not change each loop,
// does not depend on epsilon_31 parameter
//
////////////////////////////////////////////////////////////////////////////////
class StaticsGroup
{


    public:

    ////////////////////////////////////////////////////////////////////////////
    // class constructor functions
    ////////////////////////////////////////////////////////////////////////////

    StaticsGroup(const double epsilon_31_baseline, const std::string& filename);

    virtual
    ~StaticsGroup();

    ////////////////////////////////////////////////////////////////////////////
    // loop functions
    ////////////////////////////////////////////////////////////////////////////

    // print c_nEqNull, c_nEqTwo to canvas
    void
    PrintInputData();

    // removed
    void
    CreatePrintIntermediateData();

    // event loop
    void
    EventLoop();

    // statistics robustness test - apply statistics randomization
    void
    ApplyStatisticsRandomization
    (TH1D * const h_el_energy_reweight, TH1D * const h_el_energy_sum_reweight) const;

    ////////////////////////////////////////////////////////////////////////////
    // get / set flag / variable functions
    // (call before event loop)
    ////////////////////////////////////////////////////////////////////////////

    TH1D *const
    Get_h_el_energy_original() const;

    TH1D *const
    Get_h_el_energy_sum_original() const;

    double
    Get_epsilon_31_baseline() const;

    Int_t
    Get_num_bins() const;

    TTree*
    Get_tree() const;

    const
    Int_t& Get_nElectrons() const;

    const
    Double_t& Get_trueT1() const;

    const
    Double_t& Get_trueT2() const;

    const
    Double_t* const Get_el_energy_() const;

    const
    Double_t& Get_gen_weight() const;

    const
    Double_t& Get_bb_Q() const;

    const
    double& Get_psiN0() const;

    const double&
    Get_psiN2() const;

    const TH2D* const
    Get_h_nEqNull() const;

    const TH2D* const
    Get_h_nEqTwo() const;

    ////////////////////////////////////////////////////////////////////////////
    // GET / SET PROGRAM BEHAVIOUR FLAGS / VARIABLES
    ////////////////////////////////////////////////////////////////////////////

    void
    SetStatisticsRobustnessTestEnable(const bool enable);
    
    bool
    GetStatisticsRobustnessTestEnable() const;

    void
    SetStatisticsRobustnessTestIterations(const int iterations);

    int
    GetStatisticsRobustnessTestIterations() const;

    void
    SetStatisticsRobustnessTestCurrentIteration(const int current);

    int
    GetStatisticsRobustnessTestCurrentIteration() const;

    void
    StoreChi2Result(const double chi2);

    void
    StoreEpsilon31Result(const double epsilon31);

    void
    StoreEpsilon31ErrResult(const double epsilon31err);

    void
    PrintChi2Results() const;

    void
    PrintEpsilon31Results() const;

    protected:

    // epsilon_31 is not a parameter here as it is a parameter managed
    // by Minuit2 for minimization

    // epsilon_31_baseline is a parameter here
    double epsilon_31_baseline;

    // TFile and TTree for reading NEMO3 MC/data
    std::string filename;
    TFile *f;
    TTree *t;


    //std::string arg_output_filename;

    
    ////////////////////////////////////////////////////////////////////////////
    // CONSTANTS
    ////////////////////////////////////////////////////////////////////////////

    // Q value of decay, MeV
    const Double_t bb_Q;
    
    // number of bins in analysis, 100 keV
    const Int_t num_bins;

    // dimension of 2d decay rate histograms
    const Int_t dimension_xy;


    ////////////////////////////////////////////////////////////////////////////
    // DATA
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage for raw data
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;
    double psiN0;
    double psiN2;
    
    // intermediate data


    ////////////////////////////////////////////////////////////////////////////
    // EVENT LOOP DATA
    ////////////////////////////////////////////////////////////////////////////

    Int_t nElectrons;
    Double_t trueT1;
    Double_t trueT2;
    Double_t el_energy_[2];
    Double_t gen_weight;
    
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // data histograms
    TH2D *h_nEqNull;
    TH2D *h_nEqTwo;
    
    // original electron energy (data)
    TH1D *h_el_energy_original;
    
    // original summed electron energy (data)
    TH1D *h_el_energy_sum_original;

    // general weight histogram
    TH2D *h_gen_weight; // TODO: remove?


    ////////////////////////////////////////////////////////////////////////////
    // PROGRAM BEHAVIOUR FLAGS
    ////////////////////////////////////////////////////////////////////////////

    // flag to enable statistics robustness test
    // this is a diagnostics test to ensure the robustness of the fitting
    // procedure to statistical fluctuations
    // procedure details:
    // take reweight histograms, and apply statistical fluctuation
    // this means: for each bin: draw a random number from a Gaussian
    // with mean of the number of events in that bin, width of the error
    // on that bin
    // use this number as the new number of events (weight) in that bin
    bool statistics_robustness_test_enable;
    int statistics_robustness_test_iterations;
    int statistics_robustness_test_current_iteration;
    std::vector<double> statistics_robustness_test_chi2;
    std::vector<double> statistics_robustness_test_epsilon_31;
    std::vector<double> statistics_robustness_test_epsilon_31_err;
    //TH1D *h_robustness_test_chi2;

};



#endif // STATICS_GROUP
