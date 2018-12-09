

// C++ headers
#include <iostream>


// Root headers, minuit
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"


// Root headers, other
#include "TFile.h"
#include "TTree.h"

#include "TH1.h"
#include "TH2.h"

#include "TCanvas.h"


// local headers, external
#include "program_arguments.hpp"


// local headers, internal
#include "StaticsGroup.hpp"
#include "read_data.hpp"
#include "Flag.hpp"
#include "ReWeight.hpp"
#include "MinimizeFCNEpsilon31.hpp"


int main(int argc, char* argv[])
{

    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    // new program arguments processing method
    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    //pa.Add("filename", "--filename", "NewElectronNtuplizerExe_Int_ManDB_output.root");
    pa.Add("filename", "--filename", "NewElectronNtuplizerExe_Int_ManDB_output_2e.root");
    //pa.Add("filename", "--filename", "/media/ramdisk/NewElectronNtuplizerExe_Int_ManDB_output.root");
    pa.Add("epsilon-baseline", "--epsilon-baseline", "0.368");
    //pa.Add("batch_mode", "--batch-mode", "false");
    //pa.Add("300 keV energy cut", "--energy-cut", "false");
    pa.Add("energy_cut", "--energy-cut", "false");
    pa.Add("fit_subrange", "--fit-subrange", "false");
    //pa.Add("log_mode", "--log-mode", "false");
    //pa.Add("output_filename", "--output-file", "of_data.txt");
    pa.Add("syst_energy_mult", "--systematic-energy-mult", "1.0");
    pa.Add("syst_energy_add", "--systematic-energy-add", "0.1");
    // TODO: systematics
    
    
    pa.Print();
    pa.Process(argc, argv);


    // this argument is a string
    std::string help{pa.Get("help")};
    // this argument is a string
    std::string filename{pa.Get("filename")};
    std::cout << "filename=" << filename << std::endl;
    // this argument is to be converted to a double
    std::string arg_epsilon_31{pa.Get("epsilon-baseline")};
    
    // this argument is to be converted to a bool
    std::string arg_energy_cut{pa.Get("energy_cut")};
    std::string arg_fit_subrange{pa.Get("fit_subrange")};
    //std::string arg_log_mode{pa.Get("log_mode")};
    //std::string arg_output_filename{pa.Get("output_filename")};
    std::string arg_systematic_energy_mult{pa.Get("syst_energy_mult")};
    std::string arg_systematic_energy_add{pa.Get("syst_energy_add")};
    // TODO: systematics
    
    
    // process gathered argument data
    bool gen_weight_enable{false};
    if(filename != pa.GetDefault("filename"))
    {
        gen_weight_enable = true;
        std::cout << "gen_weight_enable" << std::endl;
    }
    else
    {
        std::cout << "gen_weight_enable = false" << std::endl;
    }
    // TODO: check filename exists!


    // List of primary functions of interest called for analysis:
    //
    // (Analysis.)
    // Analysis [constructor]
    // ReadData
    // InitEventLoopTree
    // RunOverEpsilonVector:
    //  -> InitEventLoopHistogram
    //  -> EventLoop
    //  -> PostProcess
    //  -> SummedEnergyFit
    //  -> SensitivityMeasurementChisquare1
    //  -> ChiSquare_BaselineNoSystematic
    //  -> PrintOutputToFile
    //  -> MakeSensitivityCanvas
    //  -> MakeChiSquareType1


        /////////////////////////
        // Analysis.Analysis() //
        /////////////////////////

    //double epsilon_31{0.368};
    double epsilon_31{0.8};
    double systematic_energy_mult{1.0}; // TODO: others




    

    const double epsilon_31_baseline{0.368};
    StaticsGroup staticsgroup(epsilon_31_baseline, filename);


    ////////////////////////////////////////////////////////////////////////////
    // MINIMIZATION EPSILON_31 
    ////////////////////////////////////////////////////////////////////////////

    MinimizeFCNEpsilon31 theFCN(staticsgroup);
    
    // initial parameters, uncertainties
    std::vector<double> init_par;
    //init_par.push_back(epsilon_31_baseline);
    init_par.push_back(epsilon_31); // temp changed to see if histo output changes as expected
    // TODO: print histos each itter
    std::vector<double> init_err;
    init_err.push_back(0.1);

    // create minimizer
    ROOT::Minuit2::VariableMetricMinimizer theMinimizer;

    // minimize
    ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, init_par, init_err);
    std::cout << "minimum: " << FCN_min << std::endl;

    // TODO: must be a better method of getting chisquare
    std::cout << FCN_min.Parameters().Vec() << std::endl;
    std::vector<double> end_par;
    const double* vec_data{FCN_min.Parameters().Vec().Data()};
    for(int c{0}; c < FCN_min.Parameters().Vec().size(); ++ c)
    {
        end_par.push_back(vec_data[c]);
    }
    std::cout << "chi2=" << theFCN.operator()(end_par) << std::endl;
    double epsilon_31_best_fit{end_par.at(0)};


    #if 0
    TH1D *h_el_energy_sum_reweight_scale = new TH1D("h_el_energy_sum_reweight_scale", "h_el_energy_sum_reweight_scale", num_bins, 0.0, 4.0);
    for(int ix{1}; ix <= h_el_energy_sum_reweight_scale->GetNbinsX(); ++ ix)
    {
        Double_t amplitude{end_par.at(0)};
        Double_t c{h_el_energy_sum_reweight->GetBinContent(ix) * amplitude};
        Double_t e{h_el_energy_sum_reweight->GetBinError(ix) * amplitude}; // TODO: scale by amplitude?
        h_el_energy_sum_reweight_scale->SetBinContent(ix, c);
        h_el_energy_sum_reweight_scale->SetBinError(ix, e);
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    h_el_energy_sum_original->SetMarkerColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(2);
    h_el_energy_sum_reweight_scale->SetMarkerColor(6);
    h_el_energy_sum_original->SetLineColor(4);
    h_el_energy_sum_reweight->SetLineColor(2);
    h_el_energy_sum_reweight_scale->SetLineColor(6);
    h_el_energy_sum_reweight->Draw("e");
    h_el_energy_sum_original->Draw("esame");
    h_el_energy_sum_reweight_scale->Draw("esame");
    c2->SaveAs("c2.png");
    #endif


    ////////////////////////////////////////////////////////////////////////////
    // MINIMIZATION
    ////////////////////////////////////////////////////////////////////////////

    /*
    MinimizeFunction theFCN(params);

    // initial parameters, uncertainties
    std::vector<double> init_par;
    init_par.push_back(0.368);
    std::vector<double> init_err;
    init_err.push_back(0.1);

    // create minimizer
    VariableMetricMinimizer theMinimizer;

    // minimize
    FunctionMinimum FCN_min = theMinimizer.minimize(theFCN, init_par, init_err);
    std::cout << "minimum: " << FCN_min << std::endl;
    */







#if 0

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    h_el_energy_sum_original->SetMarkerColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(2);
    h_el_energy_sum_original->SetLineColor(4);
    h_el_energy_sum_reweight->SetLineColor(2);
    h_el_energy_sum_reweight->Draw("e");
    h_el_energy_sum_original->Draw("esame");
    c->SaveAs("c.png");
#endif




    //  -> SensitivityMeasurementChisquare1
    //  -> ChiSquare_BaselineNoSystematic
    //  -> PrintOutputToFile
    //  -> MakeSensitivityCanvas
    //  -> MakeChiSquareType1
    
    

}
