

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
#include "read_data.hpp"
#include "Flag.hpp"
#include "ReWeight.hpp"
#include "HistogramFCN.hpp"


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

    TFile *f{nullptr};
    TTree *t{nullptr};

    Int_t num_bins{40};

    //std::string filename;
    std::string arg_output_filename; // Remove?

    
    ////////////////////////////////////////////////////////////////////////////
    // CONSTANTS
    ////////////////////////////////////////////////////////////////////////////

    const Double_t bb_Q{3.034}; // Q value of decay, MeV

    const Int_t dimension_xy{1001}; // dimension of 2d decay rate histograms


    ////////////////////////////////////////////////////////////////////////////
    // DATA
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage for raw data
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;
    double psiN0;
    double psiN2;


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
    // reweighted electron energy after reweighting by epsilon (model)
    TH1D *h_el_energy_reweight;
    
    // original summed electron energy (data)
    TH1D *h_el_energy_sum_original;
    // reweigted summed electron energy (model)
    TH1D *h_el_energy_sum_reweight;

    // general weight histogram
    TH2D *h_gen_weight; // TODO: remove?






        /////////////////////////
        // Analysis.ReadData() //
        /////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA FOR PHASE SPACE FACTORS AND DECAY RATE DATA
    ////////////////////////////////////////////////////////////////////////////

    // read
    read_data_helper("nEqNull.dat", "nEqTwo.dat", "psiN0.txt", "psiN2.txt", data_nEqNull, data_nEqTwo, psiN0, psiN2);

    ////////////////////////////////////////////////////////////////////////////
    // TODO: apply phase space factor to data
    ////////////////////////////////////////////////////////////////////////////

    // ...
    // Note: Do not do this

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    convert_data_to_histogram_format(data_nEqNull, data_nEqTwo, h_nEqNull, h_nEqTwo);
    // TODO: move code from this function including alloc. of histograms
    // into main code

    // Print graphs for h_nEqNull, h_nEqTwo
    #if 0
    TCanvas *c_nEqNull;
    c_nEqNull = new TCanvas("c_nEqNull", "", 4000, 3000);
    h_nEqNull->Draw("colz");
    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");

    TCanvas *c_nEqTwo;
    c_nEqTwo = new TCanvas("c_nEqTwo", "", 4000, 3000);
    h_nEqTwo->Draw("colz");
    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    #endif

    std::vector<std::vector<double>> data_0;
    std::vector<std::vector<double>> data_1;
    std::vector<std::vector<double>> data_2;

    ////////////////////////////////////////////////////////////////////////////
    // CREATE INTERMEDIATE DATA
    // Format: Nx3 array, each array a different value of epsilon
    ////////////////////////////////////////////////////////////////////////////

    // create data array for "complete data"
    // (with phase space factors)
    #if 0
    // epsilon_31 = 0.0
    create_data_with_phase_space_factor(data_0, 0.0, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    // epsilon_31 = 0.4
    create_data_with_phase_space_factor(data_1, 0.4, data_nEqNull, data_nEqTwo, psiN0, psiN2);
    
    // epsilon_31 = 0.8
    create_data_with_phase_space_factor(data_2, 0.8, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    TH2D *h_data_0;
    TH2D *h_data_1;
    TH2D *h_data_2;
    h_data_0 = new TH2D("h_data_0", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_data_1 = new TH2D("h_data_1", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_data_2 = new TH2D("h_data_2", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    h_data_0->SetStats(0);
    h_data_1->SetStats(0);
    h_data_2->SetStats(0);
    //h_ratio->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            h_data_0->SetBinContent(i + 1, j + 1, data_0.at(i * dimension_xy + j)[2]);
            h_data_1->SetBinContent(i + 1, j + 1, data_1.at(i * dimension_xy + j)[2]);
            h_data_2->SetBinContent(i + 1, j + 1, data_2.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            // TODO: -1 or -0 ?
            if(i + j < dimension_xy - 1)
            {
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            }
        }
    }

    TCanvas *c_data_0;
    c_data_0 = new TCanvas("c_data_0", "", 4000, 3000);
    h_data_0->Draw("colz");
    c_data_0->SaveAs("c_data_0.png");
    c_data_0->SaveAs("c_data_0.pdf");
    c_data_0->SaveAs("c_data_0.C");

    TCanvas *c_data_1;
    c_data_1 = new TCanvas("c_data_1", "", 4000, 3000);
    h_data_1->Draw("colz");
    c_data_1->SaveAs("c_data_1.png");
    c_data_1->SaveAs("c_data_1.pdf");
    c_data_1->SaveAs("c_data_1.C");

    TCanvas *c_data_2;
    c_data_2 = new TCanvas("c_data_2", "", 4000, 3000);
    h_data_2->Draw("colz");
    c_data_2->SaveAs("c_data_2.png");
    c_data_2->SaveAs("c_data_2.pdf");
    c_data_2->SaveAs("c_data_2.C");
    #endif


        //////////////////////////////////
        // Analysis.InitEventLoopTree() //
        //////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // initialize file and branch addresses
    ////////////////////////////////////////////////////////////////////////////

    // NEMO-3 data/MC read from tree
    // input file

    // read file using program arguments
    f = new TFile(filename.c_str());
    if(f->IsOpen()) std::cout << "OPEN" << std::endl;
    t = (TTree*)f->Get("NewElectronNtuplizer/NewElectronNtuplizer");
    std::cout << "f=" << f << " t=" << t << std::endl;

    t->SetBranchAddress("nElectrons", &nElectrons);
    t->SetBranchAddress("trueT1", &trueT1);
    t->SetBranchAddress("trueT2", &trueT2);
    t->SetBranchAddress("el_energy_", el_energy_);
    // new method requires weight to be saved to tree
    gen_weight = 1.0;
    if(gen_weight_enable == true)
    {
        std::cout << "gen_weight_enable" << std::endl;
        t->SetBranchAddress("gen_weight", &gen_weight);
    }
    else
    {
        std::cout << "gen_weight_disable" << std::endl;
    }

    // default to MC mode
    MODE_FLAG m_mode = MODE_FLAG::MODE_MC;
    
    // switch to enable data mode
    // defaut is MC
    #if 0
        analysis.SetModeFlag(MODE_FLAG::MODE_DATA);
    #endif





        /////////////////////////////////////
        // Analysis.RunOverEpsilonVector() //
        /////////////////////////////////////

        ///////////////////////////////////////
        // Analysis.InitEventLoopHistogram() //
        ///////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////

    h_gen_weight = new TH2D("h_gen_weight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_gen_weight->SetStats(0);


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////
    
    std::string h_name_append;

    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    h_el_energy_original = new TH1D((std::string("h_el_energy_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    //h_el_energy_original->SetStats(0);
    //h_el_energy_original->SetLineColor(2);
    //h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    /*
    h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    */
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    h_el_energy_sum_original = new TH1D((std::string("h_el_energy_sum_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    /*
    h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore
    */


    // these are re-made each iteration, TODO

    h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);

    h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore


    // TODO: I don't remember what to do about systematics, so remove for now

        //////////////////////////
        // Analysis.EventLoop() //
        //////////////////////////

        // baseline (does not depend on epsilon_31 fit parameter)

    std::cout << "Processing data (BASELINE)" << std::endl;
    //const Double epsilon_31{0.0};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        t->GetEntry(ix);
        

        // analysis only valid for 2 electron events
        if(nElectrons != 2) continue;


        // note: no systematic energy shift in this class
        //Double_t el_energy_0{el_energy_[0] * systematic_energy_mult};
        //Double_t el_energy_1{el_energy_[1] * systematic_energy_mult};
        Double_t el_energy_0{el_energy_[0]};
        Double_t el_energy_1{el_energy_[1]}; // * systematic_energy_mult



        
        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        //std::cout << "T1=" << T1 << "T2= " << T2 << std::endl;
        //std::cin.get();

        // if MC apply energy degradation (correction)
        /*
        if(m_mode == MODE_FLAG::MODE_MC)
        {
            Double_t visible_true_ratio_1{g_energy_correction->Eval(T1)};
            Double_t visible_true_ratio_2{g_energy_correction->Eval(T2)};

            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }
        */
        // if MC apply energy degradation (correction)
        if(m_mode == MODE_FLAG::MODE_MC)
        {

            Double_t visible_true_ratio_1{1.0};
            Double_t visible_true_ratio_2{1.0};
            
            /*
            if(b_energy_correction_systematic_enabled == true)
            {
                //std::cout << "energy correction systematic is enabled" << std::endl;

                // if energy correction systematic is enabled, choose
                // energy correction value depending on which subanalysis
                // class this is
                // 0 is default
                // 1 is high systematic
                // -1 is low systematic
                if(m_energy_correction_systematic_mode == 0)
                {
                    visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
                }
                else if(m_energy_correction_systematic_mode == 1)
                {
                    visible_true_ratio_1 = g_energy_correction_systematic_high->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction_systematic_high->Eval(1000.0 * T2);
                }
                else if(m_energy_correction_systematic_mode == -1)
                {
                    visible_true_ratio_1 = g_energy_correction_systematic_low->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction_systematic_low->Eval(1000.0 * T2);
                }
            }
            else
            {
            */
                //std::cout << "energy correction systematic is disabled" << std::endl;

                // if systematics for energy correction are disabled...
                /*
                visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
                visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
                */
            /*
            }
            */

            //std::cout << "visible_true_ratio = " << visible_true_ratio_1 << ", " << visible_true_ratio_2 << std::endl;
            // apply energy correction with systematics if enabled
            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }
        
        /*** SYSTEMATICS **********************************************************/

        /*
        // this if statement sorts out the logical problem of having different
        // high/low sysematic energy multipliers for the purpose of using them
        // as labels to address the SubAnalysis entries in the map inside Analysis,
        // and simultaniously allowing the systematic energy mult systematic to be
        // turned off while another systematic is on
        if(systematic_energy_mult_enable == true)
        {
            el_energy_0 = el_energy_0 * systematic_energy_mult;
            el_energy_1 = el_energy_1 * systematic_energy_mult;
        }
            
        // linear energy offset systematic
        el_energy_0 = el_energy_0 + systematic_energy_offset;
        el_energy_1 = el_energy_1 + systematic_energy_offset;
        */

        // efficiency systematic
        // TODO: can remove, as weight_efficiency = systematic_efficiency
        Double_t weight_efficiency = 1.0;
        //weight_efficiency = weight_efficiency * systematic_efficiency;
        weight_efficiency = weight_efficiency * 1.0;

        // generator weight (MC weight) multiplied by weight efficiency
        Double_t aux_weight{gen_weight};
        aux_weight = aux_weight * weight_efficiency;
        
        // TODO; what happens if the energy shift / systematics move the energy
        // out of a valid range
        // answer: nothing, reweight function depends only on T1 T2
        // TODO: should T1 and T2 be shifted by systematic?


        // TODO: energy degratation systematic

        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) continue;
            if(el_energy_1 < 0.3) continue;
        }
        */


        // NOTE: more logical to set variables
        // weight_1, weight_2
        // for baseline and reweighted (now "baseline" and "test" / "universe")
        // then fill each histogram with each different weight
        // ?
        // NOTE: why would we reweight at all, why not use the decay rates from the
        // theorists directly?

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        //Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

        // original
        h_el_energy_original->Fill(el_energy_0, 1.0 * aux_weight);
        h_el_energy_original->Fill(el_energy_1, 1.0 * aux_weight);
        
        h_el_energy_sum_original->Fill(el_energy_0 + el_energy_1, 1.0 * aux_weight);
        
        /*
        if(el_energy_0 <= el_energy_1)
        {
            h_el_energy_2d_original->Fill(el_energy_0, el_energy_1, 1.0 * aux_weight);
        }
        else
        {
            h_el_energy_2d_original->Fill(el_energy_1, el_energy_0, 1.0 * aux_weight);
        }
        */


        h_gen_weight->Fill(trueT1, trueT2, gen_weight);


    }



        //////////////////////////
        // Analysis.EventLoop() //
        //////////////////////////


    std::cout << "Processing data" << std::endl;
    Long64_t prog_c{-1};
    //const Double epsilon_31{0.0};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        t->GetEntry(ix);


        // analysis only valid for 2 electron events
        if(nElectrons != 2) continue;
        

        Double_t el_energy_0{el_energy_[0]};
        Double_t el_energy_1{el_energy_[1]}; // * systematic_energy_mult



        
        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        // if MC apply energy degradation (correction)
        if(m_mode == MODE_FLAG::MODE_MC)
        {

            Double_t visible_true_ratio_1{1.0};
            Double_t visible_true_ratio_2{1.0};
            
            /*
            if(b_energy_correction_systematic_enabled == true)
            {
                //std::cout << "energy correction systematic is enabled" << std::endl;

                // if energy correction systematic is enabled, choose
                // energy correction value depending on which subanalysis
                // class this is
                // 0 is default
                // 1 is high systematic
                // -1 is low systematic
                if(m_energy_correction_systematic_mode == 0)
                {
                    visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
                }
                else if(m_energy_correction_systematic_mode == 1)
                {
                    visible_true_ratio_1 = g_energy_correction_systematic_high->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction_systematic_high->Eval(1000.0 * T2);
                }
                else if(m_energy_correction_systematic_mode == -1)
                {
                    visible_true_ratio_1 = g_energy_correction_systematic_low->Eval(1000.0 * T1);
                    visible_true_ratio_2 = g_energy_correction_systematic_low->Eval(1000.0 * T2);
                }
            }
            else
            {
            */
                //std::cout << "energy correction systematic is disabled" << std::endl;

                // if systematics for energy correction are disabled...
                /*
                visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
                visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
                */
            /*
            }
            */

            //std::cout << "visible_true_ratio = " << visible_true_ratio_1 << ", " << visible_true_ratio_2 << std::endl;
            // apply energy correction with systematics if enabled
            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }

        
        /*** SYSTEMATICS **********************************************************/

        /*
        // this if statement sorts out the logical problem of having different
        // high/low sysematic energy multipliers for the purpose of using them
        // as labels to address the SubAnalysis entries in the map inside Analysis,
        // and simultaniously allowing the systematic energy mult systematic to be
        // turned off while another systematic is on
        if(systematic_energy_mult_enable == true)
        {
            el_energy_0 = el_energy_0 * systematic_energy_mult;
            el_energy_1 = el_energy_1 * systematic_energy_mult;
        }
            
        // linear energy offset systematic
        el_energy_0 = el_energy_0 + systematic_energy_offset;
        el_energy_1 = el_energy_1 + systematic_energy_offset;
        */

        // efficiency systematic
        // TODO: can remove, as weight_efficiency = systematic_efficiency
        Double_t weight_efficiency = 1.0;
        //weight_efficiency = weight_efficiency * systematic_efficiency;
        weight_efficiency = weight_efficiency * 1.0;

        // generator weight (MC weight) multiplied by weight efficiency
        Double_t aux_weight{gen_weight};
        aux_weight = aux_weight * weight_efficiency;
        
        // TODO; what happens if the energy shift / systematics move the energy
        // out of a valid range
        // answer: nothing, reweight function depends only on T1 T2
        // TODO: should T1 and T2 be shifted by systematic?


        // TODO: energy degratation systematic
        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) return;
            if(el_energy_1 < 0.3) return;
        }
        */


        // NOTE: more logical to set variables
        // weight_1, weight_2
        // for baseline and reweighted (now "baseline" and "test" / "universe")
        // then fill each histogram with each different weight
        // ?
        // NOTE: why would we reweight at all, why not use the decay rates from the
        // theorists directly?

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};


        // reweight
        h_el_energy_reweight->Fill(el_energy_0, weight * aux_weight);
        h_el_energy_reweight->Fill(el_energy_1, weight * aux_weight);

        h_el_energy_sum_reweight->Fill(el_energy_0 + el_energy_1, weight * aux_weight);

        /*
        if(el_energy_0 <= el_energy_1)
        {
            h_el_energy_2d_reweight->Fill(el_energy_0, el_energy_1, weight * aux_weight);
        }
        else
        {
            h_el_energy_2d_reweight->Fill(el_energy_1, el_energy_0, weight * aux_weight);
        }
        */

        // note: no systematic energy shift in this class
        //Double_t el_energy_0{el_energy_[0] * systematic_energy_mult};
        //Double_t el_energy_1{el_energy_[1] * systematic_energy_mult};
        //Double_t el_energy_0{el_energy_[0]};
        //Double_t el_energy_1{el_energy_[1]};
        // NOTE: removed, defined above, need to check


        

        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        //Double_t T1{trueT1 / bb_Q};
        //Double_t T2{trueT2 / bb_Q};
        // NOTE: removed, defined above, need to check

        // if MC apply energy degradation (correction)
        /*
        if(m_mode == MODE_FLAG::MODE_MC)
        {
            Double_t visible_true_ratio_1{g_energy_correction->Eval(T1)};
            Double_t visible_true_ratio_2{g_energy_correction->Eval(T2)};

            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }
        */

        // TODO: not sure if this should be removed
        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) continue;
            if(el_energy_1 < 0.3) continue;
        }
        */


        //h_gen_weight->Fill(trueT1, trueT2, gen_weight);


    }

        // Analysis.PostProcess() //

    //if(_batch_mode_ == false)
    {
        TCanvas* c_gen_weight;
        c_gen_weight = new TCanvas("c_gen_weight", "", 800, 600);
        c_gen_weight->SetRightMargin(0.15);
        h_gen_weight->Draw("colz");
        c_gen_weight->SaveAs("c_gen_weight.C");
        c_gen_weight->SaveAs("c_gen_weight.png");
        c_gen_weight->SaveAs("c_gen_weight.pdf");
        delete c_gen_weight;
    }


    TCanvas *c = new TCanvas("c", "c", 800, 600);
    h_el_energy_sum_original->SetMarkerColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(2);
    h_el_energy_sum_original->SetLineColor(4);
    h_el_energy_sum_reweight->SetLineColor(2);
    h_el_energy_sum_reweight->Draw("e");
    h_el_energy_sum_original->Draw("esame");
    c->SaveAs("c.png");



        ////////////////////////////////
        // Analysis.SummedEnergyFit() //
        ////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////
    // MINIMIZATION TEST
    ////////////////////////////////////////////////////////////////////////////

    // this minimze function fits the first argument * amplitude to the second argument
    // in other words, it minimizes:
    // chi2 = (h1 * amplitude - h2) / e1
    // reweight should be scaled by the amplitude to fit the original
    // the errors are taken from original
    // notes: how it worked before
    // a TF1 (function) was created with 1 free parameter for amplitude and
    // fixed parameters containing the bin content and bin edges for
    // h_el_energy_sum_reweight
    // the TF1 was fit to the histogram h_el_energy_sum_original, and therefore
    // the errors were obtained from h_el_energy_sum_original
    /*MinimizeFunctionSummedEnergy*/ HistogramFCN theFCN_amplitude(h_el_energy_sum_reweight, h_el_energy_sum_original);
    
    // initial parameters, uncertainties
    std::vector<double> init_par;
    init_par.push_back(1.0); // amplitude, start at 1.0
    std::vector<double> init_err;
    init_err.push_back(0.1);

    // create minimizer
    ROOT::Minuit2::VariableMetricMinimizer theMinimizer;

    // minimize
    ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN_amplitude, init_par, init_err);
    std::cout << "minimum: " << FCN_min << std::endl;

    // TODO: must be a better method of getting chisquare
    std::cout << FCN_min.Parameters().Vec() << std::endl;
    std::vector<double> end_par;
    const double* vec_data{FCN_min.Parameters().Vec().Data()};
    for(int c{0}; c < FCN_min.Parameters().Vec().size(); ++ c)
    {
        end_par.push_back(vec_data[c]);
    }
    std::cout << "chi2=" << theFCN_amplitude.operator()(end_par) << std::endl;
    double amplitude{end_par.at(0)};


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


    // scale the red histogram using the amplitude parameter
    h_el_energy_sum_reweight->Scale(amplitude);
    h_el_energy_reweight->Scale(amplitude);
    //h_el_energy_2d_reweight->Scale(amplitude);


    //  -> SensitivityMeasurementChisquare1
    //  -> ChiSquare_BaselineNoSystematic
    //  -> PrintOutputToFile
    //  -> MakeSensitivityCanvas
    //  -> MakeChiSquareType1
    
    

}
