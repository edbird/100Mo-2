#include "StaticsGroup.hpp"

#include <algorithm>

#include "TCanvas.h"

StaticsGroup::StaticsGroup(const double epsilon_31_baseline, const std::string& filename)
    : epsilon_31_baseline{epsilon_31_baseline} // TODO: pass to ReWeight2
    , filename{filename}
    , f{nullptr}
    , t{nullptr}
    , num_bins{40}
    , bb_Q{3.034}
    , dimension_xy{1001}
    // event loop data
    , gen_weight{1.0}
    // event loop data / histograms
    , h_nEqNull{nullptr}
    , h_nEqTwo{nullptr}
    , h_el_energy_original{nullptr}
    , h_el_energy_sum_original{nullptr}
    , h_gen_weight{nullptr}
    // optical correction
    , optical_correction_enable{true}
    , f_optical_correction{nullptr}
    , g_optical_correction{nullptr}
    , g_optical_correction_systematic_high{nullptr}
    , g_optical_correction_systematic_low{nullptr}
    , systematic_enable_optical_correction_statistical{false}
    , systematic_optical_correction_statistical_direction{0}
    , systematic_enable_optical_correction_Bi207_EC_peak{false}
    , optical_correction_parameter_a{1.0}
    , optical_correction_parameter_b{0.0}
    , optical_correction_Bi207_EC_1{481.7}
    , optical_correction_Bi207_EC_2{975.7}
    , optical_correction_measured_Bi207_EC_1{481.7}
    , optical_correction_measured_Bi207_EC_2{975.7}
    // statistical robustness test
    , statistics_robustness_test_enable{false}
    , statistics_robustness_test_iterations{1}
    , statistics_robustness_test_current_iteration{0}
    // mc/data flag
    , mc_flag{false}
    // other systematics
    , systematic_enable_energy_multiply{false}
    , systematic_energy_multiply{1.0}
    , systematic_enable_energy_add{false}
    , systematic_energy_add{0.0}
    , systematic_enable_weight_multiply{false}
    , systematic_weight_multiply{1.0}
    // cuts
    , threshold_low_energy_enable{false}
    , threshold_low_energy{0.05}
    , threshold_low_energy_sum_enable{false}
    , threshold_low_energy_sum{0.15}
    //, h_robustness_test_chi2{nullptr};
{

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


    if(0)
    {
        PrintInputData();
    }

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
    t = (TTree*)f->Get("NewElectronNtuplizer/NewElectronNtuplizer");

    t->SetBranchAddress("nElectrons", &nElectrons);
    t->SetBranchAddress("trueT1", &trueT1);
    t->SetBranchAddress("trueT2", &trueT2);
    t->SetBranchAddress("el_energy_", el_energy_);
    // new method requires weight to be saved to tree
    //gen_weight = 1.0;
    /*
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
    */

    // switch to enable data mode
    // defaut is MC
    #if 0
        analysis.SetModeFlag(MODE_FLAG::MODE_DATA);
    #endif


    ////////////////////////////////////////////////////////////////////////
    // INIT NON EPS31 DEPENDENT DATA
    ////////////////////////////////////////////////////////////////////////


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
    
    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    h_el_energy_sum_original = new TH1D((std::string("h_el_energy_sum_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    //h_el_energy_sum_original->SetStats(0);
    //h_el_energy_sum_original->SetLineColor(2);
    //h_el_energy_sum_original->SetMarkerColor(2);

    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    // same as above but re-weighted
    /*
    h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    */
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);
    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    // same as above but re-weighted
    /*
    h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore
    */


    // construct the baseline histograms
    EventLoop();

}


StaticsGroup::~StaticsGroup()
{
    std::cout << "~StaticsGroup" << std::endl;
    std::cout << "delete" << std::endl;
    delete h_el_energy_original;
    h_el_energy_original = nullptr;
    std::cout << "delete sum" << std::endl;
    delete h_el_energy_sum_original;
    h_el_energy_sum_original = nullptr;
    f->Close();
    delete f;
    f = nullptr;

    // clean optical correction data if loaded
    if(f_optical_correction != nullptr)
    {
        f_optical_correction->Close();
        delete f_optical_correction;
        //delete g_optical_correction;
        //g_optical_correction = nullptr;
        delete g_optical_correction_systematic_high;
        g_optical_correction_systematic_high = nullptr;
        delete g_optical_correction_systematic_low;
        g_optical_correction_systematic_low = nullptr;
    }

    std::cout << "done StaticsGroup" << std::endl;
}


void StaticsGroup::PrintInputData()
{

    // Print graphs for h_nEqNull, h_nEqTwo
    TCanvas *c_nEqNull;
    c_nEqNull = new TCanvas("c_nEqNull", "", 4000, 3000);
    h_nEqNull->Draw("colz");
    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");
    delete c_nEqNull;

    TCanvas *c_nEqTwo;
    c_nEqTwo = new TCanvas("c_nEqTwo", "", 4000, 3000);
    h_nEqTwo->Draw("colz");
    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    delete c_nEqTwo;

}


// create the baseline histogram
void StaticsGroup::EventLoop()
{

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
        // note name change
        Double_t T0{trueT1 / bb_Q};
        Double_t T1{trueT2 / bb_Q};

        // if MC apply energy degradation (correction)
        // only apply to MC
        if(mc_flag)
        {
            if(optical_correction_enable)
            {
                // apply energy degradation (optical correction)
                Double_t visible_true_ratio_0{1.0};
                Double_t visible_true_ratio_1{1.0};

                // if energy correction systematic is enabled, choose
                // energy correction value depending on which subanalysis
                // class this is
                // 0 is default
                // 1 is high systematic
                // -1 is low systematic
                // do not apply systematics here
                visible_true_ratio_0 = g_optical_correction->Eval(1000.0 * T0);
                visible_true_ratio_1 = g_optical_correction->Eval(1000.0 * T1);

                //std::cout << "visible_true_ratio = " << visible_true_ratio_0 << ", " << visible_true_ratio_1 << std::endl;

                // TODO this goes elsewhere
                // apply energy correction with systematics if enabled
                el_energy_0 = el_energy_0 * visible_true_ratio_0;
                el_energy_1 = el_energy_1 * visible_true_ratio_1;
            }
            else
            {
                // optical correction is diabled
                // NOOP
            }

            // TODO: other types of optical correction systematic
        
        }



        


        // generator weight (MC weight) multiplied by weight efficiency
        Double_t aux_weight{gen_weight};

        // TODO: energy degratation systematic

        // cut both electrons > 300 keV
        /*
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) continue;
            if(el_energy_1 < 0.3) continue;
        }
        */

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

#if 0
    if(0)
    {

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

    }
#endif

}



TH1D *const StaticsGroup::Get_h_el_energy_original() const
{
    return h_el_energy_original;
}

TH1D *const StaticsGroup::Get_h_el_energy_sum_original() const
{
    return h_el_energy_sum_original;
}

double StaticsGroup::Get_epsilon_31_baseline() const
{
    return epsilon_31_baseline;
}

Int_t StaticsGroup::Get_num_bins() const
{
    return num_bins;
}

TTree* StaticsGroup::Get_tree() const
{
    return t;
}

const Int_t& StaticsGroup::Get_nElectrons() const
{
    return nElectrons;
}

const Double_t& StaticsGroup::Get_trueT1() const
{
    return trueT1;
}

const Double_t& StaticsGroup::Get_trueT2() const
{
    return trueT2;
}

const Double_t* const StaticsGroup::Get_el_energy_() const
{
    return &(el_energy_[0]);
}

const Double_t& StaticsGroup::Get_gen_weight() const
{
    return gen_weight;
}

const Double_t& StaticsGroup::Get_bb_Q() const
{
    return bb_Q;
}

const double& StaticsGroup::Get_psiN0() const
{
    return psiN0;
}

const double& StaticsGroup::Get_psiN2() const
{
    return psiN2;
}

const TH2D* const StaticsGroup::Get_h_nEqNull() const
{
    return h_nEqNull;
}

const TH2D* const StaticsGroup::Get_h_nEqTwo() const
{
    return h_nEqTwo;
}

////////////////////////////////////////////////////////////////////////////////
// GET / SET PROGRAM BEHAVIOUR FLAGS / VARIABLES
////////////////////////////////////////////////////////////////////////////////

void StaticsGroup::SetStatisticsRobustnessTestEnable(const bool enable)
{
    if(statistics_robustness_test_enable == false)
    {
        if(enable == true)
        {
            statistics_robustness_test_chi2.clear();
        }
    }

    statistics_robustness_test_enable = enable;

    if(statistics_robustness_test_enable == true)
    {
        // TODO: what limits?
        //h_robustness_test_chi2 = new TH1D("h_robustness_test_chi2", "", 20, 0.0, 0.0);
    }
    else
    {
        //delete h_robustness_test_chi2;
        //h_robustness_test_chi2 = nullptr;
    }
}

bool StaticsGroup::GetStatisticsRobustnessTestEnable() const
{
    return statistics_robustness_test_enable;
}

void StaticsGroup::SetStatisticsRobustnessTestIterations(const int iterations)
{
    statistics_robustness_test_iterations = iterations;
}

int StaticsGroup::GetStatisticsRobustnessTestIterations() const
{
    return statistics_robustness_test_iterations;
}

void StaticsGroup::SetStatisticsRobustnessTestCurrentIteration(const int current)
{
    statistics_robustness_test_current_iteration = current;
}

int StaticsGroup::GetStatisticsRobustnessTestCurrentIteration() const
{
    return statistics_robustness_test_current_iteration;
}

void StaticsGroup::StoreChi2Result(const double chi2)
{
    statistics_robustness_test_chi2.push_back(chi2);
}

void StaticsGroup::StoreEpsilon31Result(const double epsilon31)
{
    statistics_robustness_test_epsilon_31.push_back(epsilon31);
}

void StaticsGroup::StoreEpsilon31ErrResult(const double epsilon31err)
{
    statistics_robustness_test_epsilon_31_err.push_back(epsilon31err);
}


void StaticsGroup::PrintChi2Results() const
{
    auto result{std::minmax_element(statistics_robustness_test_chi2.begin(), statistics_robustness_test_chi2.end())};
    double x_low{*result.first};
    double x_high{*result.second};
    Int_t _num_bins{20};
    double bin_width{(x_high - x_low) / (double)_num_bins};
    TH1D *h_statistics_robustness_test_chi2 = new TH1D("h_statistics_robustness_test_chi2", "", _num_bins, x_low, x_high + bin_width);
    std::vector<double>::const_iterator it{statistics_robustness_test_chi2.cbegin()};
    for(; it != statistics_robustness_test_chi2.cend(); ++ it)
    {
        h_statistics_robustness_test_chi2->Fill(*it);
    }

    TCanvas *c_statistics_robustness_test_chi2 = new TCanvas("c_statistics_robustness_test_chi2", "c_statistics_robustness_test_chi2", 600, 480);
    h_statistics_robustness_test_chi2->SetStats(0);
    h_statistics_robustness_test_chi2->GetXaxis()->SetTitle("#chi^{2}");
    h_statistics_robustness_test_chi2->GetYaxis()->SetTitle("Count (n=100)");
    h_statistics_robustness_test_chi2->Draw("hist");
    c_statistics_robustness_test_chi2->SaveAs("c_statistics_robustness_test_chi2.png");
    delete c_statistics_robustness_test_chi2;
    c_statistics_robustness_test_chi2 = nullptr;

    TFile *f_out = new TFile("c_statistics_robustness_test_chi2.root", "recreate");
    h_statistics_robustness_test_chi2->Write();
    f_out->Close();
    delete f_out;
    f_out = nullptr;
    
}


void StaticsGroup::PrintEpsilon31Results() const
{

    Int_t _num_bins{20};


    ////////////////////////////////////////////////////////////////////////////
    // parameter
    ////////////////////////////////////////////////////////////////////////////

    if(statistics_robustness_test_epsilon_31.size() > 0)
    {
        auto result{std::minmax_element(statistics_robustness_test_epsilon_31.begin(), statistics_robustness_test_epsilon_31.end())};
        double x_low{*result.first};
        double x_high{*result.second};
        double bin_width{(x_high - x_low) / (double)_num_bins};
        TH1D *h_statistics_robustness_test_epsilon_31 = new TH1D("h_statistics_robustness_test_epsilon_31", "", _num_bins, x_low, x_high + bin_width);
        std::vector<double>::const_iterator it{statistics_robustness_test_epsilon_31.cbegin()};
        for(; it != statistics_robustness_test_epsilon_31.cend(); ++ it)
        {
            //std::cout << *it << std::endl;
            h_statistics_robustness_test_epsilon_31->Fill(*it);
        }

        TCanvas *c_statistics_robustness_test_epsilon_31 = new TCanvas("c_statistics_robustness_test_epsilon_31", "c_statistics_robustness_test_epsilon_31", 600, 480);
        h_statistics_robustness_test_epsilon_31->SetStats(0);
        h_statistics_robustness_test_epsilon_31->GetXaxis()->SetTitle("#xi_{31}");
        h_statistics_robustness_test_epsilon_31->GetYaxis()->SetTitle("Count (n=100)");
        h_statistics_robustness_test_epsilon_31->Draw("hist");
        c_statistics_robustness_test_epsilon_31->SaveAs("c_statistics_robustness_test_epsilon_31.png");
        delete c_statistics_robustness_test_epsilon_31;
        c_statistics_robustness_test_epsilon_31 = nullptr;

        TFile *f_out = new TFile("c_statistics_robustness_test_epsilon_31.root", "recreate");
        h_statistics_robustness_test_epsilon_31->Write();
        f_out->Close();
        delete f_out;
        f_out = nullptr;

    }

    ////////////////////////////////////////////////////////////////////////////
    // error
    ////////////////////////////////////////////////////////////////////////////

    if(statistics_robustness_test_epsilon_31_err.size() > 0)
    {
        auto result_err{std::minmax_element(statistics_robustness_test_epsilon_31_err.begin(), statistics_robustness_test_epsilon_31_err.end())};
        double x_low_err{*result_err.first};
        double x_high_err{*result_err.second};
        Int_t _num_bins{20};
        double bin_width_err{(x_high_err - x_low_err) / (double)_num_bins};
        TH1D *h_statistics_robustness_test_epsilon_31_err = new TH1D("h_statistics_robustness_test_epsilon_31_err", "", _num_bins, x_low_err, x_high_err + bin_width_err);
        std::vector<double>::const_iterator it_err{statistics_robustness_test_epsilon_31_err.cbegin()};
        for(; it_err != statistics_robustness_test_epsilon_31_err.cend(); ++ it_err)
        {
            //std::cout << *it_err << std::endl;
            h_statistics_robustness_test_epsilon_31_err->Fill(*it_err);
        }

        TCanvas *c_statistics_robustness_test_epsilon_31_err = new TCanvas("c_statistics_robustness_test_epsilon_31_err", "c_statistics_robustness_test_epsilon_31_err", 600, 480);
        h_statistics_robustness_test_epsilon_31_err->SetStats(0);
        h_statistics_robustness_test_epsilon_31_err->GetXaxis()->SetTitle("#xi_{31} 1#sigma");
        h_statistics_robustness_test_epsilon_31_err->GetYaxis()->SetTitle("Count (n=100)");
        h_statistics_robustness_test_epsilon_31_err->Draw("hist");
        c_statistics_robustness_test_epsilon_31_err->SaveAs("c_statistics_robustness_test_epsilon_31_err.png");
        delete c_statistics_robustness_test_epsilon_31_err;
        c_statistics_robustness_test_epsilon_31_err = nullptr;

        TFile *f_out = new TFile("c_statistics_robustness_test_epsilon_31_err.root", "recreate");
        h_statistics_robustness_test_epsilon_31_err->Write();
        f_out->Close();
        delete f_out;
        f_out = nullptr;

    }

}


////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

void StaticsGroup::SetOpticalCorrectionEnable(const bool enable)
{
    if(enable)
    {
        optical_correction_enable = true;

        // load optical correction data
        if(f_optical_correction == nullptr)
        {
            f_optical_correction = new TFile("Correction_Birks_Cerenkov_avec_systematiques.root");
            if(f_optical_correction->IsOpen())
            {
                g_optical_correction = (TGraphErrors*)f_optical_correction->Get("Graph");

                // get number of points and allocate memory
                Int_t count{g_optical_correction->GetN()};
                Double_t *p_x{new Double_t[count]};
                Double_t *p_y_1{new Double_t[count]};
                Double_t *p_y_2{new Double_t[count]};

                // construct the systematic shifted graphs
                for(Int_t ix{0}; ix < count; ++ ix)
                {
                    Double_t x, y;
                    g_optical_correction->GetPoint(ix, x, y);
                    Double_t ye{g_optical_correction->GetErrorY(ix)};
                    p_x[ix] = x;
                    p_y_1[ix] = y + ye;
                    p_y_2[ix] = y - ye;
                }

                g_optical_correction_systematic_high = new TGraph(count, p_x, p_y_1);
                g_optical_correction_systematic_low = new TGraph(count, p_x, p_y_2);

                delete p_x;
                delete p_y_1;
                delete p_y_2;

                p_x = nullptr;
                p_y_1 = nullptr;
                p_y_2 = nullptr;

            }
            else
            {
                std::cerr << "Error: Could not open Correction_Birks_Cerenkov_avec_systematiques.root" << std::endl;
                throw "Correction_Birks_Cerenkov_avec_systematiques.root";
            }
        }

    }
    else
    {
        optical_correction_enable = false;
    }
}
