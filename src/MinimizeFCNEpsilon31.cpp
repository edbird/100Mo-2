#include "MinimizeFCNEpsilon31.hpp"


#include "TRandom3.h"


unsigned long long MinimizeFCNEpsilon31::static_iteration_counter{0};



void MinimizeFCNEpsilon31::EventLoop(const double epsilon_31, TH1D *h_el_energy_reweight, TH1D *h_el_energy_sum_reweight, const TH1D *h_el_energy_sum_original) const
{

    // create data dependant on epsilon_31

    //delete h_el_energy_reweight;
    


    // these are re-made each iteration, TODO


    //h_el_energy_sum_reweight->SetStats(0);
    //h_el_energy_sum_reweight->SetLineColor(4);
    //h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore



    // TODO: I don't remember what to do about systematics, so remove for now


    //////////////////////////
    // Analysis.EventLoop() //
    //////////////////////////

    const double &epsilon_31_baseline{staticsgroup.Get_epsilon_31_baseline()};
    const Double_t &bb_Q{staticsgroup.Get_bb_Q()};
    std::cout << "bb_Q=" << bb_Q << std::endl;

    const TH2D* const h_nEqNull{staticsgroup.Get_h_nEqNull()};
    const TH2D* const h_nEqTwo{staticsgroup.Get_h_nEqTwo()};
    const double& psiN0{staticsgroup.Get_psiN0()};
    const double& psiN2{staticsgroup.Get_psiN2()};

    TTree *t{staticsgroup.Get_tree()};

    const Int_t &nElectrons{staticsgroup.Get_nElectrons()};
    const Double_t &trueT1{staticsgroup.Get_trueT1()};
    const Double_t &trueT2{staticsgroup.Get_trueT2()};
    const Double_t* const el_energy_{staticsgroup.Get_el_energy_()};
    const Double_t &gen_weight{staticsgroup.Get_gen_weight()};


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
#if 0
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
#endif

        
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
        Double_t weight{ReWeight3(T1, T2, epsilon_31_baseline, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
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


    }


    //print_histo_text(std::cout, "h_el_energy_reweight", h_el_energy_reweight);

    /*
    TFile *fout = new TFile("fout.root", "recreate");
    h_el_energy_reweight->Write();
    h_el_energy_sum_reweight->Write();
    fout->Close();
    std::cout << "fout" << std::endl;
    std::cin.get();

    TCanvas *c = new TCanvas("c_srw", "c_srw", 800, 600);
    h_el_energy_sum_reweight->Draw("e");
    c->SaveAs("c_srw.png");
    */

}


// populate std::vector with random gaussian distributed parameters for
// use in applying to reweight data to create pseudo reweight data
// (randomly, statistically fluctuated data)
void MinimizeFCNEpsilon31::CreateRandomizedStatistics()
{
    
    // get Gaussian number, width = 1.0, mean = 0.0
    TRandom3 gen;
    gen.SetSeed(0);

    // repeat twice, need enough gaussian numbers for 2 histos
    for(int r{0}; r < 2; ++ r)
    {
        // repeat for number of bins per histo, need one random
        // gaussian number for each bin
        Int_t num_bins{staticsgroup.Get_num_bins()};
        for(Int_t ix{1}; ix <= num_bins; ++ ix)
        {
            gen_gaussian_params.push_back(gen.Gaus());
        }
    }
    
    // get min max
    auto result{std::minmax_element(gen_gaussian_params.begin(), gen_gaussian_params.end())};
    double x_low{*result.first};
    double x_high{*result.second};
    if(std::abs(x_low) > std::abs(x_high))
    {
        x_high = std::abs(x_low);
    }
    else
    {
        x_low = -std::abs(x_high);
    }
    Int_t _num_bins{20};
    double bin_width{(x_high - x_low) / (double)_num_bins};

    // name
    std::string c_gen_gaussian_params_fname("c_gen_gaussian_params_");
    c_gen_gaussian_params_fname += std::to_string(staticsgroup.GetStatisticsRobustnessTestCurrentIteration());
    c_gen_gaussian_params_fname += std::string(".png");

    TH1D *h_gen_gaussian_params = new TH1D(c_gen_gaussian_params_fname.c_str(), "", _num_bins, x_low, x_high + bin_width);
    std::vector<double>::const_iterator it{gen_gaussian_params.cbegin()};
    for(; it != gen_gaussian_params.cend(); ++ it)
    {
        h_gen_gaussian_params->Fill(*it);
    }

    TCanvas *c_gen_gaussian_params = new TCanvas(c_gen_gaussian_params_fname.c_str(), c_gen_gaussian_params_fname.c_str(), 600, 480);
    h_gen_gaussian_params->SetStats(0);
    h_gen_gaussian_params->Draw("hist");
    c_gen_gaussian_params->SaveAs(c_gen_gaussian_params_fname.c_str());
    delete c_gen_gaussian_params;
    c_gen_gaussian_params = nullptr;

    h_gen_gaussian_params->SetDirectory(0);
    delete h_gen_gaussian_params;
    h_gen_gaussian_params = nullptr;

}


void MinimizeFCNEpsilon31::ApplyRandomizedStatistics(TH1D * const h_el_energy_reweight, TH1D * const h_el_energy_sum_reweight,
        TH1D * const h_el_energy_reweight_pseudo, TH1D * const h_el_energy_sum_reweight_pseudo) const
{


    std::vector<TH1D*> histo_vec;
    std::vector<TH1D*> histo_vec_out;
    histo_vec.emplace_back(h_el_energy_reweight);
    histo_vec.emplace_back(h_el_energy_sum_reweight);
    histo_vec_out.emplace_back(h_el_energy_reweight_pseudo);
    histo_vec_out.emplace_back(h_el_energy_sum_reweight_pseudo);
    std::vector<TH1D*>::const_iterator it{histo_vec.cbegin()};
    std::vector<TH1D*>::const_iterator it_out{histo_vec_out.cbegin()};
    std::size_t index{0};
    for(; it != histo_vec.cend() && it_out != histo_vec_out.cend(); ++ it, ++ it_out)
    {
        // set histogram
        TH1D* histo{*it};
        TH1D* histo_out{*it_out};

        // iterate over histogram bins, randomize content
        for(Int_t ix{1}; ix <= histo->GetNbinsX(); ++ ix)
        {

            //const double gen_gaussian{gen->Gaus()};
            const double gen_gaussian{gen_gaussian_params.at(index)};
            ++ index;

            const double error{histo->GetBinError(ix)};
            const double value{histo->GetBinContent(ix)};

            const double new_value{value + gen_gaussian * error};

            histo_out->SetBinContent(ix, new_value);

        }
    }

}



void MinimizeFCNEpsilon31::Fit(/*const double epsilon_31,*/ TH1D *h_el_energy_sum_reweight, const TH1D *h_el_energy_sum_original, double &amplitude) const
{

    ////////////////////////////////
    // Analysis.SummedEnergyFit() //
    ////////////////////////////////


    //std::cout << "fit" << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // MINIMIZATION OF SUMMED ENERGY HISTOGRAMS (FIT)
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
    MinimizeFCNSummedEnergy theFCN_amplitude(h_el_energy_sum_reweight, h_el_energy_sum_original);
    
    // initial parameters, uncertainties
    std::vector<double> init_par;
    init_par.push_back(1.0); // amplitude, start at 1.0
    std::vector<double> init_err;
    init_err.push_back(0.1);

    // create minimizer
    ROOT::Minuit2::VariableMetricMinimizer theMinimizer;

    // minimize
    ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN_amplitude, init_par, init_err);
    //std::cout << "minimum: " << FCN_min << std::endl;

    // TODO: must be a better method of getting chisquare
    /*
    std::cout << FCN_min.Parameters().Vec() << std::endl;
    std::vector<double> end_par;
    const double* vec_data{FCN_min.Parameters().Vec().Data()};
    for(int c{0}; c < FCN_min.Parameters().Vec().size(); ++ c)
    {
        end_par.push_back(vec_data[c]);
    }
    std::cout << "chi2=" << theFCN_amplitude.operator()(end_par) << std::endl;
    //double amplitude{end_par.at(0)};
    amplitude = end_par.at(0);
    */

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



    //std::cout << "fit finished" << std::endl;



}


void MinimizeFCNEpsilon31::PrintIteration(TH1* const h_el_energy_original, TH1* const h_el_energy_reweight, const double epsilon_31, const double chi2, const double amplitude) const
{

    // print histograms each loop
    // this minimization is going wrong - I think - printing warning messages
    // this bug fixed

    // construct canvas output name
    std::string canvas_name{"c_iteration_"};
    // add iteration number appended (output filename)
    canvas_name += std::to_string(staticsgroup.GetStatisticsRobustnessTestCurrentIteration()) + std::string("_");
    canvas_name += std::to_string(static_iteration_counter);
    //canvas_name += std::string(".png"); // TODO: other outputs
    std::cout << canvas_name << std::endl;

    // create canvas
    TCanvas *canvas_iteration_output = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 640, 480);
    
    // set histogram draw properties
    h_el_energy_original->SetStats(0); // TODO: put elsewhere
    h_el_energy_original->SetMarkerColor(4); // TODO: put elsewhere
    h_el_energy_original->SetLineColor(4); // TODO: put elsewhere
    h_el_energy_reweight->SetStats(0);
    h_el_energy_reweight->SetMarkerColor(2);
    h_el_energy_reweight->SetLineColor(2);
    h_el_energy_original->SetMaximum(2.5e5);

    // legend
    TLegend *l = new TLegend(0.55, 0.60, 0.85, 0.85, "Single Electron Energy");
    l->SetBorderSize(0);
    l->AddEntry(h_el_energy_original, "Baseline", "L");
    l->AddEntry(h_el_energy_reweight, "ReWeight", "L");
    l->Draw();

    // latex object label - string
    std::string latex_string("#splitline{#xi_{31}=");
    //latex_string += std::to_string(epsilon_31);
    std::ostringstream oss;
    oss << std::fixed;
    oss << std::setprecision(3);
    oss << epsilon_31;
    std::ostringstream oss_chi2;
    oss_chi2 << std::fixed;
    oss_chi2 << std::setprecision(3);
    oss_chi2 << chi2;
    latex_string += oss.str();
    latex_string += std::string("}{#chi^{2}=");
    latex_string += oss_chi2.str();
    latex_string += std::string("}");

    // latex object label - latex label
    TLatex *latex = new TLatex(2.5, 0.8e5, latex_string.c_str());
    h_el_energy_original->Draw("E");
    h_el_energy_reweight->Draw("Esame");
    latex->Draw();

    // save canvas
    // directory
    const std::string canvas_dir("./iteration_canvas_output/");
    // output formats
    std::vector<std::string> file_ext;
    file_ext.emplace_back(".png");
    file_ext.emplace_back(".pdf");
    file_ext.emplace_back(".eps");
    file_ext.emplace_back(".C");
    std::vector<std::string>::iterator it{file_ext.begin()};
    for(; it != file_ext.end(); ++ it)
    {
        std::string fullname(canvas_dir + canvas_name + (*it));
        canvas_iteration_output->SaveAs(fullname.c_str());
    }
    //canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".png")).c_str());
    //canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".pdf")).c_str());
    //canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".eps")).c_str());
    //canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".C")).c_str());

    // save in root format
    std::string f_out_fullname(canvas_dir + canvas_name + std::string(".root"));
    TFile *f_out = new TFile(f_out_fullname.c_str(), "recreate");
    canvas_iteration_output->Write();
    h_el_energy_original->Write();
    h_el_energy_reweight->Write();

    // clean
    f_out->Close();
    delete f_out;
    f_out = nullptr;

    // clean
    delete canvas_iteration_output;
    canvas_iteration_output = nullptr;

    // write parameter output
    std::string of_param_fullname(canvas_dir + canvas_name + std::string("_params.txt"));
    std::ofstream of_param(of_param_fullname.c_str());
    of_param << "chi2," << chi2 << std::endl;
    of_param << "amplitude" << amplitude << std::endl;
    of_param.flush();
    of_param.close();

}
