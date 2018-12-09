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
        {

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

            std::cout << "operator() : epsilon_31=" << epsilon_31 << std::endl;
            //std::cin.get();



            // reweighted electron energy after reweighting by epsilon (model)
            TH1D *h_el_energy_reweight{nullptr};
            // reweigted summed electron energy (model)
            TH1D *h_el_energy_sum_reweight{nullptr};

            // allocate
            Int_t num_bins{staticsgroup.Get_num_bins()};
            std::string h_name_append;

            h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
            h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
            
            const TH1D *h_el_energy_sum_original{staticsgroup.Get_h_el_energy_sum_original()};
            /*const*/ TH1D *h_el_energy_original{staticsgroup.Get_h_el_energy_original()};
            
            // fill histograms
            EventLoop(epsilon_31, h_el_energy_reweight, h_el_energy_sum_reweight, h_el_energy_sum_original);

            // fit histograms (summed energy)
            double amplitude{1.0};
            Fit(h_el_energy_sum_reweight, h_el_energy_sum_original, amplitude);

            // scale
            // scale the red histogram using the amplitude parameter
            h_el_energy_sum_reweight->Scale(amplitude);
            h_el_energy_reweight->Scale(amplitude);
            //h_el_energy_2d_reweight->Scale(amplitude);

            // calculate chi2
            double chi2{0.0};

            // Note: need to change to use systematics chisquare measurement
            // implemented in 100-Mo (original code)
            chi2 = chi_square_test(h_el_energy_reweight, h_el_energy_original);


            // print histograms each loop
            // this minimization is going wrong - I think - printing warning messages
            // this bug fixed
            std::string canvas_name{"c_iteration_"};
            canvas_name += std::to_string(static_iteration_counter);
            //canvas_name += std::string(".png"); // TODO: other outputs
            TCanvas *canvas_iteration_output = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 640, 480);
            std::cout << canvas_name << std::endl;
            h_el_energy_original->SetStats(0); // TODO: put elsewhere
            h_el_energy_original->SetMarkerColor(4); // TODO: put elsewhere
            h_el_energy_original->SetLineColor(4); // TODO: put elsewhere
            h_el_energy_reweight->SetStats(0);
            h_el_energy_reweight->SetMarkerColor(2);
            h_el_energy_reweight->SetLineColor(2);
            h_el_energy_original->SetMaximum(2.5e5);
            TLegend *l = new TLegend(0.55, 0.60, 0.85, 0.85, "Single Electron Energy");
            l->SetBorderSize(0);
            l->AddEntry(h_el_energy_original, "Baseline", "L");
            l->AddEntry(h_el_energy_reweight, "ReWeight", "L");
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
            TLatex *latex = new TLatex(2.5, 0.8e5, latex_string.c_str());
            h_el_energy_original->Draw("E");
            h_el_energy_reweight->Draw("Esame");
            l->Draw();
            latex->Draw();
            const std::string canvas_dir("./iteration_canvas_output/");
            canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".png")).c_str());
            canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".pdf")).c_str());
            canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".eps")).c_str());
            canvas_iteration_output->SaveAs((canvas_dir + canvas_name + std::string(".C")).c_str());
            delete canvas_iteration_output;
            canvas_iteration_output = nullptr;

            // read only object, must manage all memory here
            delete h_el_energy_reweight;
            h_el_energy_reweight = nullptr;
            delete h_el_energy_sum_reweight;
            h_el_energy_sum_reweight = nullptr;


            // increment the call counter
            ++ static_iteration_counter;

            return chi2;

        }


        void EventLoop(const double epsilon_31, TH1D *h_el_energy_reweight, TH1D *h_el_energy_sum_reweight, const TH1D *h_el_energy_sum_original) const
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


        void Fit(/*const double epsilon_31,*/ TH1D *h_el_energy_sum_reweight, const TH1D *h_el_energy_sum_original, double &amplitude) const
        {

            ////////////////////////////////
            // Analysis.SummedEnergyFit() //
            ////////////////////////////////


            std::cout << "fit" << std::endl;

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



            std::cout << "fit finished" << std::endl;



        }




    protected:

        static unsigned long long static_iteration_counter;

        double error_def;

        const StaticsGroup& staticsgroup;


};


#endif // MINIMIZE_FCN_EPSILON_31_HPP
