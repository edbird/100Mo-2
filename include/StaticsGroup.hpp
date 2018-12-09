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

    StaticsGroup(const double epsilon_31_baseline, const std::string& filename)
        : epsilon_31_baseline{epsilon_31_baseline} // TODO: pass to ReWeight2
        , filename{filename}
        , f{nullptr}
        , t{nullptr}
        , num_bins{40}
        , bb_Q{3.034}
        , dimension_xy{1001}
        , h_el_energy_original{nullptr}
        , h_el_energy_sum_original{nullptr}
        , h_data_0{nullptr}
        , h_data_1{nullptr}
        , h_data_2{nullptr}
        , gen_weight{1.0}
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

        if(0)
        {
            CreatePrintIntermediateData();
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


    virtual
    ~StaticsGroup()
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
        std::cout << "done StaticsGroup" << std::endl;
    }


    void PrintInputData();

    void CreatePrintIntermediateData();

    void EventLoop()
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

    TH1D *const Get_h_el_energy_original() const
    {
        return h_el_energy_original;
    }

    TH1D *const Get_h_el_energy_sum_original() const
    {
        return h_el_energy_sum_original;
    }

    double Get_epsilon_31_baseline() const
    {
        return epsilon_31_baseline;
    }

    Int_t Get_num_bins() const
    {
        return num_bins;
    }

    TTree* Get_tree() const
    {
        return t;
    }

    const Int_t& Get_nElectrons() const
    {
        return nElectrons;
    }

    const Double_t& Get_trueT1() const
    {
        return trueT1;
    }

    const Double_t& Get_trueT2() const
    {
        return trueT2;
    }

    const Double_t* const Get_el_energy_() const
    {
        return &(el_energy_[0]);
    }

    const Double_t& Get_gen_weight() const
    {
        return gen_weight;
    }

    const Double_t& Get_bb_Q() const
    {
        return bb_Q;
    }

    const double& Get_psiN0() const
    {
        return psiN0;
    }

    const double& Get_psiN2() const
    {
        return psiN2;
    }

    const TH2D* const Get_h_nEqNull() const
    {
        return h_nEqNull;
    }

    const TH2D* const Get_h_nEqTwo() const
    {
        return h_nEqTwo;
    }



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
    std::vector<std::vector<double>> data_0;
    std::vector<std::vector<double>> data_1;
    std::vector<std::vector<double>> data_2;
    TH2D *h_data_0;
    TH2D *h_data_1;
    TH2D *h_data_2;


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

};



#endif // STATICS_GROUP
