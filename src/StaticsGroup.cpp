#include "StaticsGroup.hpp"


#include "TCanvas.h"


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


void StaticsGroup::CreatePrintIntermediateData()
{

    delete h_data_0;
    delete h_data_1;
    delete h_data_2;


    ////////////////////////////////////////////////////////////////////////////
    // CREATE INTERMEDIATE DATA
    // Format: Nx3 array, each array a different value of epsilon
    ////////////////////////////////////////////////////////////////////////////

    // create data array for "complete data"
    // (with phase space factors)
    // epsilon_31 = 0.0
    create_data_with_phase_space_factor(data_0, 0.0, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    // epsilon_31 = 0.4
    create_data_with_phase_space_factor(data_1, 0.4, data_nEqNull, data_nEqTwo, psiN0, psiN2);
    
    // epsilon_31 = 0.8
    create_data_with_phase_space_factor(data_2, 0.8, data_nEqNull, data_nEqTwo, psiN0, psiN2);

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
    delete c_data_0;

    TCanvas *c_data_1;
    c_data_1 = new TCanvas("c_data_1", "", 4000, 3000);
    h_data_1->Draw("colz");
    c_data_1->SaveAs("c_data_1.png");
    c_data_1->SaveAs("c_data_1.pdf");
    c_data_1->SaveAs("c_data_1.C");
    delete c_data_1;

    TCanvas *c_data_2;
    c_data_2 = new TCanvas("c_data_2", "", 4000, 3000);
    h_data_2->Draw("colz");
    c_data_2->SaveAs("c_data_2.png");
    c_data_2->SaveAs("c_data_2.pdf");
    c_data_2->SaveAs("c_data_2.C");
    delete c_data_2;

}
