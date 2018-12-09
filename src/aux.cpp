#include "aux.hpp"

#include <TF1.h>
//#include <TF2.h>


void print_histo_text(std::ostream& dst, const std::string& title, const TH1* const histo)
{
    std::cout << title << std::endl;
    for(Int_t i{0}; i <= histo->GetNbinsX(); ++ i)
    {
        std::cout << "i=" << i << "   " << histo->GetBinContent(i) << "   " << histo->GetBinError(i) << std::endl;
    }
}

/*
Double_t fit_function_2d(Double_t *_x, Double_t *par)
{

    Double_t amplitude{par[0]};

    Int_t number_of_bins{40};

    // get x,y
    Double_t x{_x[0]};
    Double_t y{_x[1]};

    // get histogram
    char* p{(char*)par};
    TH2D *histo{(TH2D*)(p + sizeof(Double_t))};

    Int_t index{histo->FindBin(x, y)};
    Double_t content{histo->GetBinContent(index)};

    return (content * amplitude);
}
*/
// TODO: this method for 1d fit


// should the parameter go in the "function" or the "data"?
// currently I have it in the data, not the function, which is weird?

// x[0] stores the current x value, par[0] stores the "parameter" (amplitude)
// par[1 ...] stores the data for the other histogram, however these parameters
// are constants
// NOTE: par also includes the "position of the end of the final bin"
Double_t fit_function(Double_t *x, Double_t *par)
{
    Double_t amplitude{par[0]};
    //std::cout << "amplitude=" << amplitude << std::endl;
    Int_t number_of_bins{40}; // TODO: change to global variable?

    // expand par into x and y values
    // NOTE: par also includes the "position of the end of the final bin"
    Double_t *x_values{new Double_t[number_of_bins + 1]};
    Double_t *y_values{new Double_t[number_of_bins + 1]};

    for(Int_t i{0}; i <= number_of_bins; ++ i)
    {
        Int_t ix{2 * i + 1};
        Int_t jx{2 * i + 2};

        x_values[i] = par[ix];
        y_values[i] = par[jx];

        //std::cout << "x=" << x_values[i] << " y=" << y_values[i] << std::endl;
    }

    Double_t xx{x[0]};
    
    // find bin in which x lies
    if(xx < x_values[0])
    {
        // x too small
        std::cerr << "x too small" << std::endl;
        std::cout << "x = " << xx << std::endl;
    }

    if(xx > x_values[number_of_bins])
    {
        // x too large
        std::cerr << "x too large" << std::endl;
        std::cout << "x = " << xx << std::endl;
    }

    //debug
    bool found{false};
    Int_t found_bin_index{-1};

    for(Int_t i{0}; i <= number_of_bins; ++ i)
    {
        if(x_values[i] <= xx)
        {
            if(xx < x_values[i + 1])
            {
                found = true;
                found_bin_index = i;
            }
        }
    }
    if(!found)
    {
        if(xx == x_values[number_of_bins])
        {
            found = true;
            found_bin_index = number_of_bins - 1;
        }
    }

    if(!found)
    {
        std::cerr << "error, bin not found in fit_function" << std::endl;
        std::cout << "x = " << xx << std::endl;
        std::cout << x_values[number_of_bins] << std::endl;
        throw "problem";
    }


    // get y value (content) of this bin
    // return content
    Double_t ret{amplitude * y_values[found_bin_index]};

    // clean
    delete x_values;
    delete y_values;
    x_values = nullptr;
    y_values = nullptr;

    // return
    return ret;
}


// chi square test histogram to TGraph
Double_t chi_square_test(const TGraph* const graph_fit, const TH1* const histo_data)
{
    //std::cout << "START" << std::endl;

    // number of bins
    Int_t nbx{histo_data->GetNbinsX()};
    
    // check same number of bins in each histo
    //if(nbx != histo_data->GetNbinsX())
    //    throw "bin number mismatch in function histo_histo_chisquare";

    Double_t chisquare{0.0};
    for(Int_t ix{1}; ix <= nbx; ++ ix)
    {
        //Double_t content_f{histo_fit->GetBinContent(ix)};
        Double_t x{histo_data->GetBinCenter(ix)};
        Double_t content_f{graph_fit->Eval(x)};
        Double_t content_d{histo_data->GetBinContent(ix)};
        Double_t delta{content_f - content_d};
        Double_t error_d{histo_data->GetBinError(ix)};
        if(error_d <= 0.0)
        {
            if(content_d == 0.0)
            {
                continue;
            }
            else
            {
                std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
            }
        }
        Double_t chi{std::pow(delta / error_d, 2.0)};
        //std::cout << "chi=" << chi << std::endl;
        chisquare += chi;
    }
    //chisquare /= (Double_t)histo_fit->GetNbinsX();
    return chisquare;
    
}


// histo_fit is the fit function
// histo_data is the data
// TODO: assumes ranges are the same
Double_t chi_square_test(const TH1* const histo_fit, const TH1* const histo_data, const Double_t min, const Double_t max)
{
    // number of bins
    Int_t nbx{histo_fit->GetNbinsX()};
    
    // check same number of bins in each histo
    if(nbx != histo_data->GetNbinsX())
        throw "bin number mismatch in function histo_histo_chisquare";

    Double_t chisquare{0.0};
    for(Int_t ix{1}; ix <= nbx; ++ ix)
    {
        Double_t bin_center{histo_fit->GetBinCenter(ix)};
        if(bin_center < min) continue;
        if(bin_center > max) continue;
        Double_t content_f{histo_fit->GetBinContent(ix)};
        Double_t content_d{histo_data->GetBinContent(ix)};
        Double_t delta{content_f - content_d};
        Double_t error_d{histo_data->GetBinError(ix)};
        if(error_d <= 0.0)
        {
            if(content_d == 0.0)
            {
                continue;
            }
            else
            {
                std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
            }
        }
        Double_t chi{std::pow(delta / error_d, 2.0)};
        chisquare += chi;
    }
    //chisquare /= (Double_t)histo_fit->GetNbinsX();
    return chisquare;
}


// histo_fit is the fit function
// histo_data is the data
Double_t chi_square_test(const TH1* const histo_fit, const TH1* const histo_data)
{
    // number of bins
    Int_t nbx{histo_fit->GetNbinsX()};
    
    // check same number of bins in each histo
    if(nbx != histo_data->GetNbinsX())
        throw "bin number mismatch in function histo_histo_chisquare";

    Double_t chisquare{0.0};
    for(Int_t ix{1}; ix <= nbx; ++ ix)
    {
        Double_t content_f{histo_fit->GetBinContent(ix)};
        Double_t content_d{histo_data->GetBinContent(ix)};
        Double_t delta{content_f - content_d};
        Double_t error_d{histo_data->GetBinError(ix)};
        if(error_d <= 0.0)
        {
            if(content_d == 0.0)
            {
                continue;
            }
            else
            {
                std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
            }
        }
        Double_t chi{std::pow(delta / error_d, 2.0)};
        chisquare += chi;

        //std::cout << "data=," << content_d << ",error=," << error_d << ",function=," << content_f << std::endl;
    }
    //std::cout << "chisquare=" << chisquare << std::endl;
    //chisquare /= (Double_t)histo_fit->GetNbinsX();
    return chisquare;
}


// driver for chisquare test
void driver(const TH1* const histo_fit, const TH1* const histo_data, Double_t &param)
{

    // delta param, how much to change param by each step
    Double_t delta_param{param};

    // iteration direction for delta
    //Double_t direction{1.0};

    // limit for delta_param
    Double_t delta_param_min{std::abs(0.001 * param)};
    
    // copy the data fit histogram
    TH1 *histo_fit_copy_plus = (TH1*)histo_fit->Clone();
    TH1 *histo_fit_copy = (TH1*)histo_fit->Clone();
    TH1 *histo_fit_copy_minus = (TH1*)histo_fit->Clone();
    
    // scale by amplitude
    histo_fit_copy_plus->Scale(param);
    histo_fit_copy->Scale(param);
    histo_fit_copy_minus->Scale(param);

    // get current chi square
    Double_t chi_plus{chi_square_test(histo_fit_copy_plus, histo_data)};
    Double_t chi_current{chi_square_test(histo_fit_copy, histo_data)};
    Double_t chi_minus{chi_square_test(histo_fit_copy_minus, histo_data)};

    // count iterations
    Long64_t iterations{0};

    for(;;)
    {
        std::cout << "param=" << param << " chi=" << chi_current << std::endl;
        
        // exit condition
        if(std::abs(delta_param) <= delta_param_min) break;

        // step param
        //param += direction * delta_param;
        Double_t param_plus{param + delta_param};
        Double_t param_minus{param - delta_param};

        // set histo_fit_copy
        histo_fit_copy_plus = (TH1*)histo_fit->Clone();
        histo_fit_copy = (TH1*)histo_fit->Clone();
        histo_fit_copy_minus = (TH1*)histo_fit->Clone();

        // scale by amplitude
        histo_fit_copy_plus->Scale(param_plus);
        histo_fit_copy->Scale(param);
        histo_fit_copy_minus->Scale(param_minus);

        // evaluate chisquare
        Double_t chi_plus{chi_square_test(histo_fit_copy_plus, histo_data)};
        Double_t chi_current{chi_square_test(histo_fit_copy, histo_data)};
        Double_t chi_minus{chi_square_test(histo_fit_copy_minus, histo_data)};

        // TODO: there are 6 possibilities for ordering A, B, C ?
        // some are valid and some are not
        //if(chi_plus < chi_current && chi_current < chi_minus)
        // stupid method

        // make decision based on chisquare
        if(chi_plus < chi_current)
        {
            // keep going in this direction (delta_param ok)
            if(chi_plus < chi_minus)
            {
                // chi plus is the minimum - go in this direction
                param = param_plus;
                chi_current = chi_plus;

                // this can probably never happen?
            }
            else if(chi_minus < chi_plus)
            {
                // chi minus is the minimum - go in this direction
                param = param_minus;
                chi_current = chi_minus;

                // this can probably never happen
            }
            else
            {
                // chi minus and chi plus are equal, but chi current is greater
                // this cannot happen
                throw "chi error";
            }
        }
        else if(chi_current < chi_plus)
        {
            if(chi_minus < chi_current)
            {
                // chi minus is the minimum - go in this direction
                param = param_minus;
                chi_current = chi_minus;
            }
            else if(chi_current < chi_minus)
            {
                // the minimum is somewhere around chi_current
                delta_param *= 0.8;
            }
            else
            {
                // chi_minus == chi_current
                // the minimum is somewhere around chi_current
                delta_param *= 0.8;
            }
        }
        else
        {
            // chi_plus == chi_current
            // the minimum is somewhere around chi_current
            delta_param *= 0.8;
        }

        //else
        //{
        //    if(direction > 0.0)
        //    {
        //        direction = -1.0;
        //    }
        //    else if(direction < 0.0)
        //    {
        //        direction = 1.0;
        //    }
        //
        //    // switch direction of delta param and make smaller
        //    delta_param *= 0.8;
        //}

        // set chi-square
        //chi_current = chi_new;

        ++ iterations;

    }

    std::cout << "iteration converged after " << iterations << " iterations" << std::endl;

}


// histo_fit is the fit function
// histo_data is the data
// 2d version
Double_t chi_square_test(const TH2* const histo_fit, const TH2* const histo_data)
{
    // number of bins
    Int_t nbx{histo_fit->GetNbinsX()};
    Int_t nby{histo_fit->GetNbinsY()};
    
    // check same number of bins in each histo
    if(nbx != histo_data->GetNbinsX())
        throw "bin number mismatch in function histo_histo_chisquare";
    if(nby != histo_data->GetNbinsY())
        throw "bin number mismatch in function histo_histo_chisquare";

    Double_t chisquare{0.0};
    for(Int_t iy{1}; iy <= nby; ++ iy)
    {
        for(Int_t ix{1}; ix <= nbx; ++ ix)
        {
            Double_t content_f{histo_fit->GetBinContent(ix, iy)};
            Double_t content_d{histo_data->GetBinContent(ix, iy)};
            Double_t delta{content_f - content_d};
            Double_t error_d{histo_data->GetBinError(ix, iy)};
            if(error_d <= 0.0)
            {
                if(content_d == 0.0)
                {
                    continue;
                }
                else
                {
                    std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
                }
            }
            Double_t chi{std::pow(delta / error_d, 2.0)};
            chisquare += chi;
        }
    }
    //chisquare /= (Double_t)histo_fit->GetNbinsX();
    return chisquare;
}


// histo_fit is the fit function
// histo_data is the data
// 2d version
Double_t chi_square_test(const TH2* const histo_fit, const TH2* const histo_data, const Double_t min, const Double_t max, Int_t &non_empty_bins_2d)
{

    non_empty_bins_2d = 0;

    // number of bins
    Int_t nbx{histo_fit->GetNbinsX()};
    Int_t nby{histo_fit->GetNbinsY()};
    
    // check same number of bins in each histo
    if(nbx != histo_data->GetNbinsX())
        throw "bin number mismatch in function histo_histo_chisquare";
    if(nby != histo_data->GetNbinsY())
        throw "bin number mismatch in function histo_histo_chisquare";

    Double_t chisquare{0.0};
    for(Int_t iy{1}; iy <= nby; ++ iy)
    {
        for(Int_t ix{1}; ix <= nbx; ++ ix)
        {
            Double_t bin_center_x{histo_fit->GetXaxis()->GetBinCenter(ix)};
            Double_t bin_center_y{histo_fit->GetYaxis()->GetBinCenter(iy)};
            /* note: only care about high energy region - y axis
            if(bin_center_x < min) continue;
            if(bin_center_x > max) continue;
            if(bin_center_y < min) continue;
            if(bin_center_y > max) continue;
            */
            if(bin_center_y < min) continue;
            if(bin_center_y > max) continue;
            Double_t content_f{histo_fit->GetBinContent(ix, iy)};
            Double_t content_d{histo_data->GetBinContent(ix, iy)};
            if(content_d != 0.0) ++ non_empty_bins_2d;
            Double_t delta{content_f - content_d};
            Double_t error_d{histo_data->GetBinError(ix, iy)};
            if(error_d <= 0.0)
            {
                if(content_d == 0.0)
                {
                    continue;
                }
                else
                {
                    std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_data, bin error is 0.0 but content != 0.0" << std::endl;
                }
            }
            Double_t chi{std::pow(delta / error_d, 2.0)};
            chisquare += chi;
        }
    }
    //chisquare /= (Double_t)histo_fit->GetNbinsX();
    return chisquare;
}



