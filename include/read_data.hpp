#ifndef READDATA_HPP
#define READDATA_HPP


#include "TH2.h"


#include <fstream>
#include <iostream>
#include <vector>


void read_phase_space_factor(const char* buffer, double &value)
{

    const char* p{buffer};
    while(*p != '0')
    {
        ++ p;
    }
    std::string value_string;
    while((*p != '\0') && (*p != ' '))
    {
        value_string.push_back(*p);
        ++ p;
    }

    value = std::stod(value_string);

}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<T>> v)
{
    typename std::vector<std::vector<T>>::const_iterator it{v.cbegin()};
    for(; it != v.cend(); ++ it)
    {
        os << *it << "\n";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v)
{
    typename std::vector<T>::const_iterator it{v.cbegin()};
    for(; it != v.cend(); ++ it)
    {
        os << *it << " ";
    }
    return os;
}

// non numeric characters (not 0-9 and not special)
bool is_alpha(const char ch)
{
    if(' ' <= ch && ch <= '/')
    {
        return true;
    }
    else if(':' <= ch && ch <= '~')
    {
        return true;
    }
    else return false;
}

bool is_numeric_start(const char ch)
{
    if('.' == ch)
    {
        return true;
    }
    else if('0' <= ch && ch <= '9')
    {
        return true;
    }
    else return false;
}

bool is_numeric(const char ch)
{
    
}

// TODO: look for simpler approach by searching for ',' and splitting

void read_data_2(const char* buffer, std::vector<std::vector<double>>& data, const char delimiter = ' ')
{

    data.clear();
    data.emplace_back(std::vector<double>());

    std::string word;
    const char *p{buffer};
    for(;;)
    {
        if(*p == '\n' || *p == '\0' || *p == delimiter)
        {
            double value{atof(word.c_str())};
            word.clear();
            data[data.size() - 1].emplace_back(value);

            if(*p == '\n')
            {
                data.emplace_back(std::vector<double>());
            }
            else if(*p == '\0')
            {
                break;
            }
        }
        else
        {
            word.push_back(*p);
        }
        ++ p;
    }

    /*
    const char* p{buffer};
    const char* p_end;
    while(*p != '\0')
    {
        // skim all non alpha chars
        while(is_alpha(*p)) ++ p;
        if(is_numeric_start(*p))
        {
            p_end = p;
            int number_of_periods{0};
            int number_of_exponents{0}; // TODO: this, plus, minus at start
            if(*p_end == '.') ++ number_of_periods;
            ++ p_end;
            while(is_numeric(*p_end))
            {
                if(*p_end == '.') ++ number_of_periods;
                ++ p_end;
            }
            if(number_of_periods > 1)
        }
    }
    */

}

void read_data(const char* buffer, std::vector<std::vector<double>>& data, const char delimiter = ' ')
{
    
    data.clear();
    data.emplace_back(std::vector<double>());

    //std::cout << buffer << std::endl;

    const char* p{buffer};
    char* p_end;
    while(*p != '\0')
    {
        //std::cout << *p;
        //while(*p == ' ') ++ p;
        while(*p == delimiter) ++ p;
        double value{strtod(p, &p_end)};
        if(p == p_end)
        {
            data.pop_back();
            break;
        }
        data.back().emplace_back(value);
        p = p_end;
        if(*p == '\n')
        {
            //std::cout << data.at(data.size() - 1) << std::endl;
            data.emplace_back(std::vector<double>());
        }
    }

    //std::cout << "data read" << std::endl;
    //std::cout << "dimension x: " << data.size() << std::endl;
    //for(std::size_t i{0}; i < data.size(); ++ i)
    //{
    //    std::cout << "x=" << i << " dimension y: " << data[i].size() << std::endl;
    //    std::cin.get();
    //}
    //std::cout << "dimension y: " << data.at(0).size() << std::endl;

}


// write the data
void write_data_helper_2(const std::string& filename_data,
                         const std::vector<std::vector<double>> &data,
                         const char delimiter = ' ')
{

    ////////////////////////////////////////////////////////////////////////////
    // WRITE DATA OUT TO FILE
    ////////////////////////////////////////////////////////////////////////////

    std::ofstream ofs_data(filename_data.c_str());
    //std::streampos ifs_size{ifs_data.tellg()};
    //char * buf{new char[ifs_size + 1]};
    //ifs_data.seekg(0);
    //ifs_data.read(buf, ifs_size);
    //ifs_data.close();
    //std::vector<std::vector<double>> data_temp;
    //write_data_2(buf, data_temp, delimiter);
    //write_data_2(buf, data, delimiter);
    //data = data_temp;

    std::cout << "data.size()=" << data.size() << std::endl;
    std::size_t j{1};
    if(j < data.size())
    {
        for(;;)
        {
            //std::cout << "data.at(j=" << j << ").size()=" << data.size() << std::endl;
            std::size_t i{0};
            if(i < data.at(j).size())
            {
                for(;;)
                {
                    //std::cout << "j=" << j << " i=" << i << std::endl;
                    ofs_data << data.at(j).at(i);
                    if(++ i < data.at(j).size()) ofs_data << ',';
                    else break;
                }
            }
            if(++ j < data.size()) ofs_data << '\n';
            else break;
        }
    }

    ofs_data.flush();
    ofs_data.close();


    return;

}

// read the 1D data (data in electron single and sum energy histograms from paper)
void read_data_helper_2(const std::string& filename_data,
                        std::vector<std::vector<double>> &data,
                        const char delimiter = ' ')
{

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILE
    ////////////////////////////////////////////////////////////////////////////

    std::ifstream ifs_data(filename_data.c_str(), std::ios::ate);
    std::streampos ifs_size{ifs_data.tellg()};
    std::streamoff ifs_extra{1};
    char * buf{new char[ifs_size + ifs_extra]};
    ifs_data.seekg(0);
    ifs_data.read(buf, ifs_size);
    ifs_data.close();
    std::vector<std::vector<double>> data_temp;
    //read_data_2(buf, data_temp, delimiter);
    read_data(buf, data_temp, delimiter);
    delete buf;
    buf = nullptr;
    data = data_temp;
    
    return;

}

// read the 2D data (theory data)
// including phase space factors x2 and decay rate data x2
void read_data_helper(const std::string& filename_nEqNull, const std::string& filename_nEqTwo,
                      const std::string& filename_psiN0, const std::string& filename_psiN2,
                      std::vector<std::vector<double>> &data_nEqNull_ret, std::vector<std::vector<double>> &data_nEqTwo_ret,
                      double& psiN0_ret, double& psiN2_ret)
{

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

    // ifstream objects
    //std::ifstream ifs_nEqNull("nEqNull.dat", std::ios::ate);
    std::ifstream ifs_nEqNull(filename_nEqNull.c_str(), std::ios::ate);
    //std::ifstream ifs_nEqNull("test.dat", std::ios::ate);
    //std::ifstream ifs_nEqTwo("nEqTwo.dat", std::ios::ate);
    std::ifstream ifs_nEqTwo(filename_nEqTwo.c_str(), std::ios::ate);

    //std::ifstream ifs_psiN0("psiN0.txt", std::ios::ate);
    std::ifstream ifs_psiN0(filename_psiN0.c_str(), std::ios::ate);
    //std::ifstream ifs_psiN2("psiN2.txt", std::ios::ate);
    std::ifstream ifs_psiN2(filename_psiN2.c_str(), std::ios::ate);

    // get size
    std::streampos ifs_nEqNull_size{ifs_nEqNull.tellg()};
    std::streampos ifs_nEqTwo_size{ifs_nEqTwo.tellg()};
    
    std::streampos ifs_psiN0_size{ifs_psiN0.tellg()};
    std::streampos ifs_psiN2_size{ifs_psiN2.tellg()};

    std::streamoff ifs_extra{1};

    // allocate buffers
    char* buf_nEqNull{new char[ifs_nEqNull_size + ifs_extra]};
    char* buf_nEqTwo{new char[ifs_nEqTwo_size + ifs_extra]};
    
    char* buf_psiN0{new char[ifs_psiN0_size + ifs_extra]};
    char* buf_psiN2{new char[ifs_psiN2_size + ifs_extra]};

    // seek
    ifs_nEqNull.seekg(0);
    ifs_nEqTwo.seekg(0);

    ifs_psiN0.seekg(0);
    ifs_psiN2.seekg(0);

    // read
    ifs_nEqNull.read(buf_nEqNull, ifs_nEqNull_size);
    ifs_nEqTwo.read(buf_nEqTwo, ifs_nEqTwo_size);

    ifs_psiN0.read(buf_psiN0, ifs_psiN0_size);
    ifs_psiN2.read(buf_psiN2, ifs_psiN2_size);

    // close streams
    ifs_nEqNull.close();
    ifs_nEqTwo.close();

    ifs_psiN0.close();
    ifs_psiN2.close();

    // allocate data storage
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;

    double psiN0;
    double psiN2;

    // extract phase space factors
    //std::cout << buf_psiN0 << std::endl;
    //std::cout << buf_psiN2 << std::endl;
    read_phase_space_factor(buf_psiN0, psiN0);
    read_phase_space_factor(buf_psiN2, psiN2);

    std::cout << psiN0 << std::endl;
    std::cout << psiN2 << std::endl;

    // read data
    read_data(buf_nEqNull, data_nEqNull);
    std::cout << "Finished reading " << "nEqNull.dat" << std::endl;
    read_data(buf_nEqTwo, data_nEqTwo);
    std::cout << "Finished reading " << "nEqTwo.dat" << std::endl;
    
    delete buf_nEqNull;
    delete buf_nEqTwo;

    buf_nEqNull = nullptr;
    buf_nEqTwo = nullptr;

    delete buf_psiN0;
    delete buf_psiN2;

    buf_psiN0 = nullptr;
    buf_psiN2 = nullptr;

    //std::cout << buf_nEqNull << std::endl;
    //std::cout << buf_nEqTwo << std::endl;
    
    // return data
    data_nEqNull_ret = data_nEqNull;
    data_nEqTwo_ret = data_nEqTwo;
    psiN0_ret = psiN0;
    psiN2_ret = psiN2;
    

}

void convert_data_to_histogram_format(const std::vector<std::vector<double>> &data_nEqNull,
                                      const std::vector<std::vector<double>> &data_nEqTwo,
                                      TH2D * &h_nEqNull_return,
                                      TH2D * &h_nEqTwo_return)
{

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    const Int_t dimension_xy{1001};
    // don't plot raw data
    TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_nEqNull->SetStats(0);
    h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            //if(i + j < dimension_xy - 1)
            //{
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            //}
        }
    }
    std::cout << "Finished constructing input data histograms" << std::endl;
    
    h_nEqNull_return = h_nEqNull;
    h_nEqTwo_return = h_nEqTwo;

}

void create_data_with_phase_space_factor(std::vector<std::vector<double>> &data_ret, const double epsilon,
                                         std::vector<std::vector<double>> &data_nEqNull,
                                         std::vector<std::vector<double>> &data_nEqTwo,
                                         const double psiN0, const double psiN2)
{
    
    std::vector<std::vector<double>> data;
    data.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < data.size(); ++ i)
    {
        data[i].resize(data_nEqNull[i].size());
    }

    // epsilon_31 = 
    const double epsilon_31{epsilon};
    for(std::size_t i{0}; i < data.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2; ++ j)
        {
            data[i][j] = data_nEqNull[i][j]; // NOTE: this is not used
        }
        double phase_space_factor{1.0 / (psiN0 + epsilon_31 * psiN2)};
        data[i][2] = phase_space_factor * (data_nEqNull[i][2] + epsilon_31 * data_nEqTwo[i][2]);
    }
    
    data_ret = data;
    
}

#endif
