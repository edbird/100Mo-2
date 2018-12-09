#ifndef READDATA_HPP
#define READDATA_HPP


#include "TH2.h"


#include <fstream>
#include <iostream>
#include <vector>


void read_phase_space_factor(const char* buffer, double &value);

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<T>> v);

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v);

// non numeric characters (not 0-9 and not special)
bool is_alpha(const char ch);

bool is_numeric_start(const char ch);

bool is_numeric(const char ch);

// TODO: look for simpler approach by searching for ',' and splitting

void read_data_2(const char* buffer, std::vector<std::vector<double>>& data, const char delimiter = ' ');

void read_data(const char* buffer, std::vector<std::vector<double>>& data, const char delimiter = ' ');


// write the data
void write_data_helper_2(const std::string& filename_data,
                         const std::vector<std::vector<double>> &data,
                         const char delimiter = ' ');

// read the 1D data (data in electron single and sum energy histograms from paper)
void read_data_helper_2(const std::string& filename_data,
                        std::vector<std::vector<double>> &data,
                        const char delimiter = ' ');

// read the 2D data (theory data)
// including phase space factors x2 and decay rate data x2
void read_data_helper(const std::string& filename_nEqNull, const std::string& filename_nEqTwo,
                      const std::string& filename_psiN0, const std::string& filename_psiN2,
                      std::vector<std::vector<double>> &data_nEqNull_ret, std::vector<std::vector<double>> &data_nEqTwo_ret,
                      double& psiN0_ret, double& psiN2_ret);

void convert_data_to_histogram_format(const std::vector<std::vector<double>> &data_nEqNull,
                                      const std::vector<std::vector<double>> &data_nEqTwo,
                                      TH2D * &h_nEqNull_return,
                                      TH2D * &h_nEqTwo_return);

void create_data_with_phase_space_factor(std::vector<std::vector<double>> &data_ret, const double epsilon,
                                         std::vector<std::vector<double>> &data_nEqNull,
                                         std::vector<std::vector<double>> &data_nEqTwo,
                                         const double psiN0, const double psiN2);

#endif
