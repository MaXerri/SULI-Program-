#pragma once
#ifndef CONVERSIONS_H
#define CONVERSIONS_H
#include <vector>
#include <string>

std::vector<double> toReducePrimG6(std::string dc7unsrtdist_mat_out);

std::vector<double> todc7unsrt(std::vector<double> redprimg6);

std::tuple<std::vector<double>,double> invert(std::vector<double> dc7unsrt);

std::vector<double> to_cell_param_redprim(std::vector<double> g6redprim);

bool agreement(std::vector<double> g6redprim, std::vector<double> invertedg6redprim);

bool filterNonXrays(std::vector<std::string> cell_params);

double unit_cell_vol(std::vector<double> cell_params);

#endif