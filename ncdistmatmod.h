#pragma once
#ifndef NCDISTMATMOD_H
#define NCDISTMATMOD_H
#include <vector>
#include <string>

std::string makeprimredcellmod(std::string testlattice,
                               double a, double b, double c,
                               double alpha, double beta, double gamma, double extra);

#endif