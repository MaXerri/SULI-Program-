#include "conversions.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <cmath>
#include <tuple> 

using namespace std;

/**
 * Converts to a vector of doubles
 * @param dc7unsrtdist_mat_out: string vector which is to be converted to vector<double>
 * @return vector<double>
 */
std::vector<double> toReducePrimG6(std::string dc7unsrtdist_mat_out)
{
    std::stringstream iss(dc7unsrtdist_mat_out);

    double number;
    std::vector<double> myNumbers;
    while (iss >> number)
    {
        myNumbers.push_back(number);
    }
    return myNumbers;
}

/**
 * Converts from g6 reduced primitive to dc7unsrt
 * @param redprimg6: vector<double> which has reduced primitive g6 compnenets
 * @return dc7unsrt vector
 */
std::vector<double> todc7unsrt(std::vector<double> redprimg6)
{
    std::vector<double> dc7unsrt(7);
    dc7unsrt[0] = redprimg6[0];
    dc7unsrt[1] = redprimg6[1];
    dc7unsrt[2] = redprimg6[2];
    dc7unsrt[3] = redprimg6[1] + redprimg6[2] - abs(redprimg6[3]);
    dc7unsrt[4] = redprimg6[0] + redprimg6[2] - abs(redprimg6[4]);
    dc7unsrt[5] = redprimg6[0] + redprimg6[1] - abs(redprimg6[5]);
    dc7unsrt[6] = min({redprimg6[0] + redprimg6[1] + redprimg6[2] + redprimg6[3] + redprimg6[4] + redprimg6[5],
                       redprimg6[0] + redprimg6[1] + redprimg6[2] + redprimg6[3] - redprimg6[4] - redprimg6[5],
                       redprimg6[0] + redprimg6[1] + redprimg6[2] - redprimg6[3] + redprimg6[4] - redprimg6[5],
                       redprimg6[0] + redprimg6[1] + redprimg6[2] - redprimg6[3] - redprimg6[4] + redprimg6[5]});
    return dc7unsrt;
}

/**
 * Converts from dc7unsrt back to g6
 * @param dc7unsrt: dc7unsrt vector
 * @return a tuple containing the inverted g6 values from the dc7unsrt and the tau value  
 */
std::tuple<vector<double>,double> invert(std::vector<double> dc7unsrt)
{
    std::vector<double> invertedg6unsigned(6);
    double tau;
    invertedg6unsigned[0] = dc7unsrt[0];
    invertedg6unsigned[1] = dc7unsrt[1];
    invertedg6unsigned[2] = dc7unsrt[2];
    invertedg6unsigned[3] = abs(dc7unsrt[1] + dc7unsrt[2] - dc7unsrt[3]);
    invertedg6unsigned[4] = abs(dc7unsrt[0] + dc7unsrt[2] - dc7unsrt[4]);
    invertedg6unsigned[5] = abs(dc7unsrt[0] + dc7unsrt[1] - dc7unsrt[5]);

    tau = invertedg6unsigned[0] + invertedg6unsigned[1] + invertedg6unsigned[2] -
          (invertedg6unsigned[3] + invertedg6unsigned[4] + invertedg6unsigned[5]);
    tau = round( tau * 1000000.0 ) / 1000000.0; // 6 decimal places
    if (tau == round(dc7unsrt[6] *1000000.0) / 1000000.0)
    {
        if (invertedg6unsigned[3] != 0)
        {
            invertedg6unsigned[3] = -1 * invertedg6unsigned[3];
        }
        if (invertedg6unsigned[4] != 0)
        {
            invertedg6unsigned[4] = -1 * invertedg6unsigned[4];
        }
        if (invertedg6unsigned[5] != 0)
        {
            invertedg6unsigned[5] = -1 * invertedg6unsigned[5];
        }
    }
    return {invertedg6unsigned,tau};
}

/**
 * makes the primitive reduced cell parameters
 * @param g6redprim: vector<double> which has reduced primitive g6 compnenets
 * @return reduced primitive cell parameters a b c alpha beta gamma
 */
std::vector<double> to_cell_param_redprim(std::vector<double> g6redprim)
{
    std::vector<double> cellParamRedPrim(6);
    cellParamRedPrim[0] = sqrt(g6redprim[0]);
    cellParamRedPrim[1] = sqrt(g6redprim[1]);
    cellParamRedPrim[2] = sqrt(g6redprim[2]);
    cellParamRedPrim[3] = acos(g6redprim[3] / (2 * cellParamRedPrim[1] * cellParamRedPrim[2])) * (180 / M_PI);
    cellParamRedPrim[4] = acos(g6redprim[4] / (2 * cellParamRedPrim[0] * cellParamRedPrim[2])) * (180 / M_PI);
    cellParamRedPrim[5] = acos(g6redprim[5] / (2 * cellParamRedPrim[0] * cellParamRedPrim[1])) * (180 / M_PI);

    return cellParamRedPrim;
}

/**
 * compares original g6 and inverted g6 parameters to the rounded 6th decimal point 
 * @param g6redprim: vector<double> which has original reduced primitive g6 compnenets
 * @param invertedg6redprim vector<double> which has inverted reduced primitive g6 compnenets
 * @return true if these 2 vectors are not equal and false otherwise 
 */
bool agreement(std::vector<double> g6redprim, std::vector<double> invertedg6redprim)
{
    double k = 1000000.0;

    //these next 2 statements are used to round the components of the 2 input vectors 
    transform(g6redprim.begin(), g6redprim.end(), g6redprim.begin(), [k](double &c)
              { return round(c * k) / k; });
    transform(invertedg6redprim.begin(), invertedg6redprim.end(), invertedg6redprim.begin(), [k](double &c)
              { return round(c * k) / k; });

    bool agree;
    if (g6redprim == invertedg6redprim)
    {
        agree = true;
    }
    else
    {
        agree = false;
    }
    return agree;
}

/**
 * compares original g6 and inverted g6 parameters to the rounded 6th decimal point 
 * @param cell_params: vector<string> which has original cell parameters
 * @return true if this is a non_xray structure and false otherwise 
 */
bool filterNonXrays(std::vector<std::string> cell_params)
{
    std::vector<string> not_an_xray = {"1.000", "1.000", "1.000", "90.000", "90.000", "90.000"};
    bool isNonXray;
    if (cell_params == not_an_xray)
    {
        isNonXray = true;
    }
    else
    {
        isNonXray = false;
    }
    return isNonXray;
}

/**
 * this uses the scalar triple product to calculate the general unit cell volume
 * @param cell_params: vector<double> which has original cell parameters
 * @return volume given the cell parameters
 */
double unit_cell_vol(std::vector<double> cell_params)
{
    //this uses the scalar triple product to calculate the general unit cell volume 
    double volume = cell_params[0]* cell_params[1] * cell_params[2] * sqrt(1+2*cos(cell_params[3] * M_PI / 180) \
    *cos(cell_params[4] * M_PI / 180)*cos(cell_params[5]* M_PI / 180)-pow(cos(cell_params[3]* M_PI / 180),2) \
    - pow(cos(cell_params[4]* M_PI / 180),2) -pow(cos(cell_params[5]* M_PI / 180),2));
    
    return(volume);
}