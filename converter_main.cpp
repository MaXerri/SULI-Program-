#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iterator>
#include <tuple> 

#include "conversions.h"
#include "ncdistmatmod.h"

using namespace std;

void print_vec(std::vector<double> const &a)
{
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a.at(i) << ' ';
        std::cout.flush();
    }
    std::cout << std::endl;
}

void print_vec_string(std::vector<string> const &a)
{
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a.at(i) << ' ';
        std::cout.flush();
    }
    std::cout << std::endl;
}

int main()
{
    // opens input stream for desired file which is to be in the users current directory
    ifstream input;
    ofstream total_output_file;
    ofstream outlier_file;
    input.open("/home/mxerri/SULI/crystal.txt"); // read file; specify path if not in current directory
    total_output_file.open("totalOutput.csv");
    outlier_file.open("outlier.csv");
    int outlier_count = 0;

    std::string csv_header = "idcode, spg_type, spg, a, b, c, alpha, beta, gamma, aprim, bprim, cprim, alphaprim, betaprim, gammaprim, r, s, t, u, v, w, dc7_1, dc7_2, dc7_3, dc7_4, dc7_5, dc7_6, dc7_7, tau, r_reco, s_reco, t_reco, u_reco, v_reco, w_reco, original_vol, redprimcell_vol, agreement";
    total_output_file << csv_header << endl;
    
    if (!input.is_open())
    {
        std::cout << "error opening file";
    }
    else
    {
        string line;
        int i = 0;
        while(getline(input,line))
        {
            if (i >=4) // data entries start on line 4
            {
                // formats the string into a vector
                std::stringstream ss(line);
                std::istream_iterator<std::string> begin(ss);
                std::istream_iterator<std::string> end;
                std::vector<std::string> vstrings(begin, end);
                
                // breaks down parts of vector into spg, spg_type, idcode and cell parameters
                string idcode = vstrings[0];
                string spg_type = vstrings[8];
                string spg;
                switch (vstrings.size())
                {
                case 11:
                    spg = vstrings[9];
                    break;
                case 12:
                    spg = vstrings[9] + " " + vstrings[10];
                    break;
                case 13:
                    spg = vstrings[9]+ " " + vstrings[10] + " " + vstrings[11];
                    break;
                case 14:
                    spg = vstrings[9]+ " " + vstrings[10]+ " " + vstrings[11] + " " + vstrings[12];
                    break;
                }
                //cell_params as a string vetcor and double vector
                std::vector<string> cell_params(vstrings.begin() + 2, vstrings.begin() + 8);
                std::vector<double> cell_param_double = {stof(cell_params[0]),stof(cell_params[1]), \
                stof(cell_params[2]),stof(cell_params[3]),stof(cell_params[4]),stof(cell_params[5])};

                //calls my modified ncdist function (ncdistmatmod.cpp) to output the g6 components as vector<string> type
                string dc7unsrtdist_mat_out = makeprimredcellmod(spg_type, stof(cell_params[0]),
                                                                 stof(cell_params[1]), stof(cell_params[2]),
                                                                 stof(cell_params[3]), stof(cell_params[4]),
                                                                 stof(cell_params[5]), 0.0);
                //converts dc7unsrt_mat_out vector<string> to vector<double> 
                vector<double> redprim_double = toReducePrimG6(dc7unsrtdist_mat_out);

                //these next lines call functions that perform the inversions between g6, dc7unsrt and inverted g6
                vector<double> cell_param_redprim = to_cell_param_redprim(redprim_double);
                vector<double> dc7unsrt(7);
                vector<double> invertedg6(6);
                double tau;
                dc7unsrt = todc7unsrt(redprim_double);
                tie(invertedg6,tau) = invert(dc7unsrt);

                //calls functions which check if the inverted g6 agrees with the original g6
                bool agree = agreement(redprim_double, invertedg6);
                //calls function which removes non-xray structures 
                bool p = filterNonXrays(cell_params);
                //Call functions to calculate volume:
                double redprim_cell_vol = unit_cell_vol(cell_param_redprim);
                double original_cell_vol = unit_cell_vol(cell_param_double);
                
                //comma seperated string to add to a csv 
                string complete_string = idcode + ", " + spg_type + ", " + spg + ", " + cell_params[0] + ", " 
                                        + cell_params[1] + ", " + cell_params[2] + ", " + cell_params[3] + ", " + cell_params[4] 
                                        + ", " + cell_params[5] + ", " + to_string(cell_param_redprim[0]) + ", " + to_string(cell_param_redprim[1]) 
                                        + ", " + to_string(cell_param_redprim[2]) + ", " + to_string(cell_param_redprim[3]) + ", " 
                                        + to_string(cell_param_redprim[4]) + ", " + to_string(cell_param_redprim[5]) + ", "
                                        + to_string(redprim_double[0]) + ", " + to_string(redprim_double[1]) + ", " + to_string(redprim_double[2]) 
                                        + ", " + to_string(redprim_double[3]) + ", " + to_string(redprim_double[4]) + ", " + to_string(redprim_double[5]) 
                                        + ", " + to_string(dc7unsrt[0]) + ", " + to_string(dc7unsrt[1]) +", " + to_string(dc7unsrt[2]) 
                                        + ", " + to_string(dc7unsrt[3]) + ", " + to_string(dc7unsrt[4]) + ", " + to_string(dc7unsrt[5]) 
                                        + ", " + to_string(dc7unsrt[6]) + ", " + to_string(tau) + ", " + to_string(invertedg6[0]) + ", " + to_string(invertedg6[1])
                                        + ", " + to_string(invertedg6[2]) + ", " + to_string(invertedg6[3]) +", " + to_string(invertedg6[4])
                                        + ", " + to_string(invertedg6[5]) + ", " + to_string(original_cell_vol) + ", " + to_string(redprim_cell_vol) 
                                        + ", "+  to_string(agree);

                //checks if the cell is an outlier and fills respective csv file with the above string 
                if (p == false){
                    total_output_file << complete_string << endl;

                    if (agree == false){
                        if (outlier_count == 0)
                        {
                            outlier_file << csv_header << endl;
                            outlier_count = outlier_count + 1;
                        }
                        else
                        {
                            outlier_count+= 1;
                        }
                        outlier_file << complete_string <<endl;
                    }
                }
            }
            i += 1;
        }
    }
    //general output explanation 
    if (outlier_count ==0)
    {
        std::cout << "There were no outliers found in this dataset" << endl;
    }
    else
    {
        std::cout << "There was " << outlier_count << " outlier(s) and they can be found in the output.csv file"<<endl; 
    }
    std::cout << "A .csv file of the conversion process (totalOutput.csv) was outoutted into your current directory" <<endl;
    input.close();
    total_output_file.close();
    outlier_file.close();
    return 0;
}