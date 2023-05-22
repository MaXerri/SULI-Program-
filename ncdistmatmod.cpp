#include "ncdist/G6.h"
#include "ncdist/Reducer.h"
#include "ncdist/Delone.h"
#include "ncdist/LRL_Cell.h"
#include "ncdist/LRL_Cell_Degrees.h"
#include "ncdist/NCDist.h"
#include "ncdist/S6M_SellingReduce.h"

#include "ncdistmatmod.h"

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

int donc = 0, doncsq = 0,
    dodc7 = 0, dodc7sq = 0, dodc7unsrt = 0, dodc7unsrtsq = 0,
    dodc10 = 0, dodc10sq = 0, info = 0;


/**
 * produced reduced primitive g6 representation from the cell parameters.  This code utilizes functions
 * declarations in github.com/yayahjb/ncdist and this function itself is a modified version of a function
 * defined in the file ncdist_mat.cpp 
 * @param testlattice: the bravis lattice type 
 * @param a: edge length a in angstroms 
 * @param b: edge length b in angstroms
 * @param c: edge length c in angstroms
 * @param alpha: angle alpha in degrees 
 * @param beta: angle beta in degrees 
 * @param gamma: angle gamma in degrees 
 * @return vector<string> with the g6 reduced primitive components 
 */
std::string makeprimredcellmod(std::string testlattice,
                               double a, double b, double c,
                               double alpha, double beta, double gamma, double extra)
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 redprimcell;
    G6 redprimcell_as_g6;
    D7 d7redprimcell;
    G6 d7redprimcell_as_g6;
    S6 s6redprimcell;
    G6 s6redprimcell_as_g6;
    double d7primcell[7];
    double dc13[13];
    double dc7[7];
    double dc7_unsrt[7];
    double s6primcell[6];
    double dprimcell[6];
    double dg6redprimcell[6];
    G6 dredprimcell;
    Mat66 mc;
    Mat66 m;
    Mat66 dm;
    G6 primcell;
    G6 recipcell;
    G6 reducedBase;
    G6 g6primredprobe_as_g6;
    D7 d7primredprobe_as_g6;
    S6 S6primredprobe_as_g6;
    LRL_Cell g6primredprobe;
    LRL_Cell d7primredprobe;
    LRL_Cell s6primredprobe;
    double crootvol;
    LRL_Cell rawcell(a, b, c, alpha, beta, gamma);
    int ii;
    bool ret;
    int reduced;
    if (testlattice.size() < 1)
    {
        latsym = "P";
    }
    else
    {
        latsym = testlattice.substr(0, 1);
    }
    clatsym = latsym[0];
    switch (clatsym)
    {
    case 'P':
    case 'p':
    case 'A':
    case 'a':
    case 'B':
    case 'b':
    case 'C':
    case 'c':
    case 'I':
    case 'i':
    case 'F':
    case 'f':
    case 'R':
    case 'r':
    case 'H':
    case 'h':
        CS6M_CelltoG6(rawcell, v6cell);
        CS6M_LatSymMat66(v6cell, clatsym, mc, primcell);
        break;
    case 'V':
    case 'v':
        primcell[0] = a;
        primcell[1] = b;
        primcell[2] = c;
        primcell[3] = alpha;
        primcell[4] = beta;
        primcell[5] = gamma;
        break;
    case 'D':
    case 'd':
        primcell[0] = a;
        primcell[1] = b;
        primcell[2] = c;
        primcell[3] = beta - b - c;
        primcell[4] = gamma - a - c;
        primcell[5] = extra - a - b;
        break;
    case 'S':
    case 's':
        primcell[3] = 2. * a;
        primcell[4] = 2. * b;
        primcell[5] = 2. * c;
        primcell[0] = -alpha - c - b;
        primcell[1] = -beta - c - a;
        primcell[2] = -gamma - b - a;
        break;
    default:
        std::cerr << "Unrecognized lattice symbol " << testlattice << " treated as P" << std::endl;
        latsym = "P";
        CS6M_CelltoG6(rawcell, v6cell);
        CS6M_LatSymMat66(v6cell, clatsym, mc, primcell);
        break;
    }
    dprimcell[0] = primcell[0];
    dprimcell[1] = primcell[1];
    dprimcell[2] = primcell[2];
    dprimcell[3] = primcell[3];
    dprimcell[4] = primcell[4];
    dprimcell[5] = primcell[5];
    reduced = 0;
    CS6M_G6Reduce(dprimcell, dg6redprimcell, reduced);

    if (!reduced)
    {
        for (ii = 0; ii < 6; ii++)
            redprimcell[ii] = redprimcell_as_g6[ii] = 0.;
    }
    else
    {
        for (ii = 0; ii < 6; ii++)
            redprimcell_as_g6[ii] = dg6redprimcell[ii];
    }
    CS6M_G6toD7(primcell, d7primcell);
    reduced = 0;
    CS6M_D7Reduce(d7primcell, d7redprimcell, reduced);
    if (!reduced)
    {
        for (ii = 0; ii < 6; ii++)
            d7redprimcell_as_g6[ii] = 0.;
        for (ii = 0; ii < 7; ii++)
            d7redprimcell[ii] = 0.;
    }
    else
    {
        CS6M_D7toG6(d7redprimcell, d7redprimcell_as_g6);
    }
    CS6M_G6toS6(primcell, s6primcell);
    reduced = 0;
    CS6M_S6Reduce(s6primcell, s6redprimcell, reduced);
    if (!reduced)
    {
        for (ii = 0; ii < 6; ii++)
            s6redprimcell[ii] = s6redprimcell_as_g6[ii] = 0;
    }
    else
    {
        CS6M_S6toG6(s6redprimcell, s6redprimcell_as_g6);
    }
    g6primredprobe = LRL_Cell_Degrees(redprimcell_as_g6);
    d7primredprobe = LRL_Cell_Degrees(d7redprimcell_as_g6);
    s6primredprobe = LRL_Cell_Degrees(s6redprimcell_as_g6);

    return ((std::to_string(dg6redprimcell[0])) + " " + (std::to_string(dg6redprimcell[1])) + " " +
            (std::to_string(dg6redprimcell[2])) + " " + (std::to_string(dg6redprimcell[3])) + " " +
            (std::to_string(dg6redprimcell[4])) + " " + (std::to_string(dg6redprimcell[5])));
}
