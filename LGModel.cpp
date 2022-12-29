/*
    Copyright (C) 2018  Paul E. Teichen

    The LGModel class calculates the Landau-Ginzburg
    energy (e.g. nearest-neighbor interactions).

    This program was funded by Adam P. Willard
*/
#include "LGModel.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

LGModel::LGModel()
{
    L    = lattice.L;   // length of lattice (number of sites)
    dim  = lattice.dim; // dimensionality of lattice
}

double LGModel::get_energy(double lambda, int* n, int ncell)
/* Landau Ginzburg energy
*/
{
    double eLG;
    eLG = 0.0;

    int ni,nj;
    ni = n[ncell];
    nj = 0;

    int* r;
    r = lattice.unpack_position(ncell);

    if (r[0] == (L-1))
    {
        nj = n[(int)(0*pow(L,2) + r[1]*L + r[2])];
    }
    else
    {
        nj = n[(int)((r[0] + 1)*pow(L,2) + r[1]*L + r[2])];
    }

    eLG = eLG + lambda*ni*nj;

    if (r[1] == (L-1))
    {
        nj = n[(int)(r[0]*pow(L,2) + 0*L + r[2])];
    }
    else
    {
        nj = n[(int)(r[0]*pow(L,2) + (r[1] + 1)*L + r[2])];
    }

    eLG = eLG + lambda*ni*nj;

    if (r[2] == (L-1))
    {
        nj = n[(int)(r[0]*pow(L,2) + r[1]*L + 0)];
    }
    else
    {
        nj = n[(int)(r[0]*pow(L,2) + r[1]*L + r[2]+1)];
    }

    eLG = eLG + lambda*ni*nj;

    return eLG;

}

LGModel::~LGModel()
{
}

