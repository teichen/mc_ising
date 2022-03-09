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
    L = lattice.L; // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice

    temp = 1.0; // temperature
}

double LGModel::get_energy(double lambda, int* n, int ncell, double& etot) // Landau Ginzburg energy
{
    double mu;
    mu = 0.0;

    double eLG,e2,emu;
    eLG = 0.0; e2 = 0.0; emu = 0.0;

    int ri[3];
    int rj[3];
    ri[0] = 0; ri[1] = 0; ri[2] = 0;
    rj[0] = 0; rj[1] = 0; rj[2] = 0;

    double ri_c[3];
    double rj_c[3];

    ri_c[0] = 0.0; ri_c[1] = 0.0; ri_c[2] = 0.0;
    rj_c[0] = 0.0; rj_c[1] = 0.0; rj_c[2] = 0.0;

    double r;
    r = 0.0;

    int dimR;
    dimR = 0;

    double d[3];
    d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;

    int nL = (int)pow(L,dim);

    int n0,n1; 
    n0 = 0; 
    n1 = 0; 

    int nn_counts;
    nn_counts = 0;

    int ni,nj;
    ni = n[ncell];
    nj = 0;

    int i,j,k;
    int jj,kk;

    i = 0; j = 0; k = 0;
    jj = 0; kk = 0;

    emu = mu*ni; 
    eLG = mu*ni;

    // position n -> (x*L^2+y*L+z)
    // z = i % L;
    // y = (i-z) % pow(L,2);
    // x = (i-z-y*L) % pow(L,dim);

    k = ncell % L;
    j = (int)((ncell-k)/L) % L;
    i = (int)((ncell-k-j*L)/pow(L,2)) % L;

    if (i==(L-1))
    {
        nj = n[(int)(0*pow(L,2)+j*L+k)];
    }
    else
    {
        nj = n[(int)((i+1)*pow(L,2)+j*L+k)];
    }

    e2 = lambda*ni*nj; 
    eLG = eLG + lambda*ni*nj;

    if (j==(L-1))
    {
        nj = n[(int)(i*pow(L,2)+0*L+k)];
    }
    else
    {
        nj = n[(int)(i*pow(L,2)+(j+1)*L+k)];
    }

    e2 = e2 + lambda*ni*nj; 
    eLG = eLG + lambda*ni*nj;

    if (k==(L-1))
    {
        nj = n[(int)(i*pow(L,2)+j*L+0)];
    }
    else
    {
        nj = n[(int)(i*pow(L,2)+j*L+k+1)];
    }

    e2 = e2 + lambda*ni*nj; 
    eLG = eLG + lambda*ni*nj;

    etot = eLG;

    return eLG;

}

LGModel::~LGModel()
{
}

