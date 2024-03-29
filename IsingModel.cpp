/*
    Copyright (C) 2018  Paul E. Teichen

    The IsingModel class calculates the Landau-Ginzburg
    energy (e.g. nearest-neighbor interactions).

    This program was funded by Adam P. Willard
*/
#include "IsingModel.h"
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

IsingModel::IsingModel()
{
    L    = lattice.L;   // length of lattice (number of sites)
    dim  = lattice.dim; // dimensionality of lattice
}

double IsingModel::get_energy(double lambda, int* n, int ncell)
/* local Landau Ginzburg energy including nearest-neighbor interactions
   with the local lattice cell
   Args:
           lambda (double): nearest-neighbror interaction strength
           n (int*)       : binary field
           ncell (int)    : lattice cell
*/
{
    double eLG;
    eLG = 0.0;

    int ni,nj;
    ni = n[ncell];
    nj = 0;

    int r[3];
    lattice.unpack_position(ncell, r);

    int nn_vals[3];
    lattice.nearest_neighbor_values(n, r, nn_vals);

    int i;
    for (i=0; i<dim; i++)
    {
        eLG = eLG + lambda * ni * nn_vals[i];
    }
    
    return eLG;
}

IsingModel::~IsingModel()
{
}

