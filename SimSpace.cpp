/*
    Copyright (C) 2018  Paul E. Teichen

    The SimSpace class sets the dimensions of the
    simulation volume.

    This program was funded by Adam P. Willard
*/
#include "SimSpace.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <iostream>

using namespace std;

SimSpace::SimSpace()
{
    /* Simulation volume
    */
    L   = 8;
    dim = 3;
}

double SimSpace::separation(double r1[3], double r2[3])
{
    /* calculate separation between two positions subject to
       periodic boundary conditions
       Args:
               r1 (double[3]): first real-space position
               r2 (double[3]): second real-space position
    */
    double d[3];
    d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;

    double r;
    int dimR;

    for (dimR=0; dimR<dim; dimR++)
    {
        if ((r1[dimR]-r2[dimR])>(double)(L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] - L;
        }
        else if((r1[dimR]-r2[dimR])<(double)(-L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] + L;
        }
        else
        {
            d[dimR] = r1[dimR]-r2[dimR];
        }
    }

    r = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));

    return r;
}

void SimSpace::unpack_position(int n, int r[3])
{
    /* unpack a flattened real-space lattice position
       Args:
               n (int):    flattened position, lattice cell
               r (int[3]): unpacked position
    */
    r[2] = n % L;
    r[1] = (int)((n - r[2]) / L) % L;
    r[0] = (int)((n - r[2] - r[1] * L) / pow(L, 2)) % L;
}

int SimSpace::flatten_position(int i, int j, int k)
{
    /* flatten position, position n -> (x*L^2+y*L+z)
       Args:
               i (int): x cell number
               j (int): y cell number
               k (int): z cell number
    */
    int n;
    n = (int)(i * pow(L, 2) + j * L + k);

    return n;
}

void SimSpace::nearest_neighbors(int* r, int nn[3])
{
    /* find all nearest neighbors for an input lattice position
       Args:
               r (int*)   : lattice position
               nn (int[3]): list of flattened positions for nearest-neighbors
    */
    if (r[0] == (L-1))
    {
        nn[0] =(int)(0*pow(L,2) + r[1]*L + r[2]);
    }
    else
    {
        nn[0] = (int)((r[0] + 1)*pow(L,2) + r[1]*L + r[2]);
    }
    if (r[1] == (L-1))
    {
        nn[1] = (int)(r[0]*pow(L,2) + 0*L + r[2]);
    }
    else
    {
        nn[1] = (int)(r[0]*pow(L,2) + (r[1] + 1)*L + r[2]);
    }
    if (r[2] == (L-1))
    {
        nn[2] = (int)(r[0]*pow(L,2) + r[1]*L + 0);
    }
    else
    {
        nn[2] = (int)(r[0]*pow(L,2) + r[1]*L + r[2]+1);
    }
}

void SimSpace::nearest_neighbor_values(int* n, int* nn, int nn_vals[3])
{
    /* extract field values for an input array of flattened nearest-neighbor locations
       Args:
               n (int*)        : binary field
               nn (int*)       : nearest neighbor locations
               nn_vals (int[3]): output values
    */
    nn_vals[0] = n[nn[0]];
    nn_vals[1] = n[nn[1]];
    nn_vals[2] = n[nn[2]];
}

SimSpace::~SimSpace()
{
}
