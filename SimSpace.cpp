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
    L = 8;
    dim = 3;
}

double SimSpace::separation(double r1[3], double r2[3])
{
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

int* SimSpace::unpack_position(int n)
{
    // z = i % L;
    // y = (i-z) % pow(L,2);
    // x = (i-z-y*L) % pow(L,dim);
    static int r[3];

    r[2] = n % L;
    r[1] = (int)((n - r[2]) / L) % L;
    r[0] = (int)((n - r[2] - r[1] * L) / pow(L, 2)) % L;
    
    return r;
}

int SimSpace::flatten_position(int i, int j, int k)
{
    // position n -> (x*L^2+y*L+z)
    int n;
    n = (int)(i * pow(L, 2) + j * L + k);

    return n;
}

int* SimSpace::nearest_neighbors(int* n, int* r)
{
    static int nn[3];

    if (r[0] == (L-1))
    {
        nn[0] = n[(int)(0*pow(L,2) + r[1]*L + r[2])];
    }
    else
    {
        nn[0] = n[(int)((r[0] + 1)*pow(L,2) + r[1]*L + r[2])];
    }
    if (r[1] == (L-1))
    {
        nn[1] = n[(int)(r[0]*pow(L,2) + 0*L + r[2])];
    }
    else
    {
        nn[1] = n[(int)(r[0]*pow(L,2) + (r[1] + 1)*L + r[2])];
    }
    if (r[2] == (L-1))
    {
        nn[2] = n[(int)(r[0]*pow(L,2) + r[1]*L + 0)];
    }
    else
    {
        nn[2] = n[(int)(r[0]*pow(L,2) + r[1]*L + r[2]+1)];
    }

    return nn;
}


SimSpace::~SimSpace()
{
}
