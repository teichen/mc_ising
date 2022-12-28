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

SimSpace::~SimSpace()
{
}
