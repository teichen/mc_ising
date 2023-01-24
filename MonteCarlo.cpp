/*
    Copyright (C) 2018  Paul E. Teichen

    The MonteCarlo class sets the Monte Carlo
    sampling parameters: total trial steps.

    This program was funded by Adam P. Willard
*/
#include "MonteCarlo.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;

MonteCarlo::MonteCarlo()
{
    // Pawley et al. 1984 MCRG benchmark
    //tsteps = 2e6;
    //twrite = 1e4;
    tsteps = 10;
    twrite = 10;

    //2.5e5 sweeps generate independent config. for L=64
    //5e5 sweeps for L=32
    //2e6 sweeps for L=16
    //4e6 sweeps for L=8

    srand (time(NULL)); // re-initialize random seed

    n_rn   = 10000; 
}

void MonteCarlo::randomize_sampling(int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        rn_flip[i] = (rand() % 1000) * 0.001;
    }
}

bool MonteCarlo::trial_flip(int i, double de)
{
    double rn, boltz;

    rn = rn_flip[i];

    boltz = exp(-de); // temperature = 1.0
    // boltz = 1.0/(1.0+exp(de/temp));

    if (rn <= boltz)
    {
        return true;
    }
    else
    {
        return false;
    }
}

MonteCarlo::~MonteCarlo()
{
}
