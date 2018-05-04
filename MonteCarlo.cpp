/*
    Copyright (C) 2018  Paul E. Teichen

    The MonteCarlo class sets the Monte Carlo
    sampling parameters: total trial steps.

    This program was funded by Adam P. Willard
*/
#include "MonteCarlo.h"

using namespace std;

MonteCarlo::MonteCarlo()
{
    // Pawley et al. 1984 MCRG benchmark
    tsteps = 2e6;
    twrite = 1e4;

    //2.5e5 sweeps generate independent config. for L=64
    //5e5 sweeps for L=32
    //2e6 sweeps for L=16
    //4e6 sweeps for L=8
}

MonteCarlo::~MonteCarlo()
{
}
