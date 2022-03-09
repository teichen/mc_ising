/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// Visualizer.h
#ifndef _VISUALIZER
#define _VISUALIZER

#include "SimSpace.h"
#include "MonteCarlo.h"

#include <iostream>

using namespace std;

class Visualizer
{

public:

    bool mem_test;

    Visualizer();
    Visualizer(int);

    SimSpace lattice; // define the coarse lattice
    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    MonteCarlo mc; // Monte Carlo parameterization
    int tsteps; // number of MC attempts
    int twrite; // write every twrite steps

    int* n;

    void initarrays();

    ~Visualizer(); // destructor

private:

};

#endif

