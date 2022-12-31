/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// GlauberIsing.h
#ifndef _GLAUBERISING
#define _GLAUBERISING

#include "SimSpace.h"
#include "IsingModel.h"
#include "MonteCarlo.h"

#include <iostream>

using namespace std;

class GlauberIsing
{

public:

    bool mem_test;

    GlauberIsing();
    GlauberIsing(bool,bool,double);

    // ----------------------------------------------- //

    SimSpace lattice; // define the coarse lattice
    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    double temp;

    MonteCarlo mc; // Monte Carlo parameterization
    int tsteps; // number of MC attempts
    int twrite; // write every twrite steps

    IsingModel model;

    double eLG,etot;

    int nL;

    int* n;
    int* n0;

    double e0,e0_LG;

    int n_rn;
    double rnvec[10000];
    double rn_flip[32*32*32];

    void initarrays();
    void glauber_flip(double,int);

    double eLG_i,etot_i,eLG_f,etot_f;

    void read_set_field(int*);
    void write_field(int*,string);

    ~GlauberIsing(); 

private:

};

#endif

