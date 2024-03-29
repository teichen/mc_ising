/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// MonteCarlo.h
#ifndef _MONTECARLO
#define _MONTECARLO

#include <iostream>

using namespace std;

class MonteCarlo
{

public:

    bool mem_test;

    MonteCarlo(); // build function + input

    int tsteps; // number of MC attempts
    int twrite; // write every twrite steps

    void randomize_sampling(int);
    bool trial_flip(int, double);

    int n_rn;
    double rn_flip[32*32*32];

    void initarrays();

    ~MonteCarlo(); // destructor

private:

};

#endif
