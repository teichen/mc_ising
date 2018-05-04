/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// SimSpace.h
#ifndef _SIMSPACE
#define _SIMSPACE

using namespace std;

class SimSpace
{

public:

    bool mem_test;

    SimSpace(); 

    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice
    int neighbors;

    ~SimSpace(); 

private:

};

#endif
