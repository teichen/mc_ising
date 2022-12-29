/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// SimSpace.h
#ifndef _SIMSPACE
#define _SIMSPACE

#include <iostream>

using namespace std;

class SimSpace
{

public:

    bool mem_test;

    SimSpace(); 

    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice
    int neighbors;

    double separation(double[3],double[3]);

    int* unpack_position(int);

    int flatten_position(int, int, int);

    int* nearest_neighbors(int*, int*);

    ~SimSpace(); 

private:

};

#endif
