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

    double separation(double[3],double[3]);

    void unpack_position(int, int[3]);

    int flatten_position(int, int, int);

    void nearest_neighbors(int*, int[3]);
    void nearest_neighbor_values(int*, int*, int[3]);

    ~SimSpace(); 

private:

};

#endif
