/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// IsingModel.h
#ifndef _ISINGMODEL
#define _ISINGMODEL

#include "SimSpace.h"

#include <iostream>

using namespace std;

class IsingModel
{

public:

    bool mem_test;

    IsingModel(); 

    SimSpace lattice; // define the coarse lattice
    int L;            // length of lattice (number of sites)
    int dim;          // dimensionality of lattice

    double get_energy(double,int*,int); // Landau Ginzburg energy

    ~IsingModel(); 

private:

};

#endif


