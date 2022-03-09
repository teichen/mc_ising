/*
    Copyright (C) 2018  Paul E. Teichen

    This program was funded by Adam P. Willard
*/
// LGModel.h
#ifndef _LGMODEL
#define _LGMODEL

#include "SimSpace.h"

#include <iostream>

using namespace std;

class LGModel
{

public:

    bool mem_test;

    LGModel(); 

    SimSpace lattice; // define the coarse lattice
    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    double temp;

    double get_energy(double,int*,int,double&); // Landau Ginzburg energy

    ~LGModel(); 

private:

};

#endif


