/*
    Copyright (C) 2018  Paul E. Teichen

    The SimSpace class sets the dimensions of the
    simulation volume.

    This program was funded by Adam P. Willard
*/
#include "SimSpace.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>

using namespace std;

SimSpace::SimSpace()
{
    //L = 8;
    L = 16;
    //L = 32;
    dim = 3;

}

SimSpace::~SimSpace()
{
}
