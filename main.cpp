/*
    Copyright (C) 2018  Paul E. Teichen

    This program executes Monte Carlo sampling of a
    three-dimensional Ising field.

    mu is an external potential,
    lambda tunes the nearest-neighbor interaction strength,
    restart=True signals if a restart file should be read.

    This program was funded by Adam P. Willard
*/
#include <iostream>
#include <fstream>
#include <stdlib.h>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>
#include <cstdlib>

#include "GlauberIsing.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout << "Program Running..." << endl;
    cout << "" << endl;

    if (argc != 5)
    {
        if (argv[0])
            std::cout << "Usage: " << argv[0] << " <restart> <logging> <lambda> <tsteps>" << '\n';
        else
            std::cout << "Usage: glauber_ising <restart> <logging> <lambda> <tsteps>" << '\n';
        
        exit(1);
    }

    std::stringstream restart_input(argv[1]);
    std::stringstream logging_input(argv[2]);
    std::stringstream lambda_input(argv[3]);
    std::stringstream tsteps_input(argv[4]);

    bool restart, logging;
    double lambda;
    int tsteps;

    if (!(restart_input >> restart))
        restart = false;
    if (!(logging_input >> logging))
        logging = false;
    if (!(lambda_input >> lambda))
        lambda = -1;
    if (!(tsteps_input >> tsteps))
        tsteps = -1;

    GlauberIsing glauber_ising(restart,logging,lambda,tsteps); 

    cout << "Program Exiting..." << endl;
    cout << "" << endl;

    return 0;
}
