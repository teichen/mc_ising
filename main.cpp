/*
    Copyright (C) 2018  Paul E. Teichen

    This program executes Monte Carlo sampling of a
    three-dimensional Ising field.

    mu is an external potential,
    lambda tunes the nearest-neighbor interaction strength,
    run_flag=1 signals if a restart file should be read.

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

#include "ising.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout << "Program Running..." << endl;
    cout << "" << endl;

    if (argc != 4)
    {
        if (argv[0])
            std::cout << "Usage: " << argv[0] << " <run_flag> <mu> <lambda>" << '\n';
        else
            std::cout << "Usage: isingSim <run_flag> <mu> <lambda>" << '\n';
        
        exit(1);
    }

    std::stringstream convert1(argv[1]);
    std::stringstream convert2(argv[2]);
    std::stringstream convert3(argv[3]);

    int run_flag;
    double mu,lambda;

    if (!(convert1 >> run_flag))
        run_flag = -1;
    if (!(convert2 >> mu))
        mu = -1;
    if (!(convert3 >> lambda))
        lambda = -1;

    ising mycalcs(run_flag,mu,lambda); 

    cout << "Program Exiting..." << endl;
    cout << "" << endl;

    return 0;
}
