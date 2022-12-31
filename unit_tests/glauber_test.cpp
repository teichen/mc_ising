#include <iostream>
#include <cassert>
using std::cerr;
using std::cout;
using std::endl;
#include <math.h>
#include "../GlauberIsing.h"

using namespace std;
int main()
{
    // unit testing with run time asserts
    // ising tests:
    bool restart, logging;
    double lambda;
    int tsteps;

    restart = 0;
    logging = 0;
    lambda  = 0.0;
    tsteps  = 1;

    GlauberIsing glauber_ising(restart, logging, lambda, tsteps); 

    tsteps = 100;
    glauber_ising.run(tsteps);

    //glauber_ising.lambda = 5.0;
    //cout << glauber_ising.lambda << endl;
    //cout << glauber_ising.tsteps << endl;
    //glauber_ising.tsteps = 2;
    //cout << glauber_ising.tsteps << endl;

    int i;
    int nsum = 0;
    for (i=0; i<(int)(glauber_ising.nL); i++)
    {
        nsum += glauber_ising.n[i];
    }
    cout << nsum << endl;
    // assert(dim==3); // only supporting d=3

    glauber_ising.lambda = 0.01;
    glauber_ising.run(tsteps);

    nsum = 0;
    for (i=0; i<(int)(glauber_ising.nL); i++)
    {
        nsum += glauber_ising.n[i];
    }
    cout << nsum << endl;

    return 0;
}
