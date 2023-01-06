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

    tsteps = 1000;
    glauber_ising.run(tsteps);

    int i, j;
    int nsum = 0;
    
    int ri[3];
    int nn_vals[3];

    for (i=0; i<(int)(glauber_ising.nL); i++)
    {
        glauber_ising.lattice.unpack_position(i, ri);
        glauber_ising.lattice.nearest_neighbor_values(glauber_ising.n, ri, nn_vals);

        for (j=0; j<(int)(glauber_ising.dim); j++)
        {
            nsum += glauber_ising.n[i] * nn_vals[j];
        }
    }
    assert(nsum < 100);
    assert(nsum > -100);

    glauber_ising.lambda = -10.0; // pressure to reduce interface
    glauber_ising.run(tsteps);

    nsum = 0;
    for (i=0; i<(int)(glauber_ising.nL); i++)
    {
        glauber_ising.lattice.unpack_position(i, ri);
        glauber_ising.lattice.nearest_neighbor_values(glauber_ising.n, ri, nn_vals);

        for (j=0; j<(int)(glauber_ising.dim); j++)
        {
            nsum += glauber_ising.n[i] * nn_vals[j];
        }
    }
    assert(nsum > 500);

    glauber_ising.lambda = 10.0; // pressure to increase interface
    glauber_ising.run(tsteps);

    nsum = 0;
    for (i=0; i<(int)(glauber_ising.nL); i++)
    {
        glauber_ising.lattice.unpack_position(i, ri);
        glauber_ising.lattice.nearest_neighbor_values(glauber_ising.n, ri, nn_vals);

        for (j=0; j<(int)(glauber_ising.dim); j++)
        {
            nsum += glauber_ising.n[i] * nn_vals[j];
        }
    }
    assert(nsum < -500);

    return 0;
}
