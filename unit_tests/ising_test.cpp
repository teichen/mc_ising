#include <iostream>
#include <cassert>
#include <math.h>
#include "../SimSpace.h"
#include "../IsingModel.h"

using namespace std;
int main()
{
    // unit testing with run time asserts
    // lattice tests:
    SimSpace lattice; // define the coarse lattice
    int L;            // length of lattice (number of sites)
    int dim;          // dimensionality of lattice

    L = lattice.L;     // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice

    assert(dim==3); // only supporting d=3

    double r1[3], r2[3];
    r1[0] = 0.0; r1[1] = 0.0; r1[2] = 0.0;
    r2[0] = 0.0; r2[1] = 0.0; r2[2] = 0.0;

    double d;
    d = lattice.separation(r1, r2);

    assert(d==0.0);

    r2[0] = 1.0;
    d = lattice.separation(r1, r2);
    assert(d==1.0);

    r2[0] = 0.0;
    r1[0] = L-1;
    d = lattice.separation(r1, r2);
    assert(d==1.0);

    // model tests:
    IsingModel model;

    double lambda;
    lambda = 1.0;

    double eLG,etot;
    eLG = 0.0; etot = 0.0;

    int n[(int)(pow(L,dim))];

    int i;
    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        n[i] = 1;
    }
    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        // calculate Landau Ginzburg and Gaussian energies only
        eLG = model.get_energy(lambda,n,i);

        etot = etot + eLG;

    }
    assert(etot==(double)(pow(L,dim) * dim * lambda));

    return 0;
}
