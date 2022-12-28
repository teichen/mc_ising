#include <iostream>
#include <cassert>
#include "../SimSpace.h"

using namespace std;
int main()
{
    SimSpace lattice; // define the coarse lattice
    int L; // length of lattice (number of sites)
    int dim; // dimensionality of lattice

    L = lattice.L; // length of lattice (number of sites)
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

    //assert(1==2); // run time assert
    //static_assert(1==2, "1 = 2"); // compile time assert
    //cout << "expression valid...execution continues.\n";
}
