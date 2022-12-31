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

    restart = 0;
    logging = 0;
    lambda  = 0.0;

    GlauberIsing glauber_ising(restart,logging,lambda); 

    int nsum = 0;
    for (int i=0; i<(int)(glauber_ising.nL); i++)
    {
        nsum += glauber_ising.n[i];
    }
    cout << nsum << endl;
    // assert(dim==3); // only supporting d=3

    return 0;
}
