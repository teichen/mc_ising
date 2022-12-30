#include <iostream>
#include <cassert>
using std::cerr;
using std::cout;
using std::endl;
#include <math.h>
#include "../ising.h"

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

    ising glauber_calcs(restart,logging,lambda); 

    int nsum = 0;
    for (int i=0; i<(int)(glauber_calcs.nL); i++)
    {
        nsum += glauber_calcs.n[i];
    }
    cout << nsum << endl;
    // assert(dim==3); // only supporting d=3

    return 0;
}
