/*
    Copyright (C) 2018  Paul E. Teichen

    This program constructs .xyz data files for visualization 
    in VMD from the Ising field represented on the 
    coarse-grained lattice.

    This program was funded by Adam P. Willard
*/
#include "Visualizer.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

#include <time.h>
#include <math.h>
#include <cmath>

const double PI = 3.1415926535897932384626433832795028841971693;

using namespace std;

Visualizer::Visualizer()
{
    L = lattice.L; // length of lattice (number of sites)

    dim = lattice.dim; // dimensionality of lattice

    tsteps = mc.tsteps; // number of MC attempts
    twrite = mc.twrite; // write every twrite steps

    mem_test = false;
    initarrays();

    int t;
    int i,j,k;

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    // read in data

    filenameStream << "./n_restart.dat";
    filename = filenameStream.str();

    fin.open(filename.c_str());

    j = 0;
    while( getline(fin,line) )
    {
        n[j] = atoi (line.c_str());
        j++;
    }

    fin.close();
    filenameStream.str("");

    std::ofstream vmddump;
    vmddump.open("traj.xyz", std::ios_base::app);

    // visualize in VMD (xy plane)

    t = 0;

    for (i=0; i<L; i++)
    {
        for (j=0; j<L; j++)
        {
            for (k=0; k<L; k++)
            {
                if (n[(int)(t*pow(L,dim)+(i*pow(L,2)+j*L+k))]==1)
                {
                    vmddump << "O\t" << i << "\t" << j << "\t" << k << "\n";
                }
                else
                {
                    vmddump << "N\t" << i << "\t" << j << "\t" << k << "\n";
                }
            }
        }
    }

    vmddump.close();

}

void Visualizer::initarrays()
{
    n = (int*) calloc (pow(L,dim), sizeof(int));

    mem_test = true;
}

Visualizer::~Visualizer()
{
    if(mem_test==true)
    {
    delete [] n;
    cout << "Deallocate Visualizer memory" << endl;

    }
}
