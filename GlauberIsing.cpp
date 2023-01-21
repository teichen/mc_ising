/*
    Copyright (C) 2018  Paul E. Teichen

    This program executes Monte Carlo sampling of a
    three-dimensional Ising field.

    lambda tunes the nearest-neighbor interaction strength,
    restart=True signals if a restart file should be read.

    This program was funded by Adam P. Willard
*/
#include "GlauberIsing.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <gperftools/profiler.h>

#include <algorithm>

#include <sys/time.h>
#include <unistd.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
using std::ifstream;
#include <cstring>

const double PI = 3.1415926535897932384626433832795028841971693;

// Clocking prototype
long long
timeval_diff(struct timeval *difference,
             struct timeval *end_time,
             struct timeval *start_time
            )
{

  struct timeval temp_diff;

  if(difference==NULL)
  {
    difference=&temp_diff;
  }

  difference->tv_sec =end_time->tv_sec -start_time->tv_sec ;
  difference->tv_usec=end_time->tv_usec-start_time->tv_usec;

  while(difference->tv_usec<0)
  {
    difference->tv_usec+=1000000;
    difference->tv_sec -=1;
  }

  return 1000000LL*difference->tv_sec+
                   difference->tv_usec;

}

using namespace std;

GlauberIsing::GlauberIsing()
{
}

GlauberIsing::GlauberIsing(bool restart, bool logging, double lambda, int tsteps)
{
    /* Args:
               restart (bool) : if true, read in initial boolean field
               logging (bool) : if true, print Glauber dynamics
               lambda (double): nearest-neighbor interaction strength
               tsteps (int)   : number of sampling steps
    */
    ProfilerStart("/tmp/prof.out"); // memory profiler

    restart = restart;
    logging = logging;
    lambda  = lambda;
    tsteps  = tsteps;

    L   = lattice.L;   // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice
    nL  = pow(L,dim);

    twrite = mc.twrite; // write every twrite steps
    n_rn   = 10000; 

    mem_test = false;
    initarrays();

    srand (time(NULL)+lambda); // re-initialize random seed

    if(logging)
    {
        // output header
        cout << "## System Parameters ##" << endl;
        cout << "    L=" << L << ", dim=" << dim << endl;
        cout << "" << endl;
        cout << "## Monte Carlo Parameters ##" << endl;
        cout << "    tsteps=" << tsteps << endl;
        cout << "" << endl;

        cout << "lambda=" << lambda << endl;
        cout << "" << endl;
    }
    
    run(tsteps);
}

void GlauberIsing::run(int tsteps)
{
    /* master method to calculate Glauber dynamics (Monte Carlo sampling) using
       the Ising model
    */
    struct timeval earlier;
    struct timeval later;
    struct timeval interval;

    if(gettimeofday(&earlier,NULL)) // start clock
    {
        perror("third gettimeofday()");
        exit(1);
    }

    // initial Ising field
    if(restart)
    {
        // read in data
        read_set_field(n);
    }
    else
    {
        init_field(n);
    }

    int i;
    for (i=0; i<(int)(pow(L, dim)); i++)
    {
        n0[i] = n[i];
    }

    // evaluate initial energy
    eLG   = 0.0; 
    e0_LG = eLG;

    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        eLG   = model.get_energy(lambda, n, i);
        e0_LG = e0_LG + eLG;
    }

    eLG = e0_LG;

    if(logging)
    {
        // print Sampling trajectory energies
        cout << "## Monte Carlo Trajectory ##" << endl;
        cout << "\ttrial\tE(LG)" << endl;
        cout << "\t" << 0 << "\t" << eLG << endl;
    }

    int taccept;
    taccept = 0; // number of accepted moves

    int t;

    for (t=0; t<tsteps; t++)
    {
        glauber_sweep(lambda);

        if ( ((t+1) % twrite)==0)
        {
            if(logging)
            {
                write_field(n, "n.dat");
                cout << "\t" << t+1 << "\t" << eLG << endl;
            }
        }
    }

    eLG   = 0.0; 
    e0_LG = eLG;

    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        eLG   = model.get_energy(lambda, n, i);
        e0_LG = e0_LG + eLG;
    }

    eLG = e0_LG;

    if(logging)
    {
        cout << "## Monte Carlo Trajectory ##" << endl;
        cout << "\ttrial\tE(LG)" << endl;
        cout << "\t" << 0 << "\t" << eLG << endl;
    }

    if(logging)
    {
        // write restart files
        write_field(n, "n_restart.dat");

        std::ofstream edump;
        edump.open("energy_restart.dat", std::ios_base::trunc | std::ios_base::out);

        edump << eLG << "\n" << etot << endl;

        edump.close();
    }

    if(gettimeofday(&later,NULL)) // stop clock
    {
        perror("fourth gettimeofday()");
        exit(1);
    }

    timeval_diff(&interval,&later,&earlier);

    if(logging)
    {
        cout << "t = " << interval.tv_sec << "." <<
             interval.tv_usec << " sec" << endl;
        cout << " " << endl;

        cout << "Accepted " << taccept << " Monte Carlo moves " << endl;
        cout << " " << endl;
    }
    ProfilerStop();
}

void GlauberIsing::glauber_sweep(double lambda)
{
    /* full Glauber sweep
       Args:
               lambda (double): nearest-neighbor interaction strength
    */
    int i, j, k, ncell;

    for (i=0; i<nL; i++)
    {
        rn_flip[i] = (rand() % 1000) * 0.001;
    }

    // even x/z sweep, all y
    for (i=0; i<(int)(L/2); i++)
    {
        for (j=0; j<L; j++)
        {
            for (k=0; k<(int)(L/2); k++)
            {
                ncell = lattice.flatten_position(2*i, j, 2*k);
                glauber_flip(lambda, ncell);
            }
        }
    }

    // odd x/z sweep, all y
    for (i=0; i<(int)(L/2); i++)
    {
        for (j=0; j<L; j++)
        {
            for (k=0; k<(int)(L/2); k++)
            {
                ncell = lattice.flatten_position(2*i+1, j, 2*k+1);
                glauber_flip(lambda, ncell);
            }
        }
    }

    // odd x even z sweep, all y
    for (i=0; i<(int)(L/2); i++)
    {
        for (j=0; j<L; j++)
        {
            for (k=0; k<(int)(L/2); k++)
            {
                ncell = lattice.flatten_position(2*i+1, j, 2*k);
                glauber_flip(lambda, ncell);
            }
        }
    }
   
    // even x odd z sweep, all y
    for (i=0; i<(int)(L/2); i++)
    {
        for (j=0; j<L; j++)
        {
            for (k=0; k<(int)(L/2); k++)
            {
                ncell = lattice.flatten_position(2*i, j, 2*k+1);
                glauber_flip(lambda, ncell);
            }
        }
    }
}

void GlauberIsing::glauber_flip(double lambda, int i)
{
    /* attempt to flip a local value of the Ising field
       Args:
               lambda (double): nearest-neighbor interaction strength
               i (int)        : lattice site
    */
    int j,k;
    j = 0; k = 0;

    double de,etmp;
    de = 0.0; etmp = 0.0;

    double rn,boltz;
    rn    = 0.0;
    boltz = 0.0;

    eLG_i = model.get_energy(lambda, n, i);

    // add energy of neighboring cells
    int ri[3];
    lattice.unpack_position(i, ri);

    int nn[3];
    lattice.nearest_neighbors(ri, nn);

    for (j=0; j<dim; j++)
    {
        etmp  = model.get_energy(lambda, n, nn[j]);
        eLG_i = eLG_i + etmp;
    }

    if (n[i]==(-1))
    {
        n[i] = 1;
    }
    else
    {
        n[i] = -1;
    }

    eLG_f = 0.0;
    eLG_f = model.get_energy(lambda,n,i);

    for (j=0; j<dim; j++)
    {
        etmp  = model.get_energy(lambda, n, nn[j]);
        eLG_f = eLG_f + etmp;
    }

    eLG = eLG + (eLG_f-eLG_i);
    de  = eLG_f-eLG_i;

    // acceptance criteria

    rn = rn_flip[i];

    boltz = exp(-de); // temperature = 1.0
    // boltz = 1.0/(1.0+exp(de/temp));

    if (rn <= boltz)
    {
        n0[i] = n[i];
        e0_LG = eLG;
    }
    else
    {
        n[i] = n0[i];
        eLG  = e0_LG;
    }
}

void GlauberIsing::read_set_field(int* n)
{
    /* if restarting an earlier simulation, set the initial Ising field
       by reading from disk
       Args:
               n (int): Ising field
    */
    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    filenameStream << "./n_restart.dat";
    filename = filenameStream.str();

    fin.open(filename.c_str());

    int j = 0;
    while( getline(fin,line) )
    {
        n[j] = atoi (line.c_str());
        j++;
    }

    fin.close();
    filenameStream.str("");
}

void GlauberIsing::init_field(int* n)
{
    /* random initialization of an Ising field
       Args:
               n (int): Ising field
    */
    for (int j=0; j<(int)(pow(L,dim)); j++)
    {
        //n[j] = 1; // uniform
        n[j] = 2*(rand() % 2)-1;
    }
}

void GlauberIsing::write_field(int* n, string filename)
{
    /* write out values of a boolean field
       Args:
               n (int)          : Ising field
               filename (string): file name in which to dump Ising field values
    */
    std::ofstream ndump;

    if (filename.find("restart") != std::string::npos)
    {
        ndump.open(filename, std::ios_base::trunc | std::ios_base::out);
    }
    else
    {
        ndump.open(filename, std::ios_base::app);
    }
    
    for (int i=0; i<(int)(pow(L,dim)); i++)
    {
        ndump << n[i] << "\n";
    }

    ndump.close();
}

void GlauberIsing::initarrays()
{
    n = (int*) calloc (pow(L,dim), sizeof(int));
    n0 = (int*) calloc (pow(L,dim), sizeof(int));

    mem_test = true;
}

GlauberIsing::~GlauberIsing()
{
    if(mem_test==true)
    {
    delete [] n;
    delete [] n0;

    cout << "Deallocate GlauberIsing memory" << endl;

    }
}
