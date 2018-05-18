/*
    Copyright (C) 2018  Paul E. Teichen

    This program executes Monte Carlo sampling of a
    three-dimensional Ising field.

    mu is an external potential,
    lambda tunes the nearest-neighbor interaction strength,
    run_flag=1 signals if a restart file should be read.

    This program was funded by Adam P. Willard
*/
#include "ising.h"
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>

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

#include <omp.h>

const double PI = 3.1415926535897932384626433832795028841971693;

//-----------------------------------------------------------------//
// Clocking prototype
//-----------------------------------------------------------------//

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

  /* Using while instead of if below makes the code slightly more robust. */

  while(difference->tv_usec<0)
  {
    difference->tv_usec+=1000000;
    difference->tv_sec -=1;
  }

  return 1000000LL*difference->tv_sec+
                   difference->tv_usec;

} /* timeval_diff() */

using namespace std;

ising::ising()
{
}

ising::ising(int run_flag, double mu, double lambda)
{
    L = lattice.L; // length of lattice (number of sites)
    dim = lattice.dim; // dimensionality of lattice
    nL = pow(L,dim);

    temp = 1.0; // temperature

    tsteps = mc.tsteps; // number of MC attempts
    twrite = mc.twrite; // write every twrite steps

    n_rn = 10000; 

    // ----------------------------------------------- //

    mem_test = false;
    initarrays();

    struct timeval earlier;
    struct timeval later;
    struct timeval interval;

    if(gettimeofday(&earlier,NULL)) //-----------Start clock ----------//
    {
        perror("third gettimeofday()");
        exit(1);
    }

    // lattice position working arrays
    double ri[3];
    double rj[3];

    ri[0] = 0.0; ri[1] = 0.0; ri[2] = 0.0;
    rj[0] = 0.0; rj[1] = 0.0; rj[2] = 0.0;

    // distance between lattice cells
    double r;
    r = 0.0;

    int dimR;
    dimR = 0;

    int i,j,k;
    i = 0; j = 0; k = 0;

    int di,dj,dk;
    di = 0; dj = 0; dk = 0;

    int ncell;
    ncell = 0;

    srand (time(NULL)+lambda); // re-initialize random seed

    // output header
    cout << "## System Parameters ##" << endl;
    cout << "    L=" << L << ", dim=" << dim << endl;
    cout << "" << endl;
    cout << "## Monte Carlo Parameters ##" << endl;
    cout << "    tsteps=" << tsteps << endl;
    cout << "" << endl;

    cout << "lambda=" << lambda << endl;
    cout << "" << endl;

    std::ostringstream filenameStream;
    std:string filename;

    ifstream fin;
    string line;

    std::ofstream ndump;

    // initial Ising field

    if (run_flag==0)
    {
        for (j=0; j<(int)(pow(L,dim)); j++)
        {
            //n[j] = 1;
            n[j] = 2*(rand() % 2)-1;
        }
    }
    else 
    {
        // read in data

        filenameStream << "./n_restart" << run_flag-1 << ".dat";
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
    }

    for (j=0; j<(int)(pow(L,dim)); j++)
    {
        n0[j] = n[j];
    }

    // evaluate initial energy piece-wise by volume

    etot = 0.0;
    eLG = 0.0; 

    e0 = etot;
    e0_LG = eLG;

    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        // calculate Landau Ginzburg and Gaussian energies only
        eLG = model.get_energy(lambda,n,i,etot);

        e0_LG = e0_LG + eLG;

    }

    eLG = e0_LG;

    etot = eLG;
    e0 = etot;

    nchange = 0;

    // print Sampling trajectory energies

    cout << "## Monte Carlo Trajectory ##" << endl;
    cout << "    trial nchange E(LG) E(TOT)" << endl;
    cout << "    " << 0 << "    " << nchange << "    " << 
                      eLG << "    " << etot << endl;


    int taccept;

    taccept = 0; // number of accepted moves

    e0 = etot;
    e0_LG = eLG;

    int t;
    int c;
    int cc;

    double de,etmp;
    de = 0.0; etmp = 0.0;

    double rn;

    int nvec[n_rn];

    double boltz;

    for (i=0; i<n_rn; i++)
    {
        //cvec[i] = rand() % nvols;
        rnvec[i] = (rand() % 1000)*0.001;
        nvec[i] = rand() % (int)(pow(L,dim));

    }

    for (t=0; t<tsteps; t++)
    {
        // full Glauber sweep
        for (i=0; i<nL; i++)
        {
            rn_flip[i] = (rand() % 1000)*0.001;
        }

        nchange = 0;

        for (i=0; i<(int)(L/2); i++)
        {
            for (j=0; j<L; j++)
            {
                for (k=0; k<(int)(L/2); k++)
                {
                    ncell = (2*i)*pow(L,2)+j*L+(2*k);
                    glauber_flip(lambda,ncell);
                }
            }
        }
        for (i=0; i<(int)(L/2); i++)
        {
            for (j=0; j<L; j++)
            {
                for (k=0; k<(int)(L/2); k++)
                {
                    ncell = (2*i+1)*pow(L,2)+j*L+(2*k+1);
                    glauber_flip(lambda,ncell);
                }
            }
        }
        for (i=0; i<(int)(L/2); i++)
        {
            for (j=0; j<L; j++)
            {
                for (k=0; k<(int)(L/2); k++)
                {
                    ncell = (2*i+1)*pow(L,2)+j*L+(2*k);
                    glauber_flip(lambda,ncell);
                }
            }
        }
        for (i=0; i<(int)(L/2); i++)
        {
            for (j=0; j<L; j++)
            {
                for (k=0; k<(int)(L/2); k++)
                {
                    ncell = (2*i)*pow(L,2)+j*L+(2*k+1);
                    glauber_flip(lambda,ncell);
                }
            }
        }

        if ( ((t+1) % twrite)==0)
        {
            ndump.open("n.dat", std::ios_base::app);

            for (i=0; i<(int)(pow(L,dim)); i++)
            {
                ndump << n[i] << "\n";
            }

            ndump.close();

            cout << "    " << t+1 << "    " << nchange << "    " << 
                              eLG << "    " << etot << endl;

        }
    }

    etot = 0.0;
    eLG = 0.0; 

    e0 = etot;
    e0_LG = eLG;

    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        // calculate Landau Ginzburg and Gaussian energies only
        eLG = model.get_energy(lambda,n,i,etot);

        e0_LG = e0_LG + eLG;

    }

    eLG = e0_LG;

    etot = eLG;
    e0 = etot;

    cout << "## Monte Carlo Trajectory ##" << endl;
    cout << "    trial E(LG) E(TOT)" << endl;
    cout << "    " << 0 << "    " <<
                      eLG << "    " << etot << endl;

    // write restart files
    ndump.open("n_restart.dat", std::ios_base::trunc | std::ios_base::out);

    for (i=0; i<(int)(pow(L,dim)); i++)
    {
        ndump << n[i] << "\n";
    }

    ndump.close();

    std::ofstream edump;
    edump.open("energy_restart.dat", std::ios_base::trunc | std::ios_base::out);

    edump << eLG << "\n" << etot << endl;

    edump.close();

    if(gettimeofday(&later,NULL)) //------------Stop clock --------------//
    {
        perror("fourth gettimeofday()");
        exit(1);
    }

    timeval_diff(&interval,&later,&earlier);

    cout << "t = " << interval.tv_sec << "." <<
         interval.tv_usec << " sec" << endl;
    cout << " " << endl;

    cout << "Accepted " << taccept << " Monte Carlo moves " << endl;
    cout << " " << endl;

}

void ising::glauber_flip(double lambda, int i)
{
    int j,k;
    j = 0; k = 0;

    double de,etmp;
    de = 0.0; etmp = 0.0;

    double ri[3];
    double rj[3];

    ri[0] = 0.0; ri[1] = 0.0; ri[2] = 0.0;
    rj[0] = 0.0; rj[1] = 0.0; rj[2] = 0.0;

    int cc;
    cc = 0;

    double r;
    r = 0.0;

    double rn,boltz;
    rn = 0.0;
    boltz = 0.0;

    eLG_i = 0.0;

    // calculate Landau Ginzburg and Gaussian energies only
    eLG_i = model.get_energy(lambda,n,i,etot_i);

    // add energy of backward neighboring cells

    // position n -> (x*L^2+y*L+z)
    // z = i % L;
    // y = (i-z) % pow(L,2);
    // x = (i-z-y*L) % pow(L,dim);

    ri[2] = (double)(i % L);
    ri[1] = (double)(((i-(i % L))/L) % L);
    ri[0] = (double)((int)((i-(i % L)-(((i-(i % L))/L) % L)*L)/pow(L,2)) % L);

    int ncell;

    if (ri[0]==0)
    {
        ncell = (int)((L-1)*pow(L,2)+ri[1]*L+ri[2]);
    }
    else
    {
        ncell = (int)((ri[0]-1)*pow(L,2)+ri[1]*L+ri[2]);
    }
    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_i = eLG_i + etmp;

    if (ri[1]==0)
    {
        ncell = (int)(ri[0]*pow(L,2)+(L-1)*L+ri[2]);
    }
    else
    {
        ncell = (int)(ri[0]*pow(L,2)+(ri[1]-1)*L+ri[2]);
    }

    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_i = eLG_i + etmp;

    if (ri[2]==0)
    {
        ncell = (int)(ri[0]*pow(L,2)+ri[1]*L+(L-1));
    }
    else
    {
        ncell = (int)(ri[0]*pow(L,2)+ri[1]*L+(ri[2]-1));
    }
    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_i = eLG_i + etmp;

    if (n[i]==(-1))
    {
        n[i] = 1;
    }
    else
    {
        n[i] = -1;
    }

    eLG_f = 0.0;

    eLG_f = model.get_energy(lambda,n,i,etot_f);

    if (ri[0]==0)
    {
        ncell = (int)((L-1)*pow(L,2)+ri[1]*L+ri[2]);
    }
    else
    {
        ncell = (int)((ri[0]-1)*pow(L,2)+ri[1]*L+ri[2]);
    }
    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_f = eLG_f + etmp;

    if (ri[1]==0)
    {
        ncell = (int)(ri[0]*pow(L,2)+(L-1)*L+ri[2]);
    }
    else
    {
        ncell = (int)(ri[0]*pow(L,2)+(ri[1]-1)*L+ri[2]);
    }
    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_f = eLG_f + etmp;

    if (ri[2]==0)
    {
        ncell = (int)(ri[0]*pow(L,2)+ri[1]*L+(L-1));
    }
    else
    {
        ncell = (int)(ri[0]*pow(L,2)+ri[1]*L+(ri[2]-1));
    }
    etmp = model.get_energy(lambda,n,ncell,etmp);
    eLG_f = eLG_f + etmp;

    eLG = eLG + (eLG_f-eLG_i);
    etot = eLG;

    de = eLG_f-eLG_i;

    // acceptance criteria

    //rn = rnvec[rand() % n_rn];
    rn = rn_flip[i];

    boltz = exp( -de/temp );
    // boltz = 1.0/(1.0+exp(de/temp));

    if (rn <= boltz)
    {
        n0[i] = n[i];

        e0 = etot;

        e0_LG = eLG;

        nchange = 1;
    }
    else
    {
        n[i] = n0[i];

        etot = e0;

        eLG = e0_LG;
    }

}

double ising::separation(double r1[3], double r2[3])
{
    double d[3];
    d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;

    double r;
    int dimR;

    for (dimR=0; dimR<dim; dimR++)
    {
        if ((r1[dimR]-r2[dimR])>(double)(L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] - L;
        }
        else if((r1[dimR]-r2[dimR])<(double)(-L/2))
        {
            d[dimR] = r1[dimR]-r2[dimR] + L;
        }
        else
        {
            d[dimR] = r1[dimR]-r2[dimR];
        }
    }

    r = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));

    return r;
}

void ising::initarrays()
{
    n = (int*) calloc (pow(L,dim), sizeof(int));
    n0 = (int*) calloc (pow(L,dim), sizeof(int));

    mem_test = true;
}

ising::~ising()
{
    if(mem_test=true)
    {
    delete [] n;
    delete [] n0;

    cout << "Deallocate ising memory" << endl;

    }
}
