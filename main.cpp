#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <omp.h>

#include "node.h"
#include "mathfunc.h"
#include "fstream.h"

using namespace std;

std::ofstream nsrcFile, ntrgFile, ZFile, ZtrgFile, mpoleFile;
std::ofstream potFile, fldFile;

int main(int argc, char *argv[]) {

    // Initialize parameters; eventually these will be read from file
    int Nparts = 100000; // Number of particles
    int Ntargs = 10000; // Number of targets (field points)
    int Maxparts = 1000; // maximum number of particles allowed in any box
    int Iprec = 0;  // error tolerance flag
    int srcDist = 0; // 0 - uniform; 1 - Gaussian; etc...
    int trgDist = 1; // 0 - uniform 2-D; 1 - uniform along x-axis; 2 - uniform along y-axis 
    int potFlag = 2; // 0 - don't eval pot; 1 - eval pot at srcs; 2 - eval pot at trgs
    int fldFlag = 0; 

    double epsFmm; int p; // error tolerance parameter & # of terms in truncation
    switch (Iprec) {
        case -2 : epsFmm = 0.5E-0; break;
        case -1 : epsFmm = 0.5E-1; break;
        case 0 : epsFmm = 0.5E-2; break;
        case 1 : epsFmm = 0.5E-3; break;
        case 2 : epsFmm = 0.5E-6; break;
        case 3 : epsFmm = 0.5E-9; break;
        case 4 : epsFmm = 0.5E-12; break;
        case 5 : epsFmm = 0.5E-15; break;   
    }
    p = ceil(log2(1/epsFmm));

    double startTime = omp_get_wtime();
    double currTime;
    double u1, u2, x, y;
    
    // Initialize source distribution
    std::vector<std::complex<double> > sources(Nparts);
    std::vector<double> charges(Nparts);

    for ( int i=0; i<Nparts; i++ ) {
        u1 = rand()/(double)RAND_MAX;
        u2 = rand()/(double)RAND_MAX;
        if (srcDist == 0) {
            x = u1*10-5;
            y = u2*10-5;
        } else if (srcDist == 1) {
            x = sqrt(-2*log(u2))*cos(2*M_PI*u1);
            y = sqrt(-2*log(u2))*sin(2*M_PI*u1);
        }
        sources[i] = x+I*y;
        charges[i] = 1.0;
    }

    // Initialize the master box
    std::vector<int> isources(Nparts); // index of all particles in sources array
    double xmin = real(sources[0]);
    double ymin = imag(sources[0]);
    double xmax = xmin;
    double ymax = ymax;

    for ( int i=0; i<Nparts; i++ ){
        isources[i] = i;

        if (real(sources[i]) < xmin) xmin = real(sources[i]);
        if (imag(sources[i]) < ymin) ymin = imag(sources[i]);
        if (real(sources[i]) > xmax) xmax = real(sources[i]);
        if (imag(sources[i]) > ymax) ymax = imag(sources[i]);
    }

    double size = xmax - xmin;
    double sizey = ymax - ymin;
    if (sizey > size) size = sizey;
    double centerX = (xmin+xmax)/2;
    double centerY = (ymin+ymax)/2;
    std::complex<double> center = centerX+I*centerY;
    double mind = size/(double)sqrt(Nparts);
    // cout << mind << endl;

    // Initialize target distribution
    std::vector<std::complex<double> > targets(Ntargs);
    std::vector<int> itargets(Ntargs);
    double delta = size/(double)Ntargs;

    for ( int i=0; i<Ntargs; i++ ) {
        x = trgDist == 1 ? xmin + i*delta : 0.0;
        y = trgDist == 2 ? ymin + i*delta : 1.0;
        targets[i] = x+I*y;
        itargets[i] = i;
    }
  
    cout << "Initialized domain... " << omp_get_wtime() - startTime << "s" << endl;
    currTime = omp_get_wtime();  

    node2D* treefmm2D = new node2D( sources, targets, Maxparts, isources, itargets, 0, 0, size, center );
    // treefmm2D->evalInode(0);
    cout << "Created FMM2D tree structure... " << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();
 
    treefmm2D->evalCoeffMpole(sources, charges, p);
    cout << "Evaluated multipole expansions... " << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();
 
    treefmm2D->evalIList();
    cout << "Determined interaction lists..." << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();

    treefmm2D->evalCoeffLocalExp(p);
    cout << "Evaluated local expansions... " << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();
 
    treefmm2D->evalCoeffLocalExpSum(p);
    cout << "Summed local expansions... " << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();
 
    potFile.open("out/pot.bin");
    treefmm2D->fprintPot( sources, targets, charges, potFlag, mind );
    potFile.close();
    cout << "Evaluated potentials... " << omp_get_wtime() - currTime << "s" << endl;
    currTime = omp_get_wtime();
 
    nsrcFile.open("out/nsrc.csv");    
    ZFile.open("out/Z.bin", std::ios::binary);
    treefmm2D->fprintList(sources);
    nsrcFile.close();
    ZFile.close();

    /*ntrgFile.open("out/ntrg.csv");    
    ZtrgFile.open("out/Ztrg.bin", std::ios::binary);
    treefmm2D->fprintZtrg(targets);
    ntrgFile.close();
    ZtrgFile.close();

    mpoleFile.open("out/mpole.csv");
    treefmm2D->fprintMpole(); 
    mpoleFile.close();*/

    cout << "p: " << p << endl;
    cout << "Number of boxes: " << treefmm2D->numNodes() << endl;
    cout << "Time elapsed: " << currTime - startTime << "s" << endl;

    delete treefmm2D;

    return 0;
}
