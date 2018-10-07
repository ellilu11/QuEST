#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <cmath>
#include <ctime>

#include "node.h"
#include "mathfunc.h"
#include "fstream.h"

using namespace std;

std::ofstream npartsFile;
std::ofstream ZFile;
std::ofstream mpoleFile;
std::ofstream potFile;
std::ofstream fldFile;

int main(int argc, char *argv[]) {

    // Initialize parameters; eventually these will be read from file
    int Nparts = 100000; // Number of particles
    int Ntargs = 100; // Number of targets (field points)
    int Maxparts = 1000; // maximum number of particles allowed in any box
    int Iprec = 0;  // error tolerance flag
    int srcDist = 1; // 0 - uniform; 1 - Gaussian
    int trgDist = 1; // 0 - uniform 2-D
                     // 1 - uniform along x-axis
                     // 2 - uniform along y-axis 
    int evalPotSrc = 0; 
    int evalPotTrg = 1;
    int evalFldSrc = 0;
    int evalFldTrg = 1;  

    time_t startTime, endTime;
    time(&startTime);

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

    // Initialize target distribution
    std::vector<std::complex<double> > targets(Ntargs);

    if ( evalPotTrg || evalFldTrg ) {
        if ( trgDist == 1 || trgDist == 2 ) {
            double delta = size/(double)Ntargs;

            for ( int i=0; i<Ntargs; i++ ) {
                x = xmin + i*delta;
                y = 0;
                targets[i] = x+I*y;    
            }
        }  
    }

    node2D* treefmm2D = new node2D(sources, targets, Maxparts, isources, 0, 1, size, center );
    // treefmm2D->evalInode(0); // assign inodes;
    treefmm2D->evalCoeffMpole(sources, charges, p); // calculate and merge multipole coefficients
    treefmm2D->evalIList();
    treefmm2D->evalCoeffLocalExp(p);
    treefmm2D->evalCoeffLocalExpSum(p);
    
    // write results to file
    npartsFile.open("out/nparts.csv");    
    ZFile.open("out/Z.bin", std::ios::binary);
    treefmm2D->fprintZ(sources);
    npartsFile.close();
    ZFile.close();

    mpoleFile.open("out/mpole.csv");
    treefmm2D->fprintMpole(); 
    mpoleFile.close();

    time(&endTime);
    cout << "p: " << p << endl;
    cout << "Number of boxes: " << treefmm2D->numNodes() << endl;
    cout << "Time elapsed: " << endTime - startTime << "s" << endl;

    delete treefmm2D;

    return 0;
}
