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

int main(int argc, char *argv[]) {

    // Initialize parameters; eventually these will be read from file
    int Nparts = 100000; // total number of particles
    int Ntargs = 100; // Number of targets (field points)
    int Maxparts = 1000; // maximum number of particles allowed in any box
    int Iprec = 0;  // error tolerance flag
    int distFlag = 1;   // 0 - uniform, 1 - Gaussian
    
    time_t startTime, endTime;
    time(&startTime);

    // Initialize particle distribution
    std::vector<std::complex<double> > sources(Nparts);
    std::vector<double> charges(Nparts);

    double u1, u2, xsources, ysources;
    for (int i=0;i<Nparts;i++){
        u1 = rand()/(double)RAND_MAX;
        u2 = rand()/(double)RAND_MAX;
        if (distFlag == 0) {
            xsources = u1*10-5;
            ysources = u2*10-5;
        } else if (distFlag == 1) {
            xsources = sqrt(-2*log(u2))*cos(2*M_PI*u1);
            ysources = sqrt(-2*log(u2))*sin(2*M_PI*u1);
        }
        sources[i] = xsources+I*ysources;
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

    node2D* treefmm2D = new node2D(sources, Maxparts, isources, 0, 1, size, center ); // construct tree structure
    treefmm2D->evalInode(0); // assign inodes;
    treefmm2D->evalCoeffMpole(sources, charges, p); // calculate and merge multipole coefficients
    treefmm2D->evalMasterList();
    treefmm2D->evalIList(treefmm2D);
    // treefmm2D->evalCoeffLocalExp(treefmm2D, p);

    // write results to file
    npartsFile.open("out/nparts.csv");    
    ZFile.open("out/Z.bin", std::ios::binary);
    treefmm2D->fprintZ(sources);
    npartsFile.close();
    ZFile.close();

    std::ofstream mpoleFile; 
    mpoleFile.open("out/mpole.csv");
    for ( int i=0; i<treefmm2D->coeffMpoleList.size(); i++ ){
        for ( int j=0; j<p+1; j++ )
            mpoleFile << treefmm2D->coeffMpoleList[i][j] << ",";
        mpoleFile << endl;
    }    
    mpoleFile.close();

    time(&endTime);
    cout << "p: " << p << endl;
    cout << "Number of boxes: " << treefmm2D->numNodes() << endl;
    cout << "Time elapsed: " << endTime - startTime << "s" << endl;

    delete treefmm2D;

    return 0;
}
