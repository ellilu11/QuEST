#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include "node.h"
#include "fstream.h"

using namespace std;

void node2D::fprintZ(std::vector<std::complex<double> > Z){
    if (this != NULL) {

        ::nsrcFile << iz.size() << endl;
        
        double Zreal, Zimag;
        for ( int i=0; i<iz.size(); i++ ){
            Zreal = real(Z[iz[i]]);
            Zimag = imag(Z[iz[i]]);    
            ::ZFile.write((char*)&Zreal, sizeof(double));
            ::ZFile.write((char*)&Zimag, sizeof(double));
            // ::ZFile << real(Z[iz[i]]) << "," << imag(Z[iz[i]]) << endl;
        }

        for ( int j=0; j<4; j++ ) child[j]->fprintZ(Z);
    }
}

void node2D::fprintZtrg(std::vector<std::complex<double> > Ztrg){
    if (this != NULL) {

        ::ntrgFile << iztrg.size() << endl;

        double Zreal, Zimag;
        for ( int i=0; i<iztrg.size(); i++ ){
            Zreal = real(Ztrg[iztrg[i]]);
            Zimag = imag(Ztrg[iztrg[i]]);
            ::ZtrgFile.write((char*)&Zreal, sizeof(double));
            ::ZtrgFile.write((char*)&Zimag, sizeof(double));
        }

        for ( int j=0; j<4; j++ ) child[j]->fprintZtrg(Ztrg);
    }
}

void node2D::fprintMpole(){
    if (this != NULL) {

        for ( int i=0; i<coeffMpole.size(); i++ ){
            ::mpoleFile << coeffMpole[i];
        }
        ::mpoleFile << endl;

        for ( int j=0; j<4; j++ ) child[j]->fprintMpole();
    }
}

void node2D::fprintPot( std::vector<std::complex<double> > Z, std::vector<std::complex<double> > Ztrg, 
    std::vector<double> Q, int flag, double mind ){
    if ( isLeaf ) {
        double Zreal, Zimag, potreal, potimag, potabs;
        if ( flag == 1 ) {
            std::vector<std::complex<double> > pot = evalPotSrc( Z, Q, mind ); 
            for ( int i=0; i<iz.size(); i++ ){
                // potAbs = abs(pot[i]);
                Zreal = real(Z[iz[i]]);
                Zimag = imag(Z[iz[i]]);
                potreal = real(pot[i]);
                potimag = imag(pot[i]);
                ::potFile.write((char*)&Zreal, sizeof(double));
                ::potFile.write((char*)&Zimag, sizeof(double));
                ::potFile.write((char*)&potreal, sizeof(double));
                ::potFile.write((char*)&potimag, sizeof(double));
                // ::potFile << Re(Z[iz[i]]) << "," << Im(Z << "," << pot[i] << endl;
            } 
        } else if ( flag == 2 ) {
            std::vector<std::complex<double> > pot = evalPotTrg( Z, Ztrg, Q, mind );
            for ( int i=0; i<iztrg.size(); i++ ){
                Zreal = real(Ztrg[iztrg[i]]);
                Zimag = imag(Ztrg[iztrg[i]]);
                potreal = real(pot[i]);
                potimag = imag(pot[i]);
                potabs = abs(pot[i]);
                ::potFile.write((char*)&Zreal, sizeof(double));
                ::potFile.write((char*)&Zimag, sizeof(double));
                ::potFile.write((char*)&potreal, sizeof(double));
                ::potFile.write((char*)&potimag, sizeof(double));
                ::potFile.write((char*)&potabs, sizeof(double));
                // ::potFile << Ztrg[iztrg[i]] << "," << pot[i] << endl;
            } 
        }
    } else 
        for ( int j=0; j<4; j++ )
            child[j]->fprintPot( Z, Ztrg, Q, flag, mind );
    

}

void node2D::fprintList( std::vector<std::complex<double> > Z ){

    int s = 0;
    std::vector<int> iiz;

    if ( lvl > 2 ) {
    std::vector<node2D*> list = findNeighboursSlow(0);

    double Zreal, Zimag;
    for ( int i=0; i<list.size(); i++ ){
        iiz = list[i]->iz;
        s += iiz.size();
        // cout << i << "," << s << endl;
        for ( int j=0; j<iiz.size(); j++ ){
            Zreal = real(Z[iiz[j]]);
            Zimag = imag(Z[iiz[j]]);
            ::ZFile.write((char*)&Zreal, sizeof(double));
            ::ZFile.write((char*)&Zimag, sizeof(double));
        }
    }
    ::nsrcFile << s << endl;
    }

    if ( !isLeaf ) for ( int j=0; j<4; j++ ) child[j]->fprintList(Z);
    
}
