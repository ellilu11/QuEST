#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#include "node.h"
#include "fstream.h"

using namespace std;

// counts number of subnodes of a node (including self)
int node2D::numNodes(){
    if ( isLeaf ) return 1;
    else return 1 + child[0]->numNodes() + child[1]->numNodes() +
           child[2]->numNodes() + child[3]->numNodes();
}

int node2D::isNeighbour(node2D* other){
    double dist = abs(center - other->center);
    return (lvl == other->lvl && dist <= size*sqrt(2));
}

// returns array of pointers to neighbour nodes
std::vector<node2D*> node2D::findNeighbours(){

    std::vector<node2D*> nList;
    node2D* gparent = this->parent->parent;

    for ( int j=0; j<4; j++ ){
        // find the 3 sibling neighbours
        if ( order != j )
            nList.push_back(parent->child[j]);

        // find the 5 non-sibling neighbours
        if ( parent->order != j && !gparent->child[j]->isLeaf ){
            for ( int k=0; k<4; k++ ) {
                if (isNeighbour(gparent->child[j]->child[k]))
                    nList.push_back(gparent->child[j]->child[k]);
            }
        }
    }

    return nList;
}

std::vector<std::complex<double> > 
    node2D::shiftLocalExp(std::complex<double> z0){
    
    int p = coeffLocalExp.size()-1;
    std::vector<std::complex<double> > B = coeffLocalExp;

    for ( int j=0; j<p; j++ )
        for ( int k=p-j-1; k<p; k++ )
            B[k] = B[k] - z0*B[k+1];

    return B;
}

/*
// returns pointer to ith subnode of a given node
node2D* node2D::findNode(int i){
    if (this->inode == i)
        return this;
    else if (i < child[1]->inode)
        return child[0]->findNode(i);
    else if (i < child[2]->inode)
        return child[1]->findNode(i);
    else if (i < child[3]->inode)
        return child[2]->findNode(i);
    else
        return child[3]->findNode(i);
}*/

// for a leaf node, evaluates the potential at all source points within
std::vector<std::complex<double> > node2D::evalPotSrc(std::vector<std::complex<double> > Z,  std::vector<double> Q){

    std::vector<std::complex<double> > phi(iz.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighbours();
    double r;
    int ii;

    for ( int i=0; i<iz.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Z[iz[i]]-center,l);

        // evaluate direct (near field) contributions
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                r = abs(Z[iz[i]] - Z[ii]);
                phi[i] += Q[ii]/r;
            }
        }
    }

    return phi;
}

// for a leaf node, evaluates the potential at all target points within
std::vector<std::complex<double> > node2D::evalPotTrg(
    std::vector<std::complex<double> > Z, std::vector<std::complex<double> > Ztrg, std::vector<double> Q){

    std::vector<std::complex<double> > phi(iztrg.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighbours();
    double r;
    int ii;

    for ( int i=0; i<iztrg.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Ztrg[iztrg[i]]-center,l);

        // evaluate direct (near field) contributions
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                r = abs(Ztrg[iztrg[i]] - Z[ii]);
                phi[i] += Q[ii]/r;
            }
        }
    }

    return phi;
}

void node2D::fprintZ(std::vector<std::complex<double> > Z){
    if (this != NULL) {

        ::nsrcFile << iz.size() << endl;
        for ( int i=0; i<iz.size(); i++ ){
            ::ZFile.write((char*)&real(Z[iz[i]]), sizeof(double));
            ::ZFile.write((char*)&imag(Z[iz[i]]), sizeof(double));
            // ::ZFile << real(Z[iz[i]]) << "," << imag(Z[iz[i]]) << endl;
        }

        for ( int j=0; j<4; j++ ) child[j]->fprintZ(Z);
    }
}

void node2D::fprintZtrg(std::vector<std::complex<double> > Ztrg){
    if (this != NULL) {

        ::ntrgFile << iztrg.size() << endl;
        for ( int i=0; i<iztrg.size(); i++ ){
            ::ZtrgFile.write((char*)&real(Ztrg[iztrg[i]]), sizeof(double));
            ::ZtrgFile.write((char*)&imag(Ztrg[iztrg[i]]), sizeof(double));
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
    std::vector<double> Q, int flag ){
    if ( isLeaf ) {
        if ( flag == 1 ) {
            std::vector<std::complex<double> > pot = evalPotSrc( Z, Q ); 
            for ( int i=0; i<iz.size(); i++ ){
                ::potFile.write((char*)&real(Z[iz[i]]), sizeof(double));
                ::potFile.write((char*)&imag(Z[iz[i]]), sizeof(double));
                ::potFile.write((char*)&real(pot[i]), sizeof(double));
                // ::potFile << Re(Z[iz[i]]) << "," << Im(Z << "," << pot[i] << endl;
            } 
        } else if ( flag == 2 ) {
            std::vector<std::complex<double> > pot = evalPotTrg( Z, Ztrg, Q ); 
            for ( int i=0; i<iztrg.size(); i++ ){
                ::potFile.write((char*)&real(Ztrg[iztrg[i]]), sizeof(double));
                ::potFile.write((char*)&imag(Ztrg[iztrg[i]]), sizeof(double));
                ::potFile.write((char*)&real(pot[i]), sizeof(double));
                // ::potFile << Ztrg[iztrg[i]] << "," << pot[i] << endl;
            } 
        }
    } else 
        for ( int j=0; j<4; j++ )
            child[j]->fprintPot( Z, Ztrg, Q, flag );
    

}
