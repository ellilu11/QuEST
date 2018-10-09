#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#include "node.h"
#include "fstream.h"
#include "mathfunc.h"

using namespace std;

// counts number of subnodes of a node (including self)
int node2D::numNodes(){
    if ( isLeaf ) return 1;
    else return 1 + child[0]->numNodes() + child[1]->numNodes() +
           child[2]->numNodes() + child[3]->numNodes();
}

int node2D::isNeighbour(node2D* other){
    double dist = abs(center - other->center);
    return (dist <= (1+pow(2,other->lvl-lvl))*size*sqrt(2.0)/2.0);
    // return (lvl == other->lvl && dist <= size*sqrt(2.0));
}

// returns array of pointers to neighbour nodes
std::vector<node2D*> node2D::findNeighbours(){

    std::vector<node2D*> nList;
    node2D* gparent;
    node2D* ggparent;

    if ( lvl > 1 ) gparent = this->parent->parent;
    if ( lvl > 2 ) ggparent = this->parent->parent->parent;

    for ( int j=0; j<4; j++ ){

        // find the 3 sibling neighbours
        if ( lvl > 0 && order != j )
            nList.push_back(parent->child[j]);

        // find the 5 non-sibling neighbours
        if ( lvl > 1 && parent->order != j ) { 
            if (!gparent->child[j]->isLeaf ) {
            for ( int k=0; k<4; k++ ) {
                if (isNeighbour(gparent->child[j]->child[k]))
                    nList.push_back(gparent->child[j]->child[k]);
            }
            } else {
                if (isNeighbour(gparent->child[j]))
                    nList.push_back(gparent->child[j]);
            }
        }
        // if the node is not an "interior" node of its gparent, have to look at ggparent
        if ( lvl > 2 && abs(center - gparent->center) > size*sqrt(2.0)/2.0 && gparent->order != j ){
           if ( !ggparent->child[j]->isLeaf ){
           for ( int k=0; k<4; k++ ){
               if ( !ggparent->child[j]->child[k]->isLeaf ){
                    for ( int l=0; l<4; l++ ){
                        if (isNeighbour(ggparent->child[j]->child[k]->child[l]))
                            nList.push_back(ggparent->child[j]->child[k]->child[l]);
                    }
               } else {
                    if (isNeighbour(ggparent->child[j]->child[k]))
                            nList.push_back(ggparent->child[j]->child[k]);
               }
           } 
           } else {
                if (isNeighbour(ggparent->child[j]))
                    nList.push_back(ggparent->child[j]);
           }
        }
    }         

    return nList;
}

std::vector<std::complex<double> > 
    node2D::shiftLocalExp(std::complex<double> z0){
    
    int p = coeffLocalExp.size()-1;
    std::vector<std::complex<double> > B = coeffLocalExp;
    // for ( int i=0; i<=p; i++ ) B[i] = coeffLocalExp[i];

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
std::vector<std::complex<double> > node2D::evalPotSrc(std::vector<std::complex<double> > Z,  std::vector<double> Q, double mind){

    std::vector<std::complex<double> > phi(iz.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighbours();
    double d;
    int ii;

    for ( int i=0; i<iz.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Z[iz[i]]-center,l);

        // evaluate direct (near field) contributions
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                d = abs(Z[iz[i]] - Z[ii]);
                if (d > mind) phi[i] += Q[ii]/d;
            }
        }
    }

    return phi;
}

// for a leaf node, evaluates the potential at all target points within
std::vector<std::complex<double> > node2D::evalPotTrg(
    std::vector<std::complex<double> > Z, std::vector<std::complex<double> > Ztrg, std::vector<double> Q, double mind){

    std::vector<std::complex<double> > phi(iztrg.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighbours();
    double d;
    int ii;

    for ( int i=0; i<iztrg.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Ztrg[iztrg[i]]-center,l);

        // evaluate direct (near field) contributions
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                d = abs(Ztrg[iztrg[i]] - Z[ii]);
                if (d > mind) phi[i] += Q[ii]/d;
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
    std::vector<double> Q, int flag, double mind ){
    if ( isLeaf ) {
        double potAbs;
        if ( flag == 1 ) {
            std::vector<std::complex<double> > pot = evalPotSrc( Z, Q, mind ); 
            for ( int i=0; i<iz.size(); i++ ){
                // potAbs = abs(pot[i]); 
                ::potFile.write((char*)&real(Z[iz[i]]), sizeof(double));
                ::potFile.write((char*)&imag(Z[iz[i]]), sizeof(double));
                ::potFile.write((char*)&real(pot[i]), sizeof(double));
                ::potFile.write((char*)&imag(pot[i]), sizeof(double));
                // ::potFile << Re(Z[iz[i]]) << "," << Im(Z << "," << pot[i] << endl;
            } 
        } else if ( flag == 2 ) {
            std::vector<std::complex<double> > pot = evalPotTrg( Z, Ztrg, Q, mind );
            for ( int i=0; i<iztrg.size(); i++ ){
                potAbs = abs(pot[i]);
                ::potFile.write((char*)&real(Ztrg[iztrg[i]]), sizeof(double));
                ::potFile.write((char*)&imag(Ztrg[iztrg[i]]), sizeof(double));
                ::potFile.write((char*)&real(pot[i]), sizeof(double));
                ::potFile.write((char*)&imag(pot[i]), sizeof(double));
                ::potFile.write((char*)&potAbs, sizeof(double));
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
    std::vector<node2D*> list = findNeighbours();

    for ( int i=0; i<list.size(); i++ ){
        iiz = list[i]->iz;
        s += iiz.size();
        // cout << i << "," << s << endl;
        for ( int j=0; j<iiz.size(); j++ ){
            ::ZFile.write((char*)&real(Z[iiz[j]]), sizeof(double));
            ::ZFile.write((char*)&imag(Z[iiz[j]]), sizeof(double));
        }
    }
    ::nsrcFile << s << endl;
    }

    if ( !isLeaf ) for ( int j=0; j<4; j++ ) child[j]->fprintList(Z);
    
}
