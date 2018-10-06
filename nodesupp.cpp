#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

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

void node2D::fprintZ(std::vector<std::complex<double> > Z){
    if (this != NULL) {
        // cout << "NEW NODE w/ " << iz.size() << " particles centered at " << center << endl;

        ::npartsFile << iz.size() << endl;
        for ( int i=0; i<iz.size(); i++ ){
            ::ZFile.write((char*)&real(Z[iz[i]]), sizeof(double));
            ::ZFile.write((char*)&imag(Z[iz[i]]), sizeof(double));
            // ::ZFile << real(Z[iz[i]]) << "," << imag(Z[iz[i]]) << endl;
        }

        for ( int j=0; j<4; j++ ) child[j]->fprintZ(Z);
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

