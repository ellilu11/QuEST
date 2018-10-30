#include <iostream>
#include <vector>
#include <complex>
#include <unordered_set>
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

int node2D::isColleague(node2D* other){
    double dist = abs(center - other->center);
    return (this!=other && lvl == other->lvl && dist <= size*sqrt(2.0)+1.0E-03);
}

/*
int node2D::isNeighbour(node2D* other){
    double dist = abs(center - other->center);

    return (this!=other && other->isLeaf 
        && dist <= (1.0+pow(2.0,lvl-other->lvl))*size*sqrt(2.0)/2.0+1.0E-03);
}

// flag == 0 : find neighbours (adjacent nodes of any level)
// flag != 0 : find colleagues (adjacent nodes of same level)
std::vector<node2D*> node2D::findNeighboursSlow(int flag){
        
    // find the master node
    node2D* master = this; 
    while ( master->lvl != 0 ) master = master->parent;

    // cout << "New box!" << endl;

    // loop through all subnodes of master node, finding those which are neighbours
    std::vector<node2D*> nList;
    for ( int i=0; i<master->numNodes(); i++ )
        if ( flag == 0 ){
            if (isNeighbour(master->findNode(i)))
                nList.push_back(master->findNode(i));
        } else {
            if (isColleague(master->findNode(i)))
                nList.push_back(master->findNode(i));
        }
   
    // cout << nList.size() << endl;
 
    return nList;

}*/

// find the neighbour of larger or equal size in a given direction
node2D* node2D::findNeighbourGeq(int dir){

    if ( lvl == 0 ) return NULL;

    node2D* node;
    if ( dir == 1 ) { // south
        if (order == 2) return parent->child[0];
        if (order == 3) return parent->child[1];

        node = parent->findNeighbourGeq(dir);
        if (node == NULL) return NULL;
        if (node->isLeaf) return node;
 
        if (order == 0) return node->child[2];
        if (order == 1) return node->child[3];

    } else if ( dir == 6 ) { // north
        if (order == 0) return parent->child[2];
        if (order == 1) return parent->child[3];
        
        node = parent->findNeighbourGeq(dir);
        if (node == NULL) return NULL;
        if (node->isLeaf) return node;
 
        if (order == 2) return node->child[0];
        if (order == 3) return node->child[1];        

    } else if ( dir == 3 ) { // west
        if (order == 1) return parent->child[0];
        if (order == 3) return parent->child[2];
        
        node = parent->findNeighbourGeq(dir);
        if (node == NULL) return NULL;
        if (node->isLeaf) return node;

        if (order == 0) return node->child[1];
        if (order == 2) return node->child[3];        
    
    } else if ( dir == 4 ) { // east
        if (order == 0) return parent->child[1];
        if (order == 2) return parent->child[3];
        
        node = parent->findNeighbourGeq(dir);
        if (node == NULL) return NULL;
        if (node->isLeaf) return node;

        if (order == 1) return node->child[0];
        if (order == 3) return node->child[2];        
    
    } else if ( dir == 0 ) { // southwest
        if (order == 3) return parent->child[0];
        
        if (order == 0) {
            node = parent->findNeighbourGeq(dir);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[3];
        } else if (order == 1){
            node = parent->findNeighbourGeq(1);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[2];
        } else if (order == 2){
            node = parent->findNeighbourGeq(3);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[1];
        }
    } else if ( dir == 2 ) { // southeast
        if (order == 2) return parent->child[1];
        
        if (order == 1) {
            node = parent->findNeighbourGeq(dir);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[2];
        } else if (order == 0){
            node = parent->findNeighbourGeq(1);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[3];
        } else if (order == 3){
            node = parent->findNeighbourGeq(4);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[0];
        }
    } else if ( dir == 5 ) { // northwest
        if (order == 1) return parent->child[2];
        
        if (order == 2) {
            node = parent->findNeighbourGeq(dir);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[1];
        } else if (order == 0){
            node = parent->findNeighbourGeq(3);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[3];
        } else if (order == 3){
            node = parent->findNeighbourGeq(6);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[0];
        }
    } else if ( dir == 7 ) { // northeast
        if (order == 0) return parent->child[3];
        
        if (order == 3) {
            node = parent->findNeighbourGeq(dir);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[0];
        } else if (order == 1){
            node = parent->findNeighbourGeq(4);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[2];
        } else if (order == 2){
            node = parent->findNeighbourGeq(6);
            if (node == NULL) return NULL;
            if (node->isLeaf) return node;
            return node->child[1];
        }
    }
    return NULL;
}

// find neighbours of any size in a given direction
std::vector<node2D*> node2D::findNeighboursDir(int dir){
    std::vector<node2D*> nList;
    std::vector<node2D*> nList1, nList2;

    if ( isLeaf ) nList.push_back(this); return nList;
    
    if ( dir == 1 ){ // south
        nList1 = child[2]->findNeighboursDir(dir);
        nList2 = child[3]->findNeighboursDir(dir);
    } else if ( dir == 6 ){ // north
        nList1 = child[0]->findNeighboursDir(dir);
        nList2 = child[1]->findNeighboursDir(dir);
    } else if ( dir == 3 ){ // west
        nList1 = child[1]->findNeighboursDir(dir);
        nList2 = child[3]->findNeighboursDir(dir);
    } else if ( dir == 4 ){ // east
        nList1 = child[0]->findNeighboursDir(dir);
        nList2 = child[2]->findNeighboursDir(dir);
    } else if ( dir == 0 ) // southwest
        nList1 = child[3]->findNeighboursDir(dir);
      else if ( dir == 2 ) // southeast
        nList1 = child[2]->findNeighboursDir(dir);
      else if ( dir == 5 ) // northwest
        nList1 = child[1]->findNeighboursDir(dir);
      else if ( dir == 7 ) // northeast
        nList1 = child[0]->findNeighboursDir(dir);

    nList.reserve( nList1.size() + nList2.size() );
    nList.insert( nList.end(), nList1.begin(), nList1.end() );
    nList.insert( nList.end(), nList2.begin(), nList2.end() );
 
    return nList;
}

// find neighbours (adjacent nodes of any size) of a node
std::vector<node2D*> node2D::findNeighbours(){
    
    std::vector<node2D*> nList1, nList;
    
    for ( int i = 0; i < 8; i++ ){
        if ( findNeighbourGeq(i) != NULL ){
            nList1 = findNeighbourGeq(i)->findNeighboursDir(i);
            nList.reserve( nList.size() + nList1.size() );
            nList.insert( nList.end(), nList1.begin(), nList1.end() );
        }
    }
    
    // Need to account for duplicate nodes
    std::unordered_set<node2D*> nSet;
    nSet.insert( nList.begin(), nList.end() );
    nList.clear();
    nList.insert( nList.end(), nSet.begin(), nSet.end() );

    /*int dups = 0;
    for ( int i = 0; i < nList.size(); i++ )
        for ( int j = 0; j < nList.size(); j++ )
            if ( i != j && nList[i] == nList[j] ) dups++;

    cout << dups << " ";*/

    return nList;
}

// find colleagues (adjacent nodes of same size) of a node
std::vector<node2D*> node2D::findColleagues(){
    
    std::vector<node2D*> nList;
    
    for ( int i = 0; i < 8; i++ )
        if ( findNeighbourGeq(i) != NULL )
            if ( lvl == findNeighbourGeq(i)->lvl )
                nList.push_back(findNeighbourGeq(i));

    // Need to account for duplicate nodes
    std::unordered_set<node2D*> nSet;
    nSet.insert( nList.begin(), nList.end() );
    nList.clear();
    nList.insert( nList.end(), nSet.begin(), nSet.end() );

    return nList;
}

std::vector<std::complex<double> > 
    node2D::shiftLocalExp(std::complex<double> z0){
    
    int p = coeffLocalExp.size()-1;
    std::vector<std::complex<double> > B = coeffLocalExp;

    for ( int j=0; j<=p-1; j++ )
        for ( int k=p-j-1; k<=p-1; k++ )
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


/* Need to update this to match evalPotTrg 
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
}*/

// for a leaf node, evaluates the potential at all target points within
std::vector<std::complex<double> > node2D::evalPotTrg(
    std::vector<std::complex<double> > Z, std::vector<std::complex<double> > Ztrg, std::vector<double> Q, double mind){

    std::vector<std::complex<double> > phi(iztrg.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighbours();
    double d;
    int ii;

    // cout << nList.size() << " ";

    for ( int i=0; i<iztrg.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Ztrg[iztrg[i]]-center,l);
        
        // evaluate direct (near field) contributions, excluding interactions in node proper
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                d = abs(Ztrg[iztrg[i]] - Z[ii]);
                if (d > mind)
                phi[i] += Q[ii]/d;
            }
        }

        // finally, evaluate interactions in the node proper
        for ( int k=0; k<iz.size(); k++ ){
            d = abs(Ztrg[iztrg[i]] - Z[iz[k]]);
            // if (d > mind) 
            phi[i] += Q[iz[k]]/d;
        }
    }

    return phi;
}
