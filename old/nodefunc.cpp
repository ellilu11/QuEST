#include <iostream>
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

int node2D::isColleague(node2D* other){
    double dist = abs(center - other->center);
    return (this!=other && lvl == other->lvl && dist <= size*sqrt(2.0)+1.0E-03);
}

int node2D::isNeighbour(node2D* other){
    double dist = abs(center - other->center);
    // cout << lvl << " " << other->lvl << " " << isLeaf << " " << other->isLeaf << " " << dist << " " << ((1.0+pow(2.0,lvl-other->lvl))*size*sqrt(2.0)/2.0+1.0E-03)
    //     << " " << (this!=other && other->isLeaf && dist <= (1.0+pow(2.0,lvl-other->lvl))*size*sqrt(2.0)/2.0+1.0E-03) << endl;

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

}

/* THESE DO NOT WORK!
// returns array of pointers to colleagues (adjacent nodes of same level)
std::vector<node2D*> node2D::findColleagues(){

    std::vector<node2D*> nList;
    node2D* gparent;
    node2D* ggparent;

    if ( lvl > 1 ) gparent = this->parent->parent;
    if ( lvl > 2 ) ggparent = this->parent->parent->parent;

    for ( int j=0; j<4; j++ ){

        // find the 3 sibling colleagues
        if ( lvl > 0 && order != j )
            nList.push_back(parent->child[j]);

        // find the 5 non-sibling colleagues
        if ( lvl > 1 && parent->order != j ) { 
            if (!gparent->child[j]->isLeaf ) {
                for ( int k=0; k<4; k++ ) {
                    if (isColleague(gparent->child[j]->child[k]))
                        nList.push_back(gparent->child[j]->child[k]);
                }
            }
        }
        // if the node is not an "interior" node of its gparent, have to look at ggparent
        if ( lvl > 2 && gparent->order != j ) { 
        // if ( lvl > 2 && abs(center - gparent->center) > size*sqrt(2.0)/2.0 && gparent->order != j ){
            if ( !ggparent->child[j]->isLeaf ) {
                for ( int k=0; k<4; k++ ){
                    if ( !ggparent->child[j]->child[k]->isLeaf ){
                        for ( int l=0; l<4; l++ ){
                            if (isColleague(ggparent->child[j]->child[k]->child[l]))
                                nList.push_back(ggparent->child[j]->child[k]->child[l]);
                        }
                    }
                } 
            }
        }
    }         

    return nList;
}

// returns array of pointers to neighbours (adjacent nodes of any level)
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
        if ( lvl > 2 && gparent->order != j ) { 
        // if ( lvl > 2 && abs(center - gparent->center) > size*sqrt(2.0)/2.0 && gparent->order != j ){
            if ( !ggparent->child[j]->isLeaf ) {
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
}*/

std::vector<std::complex<double> > 
    node2D::shiftLocalExp(std::complex<double> z0){
    
    int p = coeffLocalExp.size()-1;
    std::vector<std::complex<double> > B = coeffLocalExp;

    for ( int j=0; j<=p-1; j++ )
        for ( int k=p-j-1; k<=p-1; k++ )
            B[k] = B[k] - z0*B[k+1];

    return B;
}

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
}

// for a leaf node, evaluates the potential at all source points within
std::vector<std::complex<double> > node2D::evalPotSrc(std::vector<std::complex<double> > Z,  std::vector<double> Q, double mind){

    std::vector<std::complex<double> > phi(iz.size());
    int p = coeffLocalExp.size()-1;
    std::vector<node2D*> nList = findNeighboursSlow(0);
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
    std::vector<node2D*> nList = findNeighboursSlow(0);
    double d;
    int ii;

    cout << nList.size() << " ";

    for ( int i=0; i<iztrg.size(); i++ ) {

        // evaluate multipole (far field) contributions
        for ( int l=0; l<p; l++ )
            phi[i] += coeffLocalExp[l]*pow(Ztrg[iztrg[i]]-center,l);
        
        // evaluate direct (near field) contributions, excluding interactions in node proper
        for ( int j=0; j<nList.size(); j++ ){
            for ( int k=0; k<nList[j]->iz.size(); k++ ){
                ii = nList[j]->iz[k];
                d = abs(Ztrg[iztrg[i]] - Z[ii]);
                // if (d > mind)
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
