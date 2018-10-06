#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include "node.h"
#include "mathfunc.h"

using namespace std;

node2D::node2D(std::vector<std::complex<double> > Z, int maxparts,
    std::vector<int> iz, int lvl, int order, double size, std::complex<double> center){

    this->iz = iz;
    this->lvl = lvl;
    this->order = order;
    this->size = size;
    this->center = center;
    for ( int j=0; j<4; j++ ) child.push_back(NULL); // allocate memory for the 4 child nodes

    int nparts = iz.size();
    if ( nparts > maxparts ){ // more particles than allowed in this box; continue to subdivide

        isLeaf = 0;
        std::vector<std::vector<int> > izChild(4);

        for ( int i=0; i<iz.size(); i++ ){
            if ((real(Z[iz[i]]) <= real(center)) && (imag(Z[iz[i]]) <= imag(center))){
                izChild[0].push_back(iz[i]);
            } else if (real(Z[iz[i]]) > real(center) && imag(Z[iz[i]]) <= imag(center)){
                izChild[1].push_back(iz[i]);
            } else if (real(Z[iz[i]]) <= real(center) && imag(Z[iz[i]]) > imag(center)){
                izChild[2].push_back(iz[i]);
            } else if (real(Z[iz[i]]) > real(center) && imag(Z[iz[i]]) > imag(center)){
                izChild[3].push_back(iz[i]);
            } 
        }

        double centerXleft = real(center)-size/4;
        double centerXright = real(center)+size/4;
        double centerYdown = imag(center)-size/4;
        double centerYup = imag(center)+size/4;

        std::vector<std::complex<double> > centerChild(4);        
        centerChild[0] = centerXleft+I*centerYdown;
        centerChild[1] = centerXright+I*centerYdown;
        centerChild[2] = centerXleft+I*centerYup;
        centerChild[3] = centerXright+I*centerYup;

        for ( int j=0; j<4; j++ ){
            child[j] = new node2D( Z, maxparts, izChild[j], lvl+1, j, size/2, centerChild[j]);
            child[j]->parent = this;
        }
        // ~iz();
    } else // stop subdiving
        isLeaf = 1;   
}

// assigns inodes for a node and all its subnodes
void node2D::evalInode(int i){
    
    this->inode = i;
    if (isLeaf) return;
    
    child[0]->inode = i+1;
    for ( int j=1; j<4; j++ )
        child[j]->inode = child[j-1]->inode + child[j-1]->numNodes();
   
    for ( int j=0; j<4; j++ )
        child[j]->evalInode(child[j]->inode); 

}

// for any node, assigns the truncated multipole coefficients due to all charges within
void node2D::evalCoeffMpole(std::vector<std::complex<double> > Z, std::vector<double> Q, int p){

    std::vector<std::complex<double> > B(p+1);
 
    // pre-calculate binomial coefficients
    std::vector<std::vector<int> > binomCoeffs(p, std::vector<int>(p));
    for ( int l=1; l<=p; l++ )
        for ( int k=1; k<=l; k++ )
            binomCoeffs[l-1][k-1] = binomCoeff(l-1,k-1);

    if ( isLeaf ){ 
        std::vector<std::complex<double> > A(p+1);
        std::complex<double> z0 = center - parent->center;
        int nparts = iz.size();

        // compute expansion coefficients a_k about center z_0
        for ( int i=0; i<nparts; i++ ) {
            A[0] += Q[iz[i]];
            for ( int k=1; k<=p; k++ )
                A[k] -= Q[iz[i]]*pow(Z[iz[i]],k)/(double)k;
        }
        
        // compute shifted expansion coefficients b_k about origin
        B[0] = A[0];
        for ( int l=1; l<=p; l++ ){
            B[l] = -A[0]*pow(center,l)/(double)l;
            for ( int k=1; k<=l; k++ )
                B[l] += A[k]*pow(z0,l-k)*(double)binomCoeffs[l-1][k-1];
        }
    } else {
        for ( int j=0; j<4; j++ )
            child[j]->evalCoeffMpole(Z, Q, p);

        std::vector<std::vector<std::complex<double> > > A(4);
        std::vector<std::complex<double> > z0(4);
        for ( int j=0; j<4; j++ ){
            A[j] = child[j]->coeffMpole;
            z0[j] = child[j]->center - center;
        }

        for ( int j=0; j<3; j++ ){
            B[0] += A[j][0];
            for ( int l = 1; l<=p; l++ ){
                B[l] += -A[j][0]*pow(z0[j],l)/(double)l;
                for ( int k=1; k<=l; k++ )
                    B[l] += A[j][k]*pow(z0[j],l-k)*(double)binomCoeffs[l-1][k-1];
            }
        }
    }

    this->coeffMpole = B;
}

/*
void node2D::evalMasterList(){

    for ( int i=0; i < numNodes(); i++ ){
        coeffMpoleList[i] = findNode(i)->coeffMpole; 
        centerList[i] = findNode(i)->center;
    }
}*/

// Assigns a node's interaction list as indices to other nodes
void node2D::evalIList(){

    if (lvl > 2){
        // find parent's neighbours
        node2D* gparent = this->parent->parent;
        node2D* ggparent = gparent->parent;
        std::vector<node2D*> nList;
    
        for ( int j=0; j<4; j++ ){
            // find the 3 neighbours which are the parent's siblings
            if ( parent->order != j )
                nList.push_back(gparent->child[j]); 

            // find the 5 neighbours which are not the parent's siblings
            if ( gparent->order != j && !ggparent->child[j]->isLeaf ){
                for ( int k=0; k<4; k++ ) {
                    if (parent->isNeighbour(ggparent->child[j]->child[k]))
                        nList.push_back(ggparent->child[j]->child[k]);
                }
            }
        }
    
        // extract interaction list from parent's neighbours;
        for ( int j=0; j<nList.size(); j++){
            if ( !nList[j]->isLeaf ){
                for ( int k=0; k<4; k++){
                    if (!isNeighbour(nList[j]->child[k]))
                        ilist.push_back(nList[j]->child[k]);
                }
            }
        }
    } 

    if ( !isLeaf ) 
        for ( int j=0; j<4; j++ )
            child[j]->evalIList();
}

// Finds the local expansion coefficients due to all particles in a node's interaction list
void node2D::evalCoeffLocalExp( int p ){
    
    std::vector<std::vector<int> > binomCoeffs(p, std::vector<int>(p));
    for ( int l=1; l<=(2*p-1); l++ )
        for ( int k=1; k<=l; k++ )
            binomCoeffs[l-1][k-1] = binomCoeff(l-1,k-1);
  
    std::vector<std::vector<std::complex<double> > > Bs;

    std::vector<std::complex<double> > B(p+1);
    std::complex<double> z0;
    for ( int i=0; i<ilist.size(); i++ ){
        z0 = ilist[i]->center - center;
        Bs[i] = ilist[i]->coeffMpole;

        B[0] = Bs[i][0]*log(-z0);
        for ( int k=1; k<=p; k++ )
            B[0] += Bs[i][k]*pow(-1,k)/pow(z0,k);

        for ( int l=1; l<=p; l++ ){
            B[l] = -Bs[i][0]/((double)l*pow(z0,l));
            for ( int k=1; k<=p; k++ )
                B[l] += Bs[i][k]*pow(-1,k)/pow(z0,l+k)*(double)binomCoeffs[l+k-1][k-1];
        }
    } 
    this->coeffLocalExp = B;

    if ( !isLeaf )
        for ( int j=0; j<4; j++ )
            child[j]->evalCoeffLocalExp(p);
}

/*
void node2D::evalCoeffLocalExpSum(int p){
    
    std::vector<std::complex<double> > A = this->coeffLocalExp;
    std::vector<std::vector<std::complex<double> > > B(4); 
    std::vector<std::complex<double> > z0(4);
    z0[0] = child[0]->center - center;
    z0[1] = child[1]->center - center;
    z0[2] = child[2]->center - center;
    z0[3] = child[3]->center - center;

    for ( int j=0; j<4; j++) {
        B[j] = shiftLocalExp(A, z0[j]);
        for ( int k=0; k<p; k++)
            A[k] += B[j][k];
    }
    this->coeffLocalExp = A;

    is ( !isLeaf ){
        child[0]->evalCoeffLocalExpSum(p);
        child[1]->evalCoeffLocalExpSum(p);
        child[2]->evalCoeffLocalExpSum(p);
        child[3]->evalCoeffLocalExpSum(p);
    }
}*/

