#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include "node.h"
#include "mathfunc.h"

using namespace std;

node2D::node2D(std::vector<std::complex<double> > Z, std::vector<std::complex<double> > Ztrg, int maxparts,
    std::vector<int> iz, std::vector<int> iztrg, int lvl, int order, double size, std::complex<double> center){

    this->iz = iz;
    this->iztrg = iztrg;
    this->lvl = lvl;
    this->order = order;
    this->size = size;
    this->center = center;
    child.reserve(4);
    // for ( int j=0; j<4; j++ ) child.push_back(NULL); // allocate memory for the 4 child nodes

    if ( iz.size() > maxparts ){ // more particles than allowed in this box; continue to subdivide

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

        // cout << iztrg.size() << endl;
        std::vector<std::vector<int> > iztrgChild(4);
        for ( int i=0; i<iztrg.size(); i++ ){
            if ((real(Ztrg[iztrg[i]]) <= real(center)) && (imag(Ztrg[iztrg[i]]) <= imag(center))){
                iztrgChild[0].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) > real(center) && imag(Ztrg[iztrg[i]]) <= imag(center)){
                iztrgChild[1].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) <= real(center) && imag(Ztrg[iztrg[i]]) > imag(center)){
                iztrgChild[2].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) > real(center) && imag(Ztrg[iztrg[i]]) > imag(center)){
                iztrgChild[3].push_back(iztrg[i]);
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
            child[j] = new node2D( Z, Ztrg, maxparts, izChild[j], iztrgChild[j], lvl+1, j, size/2, centerChild[j]);
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
                A[k] -= Q[iz[i]]*pow(Z[iz[i]]-parent->center,k)/(double)k;
                // A[k] -= Q[iz[i]]*pow(Z[iz[i]],k)/(double)k;
        }
        
        // compute shifted expansion coefficients b_k about origin
        B[0] = A[0];
        for ( int l=1; l<=p; l++ ){
            B[l] = -A[0]*pow(z0,l)/(double)l;
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

        for ( int j=0; j<4; j++ ){
            B[0] += A[j][0];
            for ( int l=1; l<=p; l++ ){
                B[l] += -A[j][0]*pow(z0[j],l)/(double)l;
                for ( int k=1; k<=l; k++ )
                    B[l] += A[j][k]*pow(z0[j],l-k)*(double)binomCoeffs[l-1][k-1];
            }
        }
    }

    this->coeffMpole = B;
}

// Assigns a node's interaction list as indices to other nodes
void node2D::evalIList(){

    if (lvl > 2){
        // find parent's colleagues
        std::vector<node2D*> nList = parent->findColleagues();
    
        // extract interaction list from parent's neighbours;
        for ( int j=0; j<nList.size(); j++){
            if ( !nList[j]->isLeaf ){
                for ( int k=0; k<4; k++){
                    if (!isColleague(nList[j]->child[k]))
                        ilist.push_back(nList[j]->child[k]);
                }
            }
        }
    } 

    if ( !isLeaf ) 
        for ( int j=0; j<4; j++ )
            child[j]->evalIList();
}

/* Finds the local expansion coefficients due to all particles in a node's interaction list
void node2D::evalCoeffLocalExp( int p ){
 
    std::vector<std::vector<int> > binomCoeffs(2*p, std::vector<int>(p));
    for ( int l=1; l<=p; l++ )
        for ( int k=1; k<=p; k++ )
            binomCoeffs[l+k-1][k-1] = binomCoeff(l+k-1,k-1);
  
    std::vector<std::vector<std::complex<double> > > Bs(ilist.size());
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
}*/

// Assigns the local expansion coefficients of all members of a node's interaction list
// Also, allocates memory for node's own local expansion coefficients
void node2D::evalCoeffLocalExp( int p ){

 
    std::vector<std::vector<int> > binomCoeffs(2*p, std::vector<int>(p));
    for ( int l=1; l<=p; l++ )
        for ( int k=1; k<=p; k++ )
            binomCoeffs[l+k-1][k-1] = binomCoeff(l+k-1,k-1);
  
    std::vector<std::complex<double> > A = this->coeffMpole;
    std::vector<std::complex<double> > B(p+1);
    std::complex<double> z0;

    this->coeffLocalExp = B;
    for ( int i=0; i<ilist.size(); i++ ){
        z0 = (center - ilist[i]->center);

        B[0] = A[0]*log(-z0);
        for ( int k=1; k<=p; k++ )
            B[0] += A[k]*pow(-1,k)/pow(z0,k);

        for ( int l=1; l<=p; l++ ){
            B[l] = -A[0]/((double)l*pow(z0,l));
            for ( int k=1; k<=p; k++ )
                B[l] += A[k]*pow(-1,k)/pow(z0,l+k)*(double)binomCoeffs[l+k-1][k-1];
        }
        ilist[i]->coeffLocalExp = B;
    }

    if ( !isLeaf )
        for ( int j=0; j<4; j++ )
            child[j]->evalCoeffLocalExp(p);
}

// for a non-leaf, add the local expansion coefficients to those of its children
void node2D::evalCoeffLocalExpSum(int p){
    
    std::vector<std::vector<std::complex<double> > > B(4, std::vector<std::complex<double> >(p+1)); 
    std::vector<std::complex<double> > z0(4);

    if ( isLeaf ) return;
    else {
        for ( int j=0; j<4; j++ ){
            z0[j] = -(child[j]->center - center);
            B[j] = shiftLocalExp(z0[j]);
            for ( int k=0; k<p; k++ )
                child[j]->coeffLocalExp[k] += B[j][k];
        }

        for ( int j=0; j<4; j++ )
            child[j]->evalCoeffLocalExpSum(p);
    }
}

