#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

class node2D {
private:
    int lvl; // level of subdivision of this node
    int inode; // index of this node among all nodes in the tree
    int order; // index of this node among its siblings
    int isLeaf;

    double size; // length of largest dimension of this node
    std::complex<double> center; // coordinates of center of this node relative to origin 

    std::vector<int> iz; // index of all particles in Z living in this node
    std::vector<std::complex<double> > coeffMpole; // multipole coefficients due to all particles in this node
    std::vector<std::complex<double> > coeffLocalExp; // local expansion coefficients due to all particles in i-list
    std::vector<node2D*> ilist; // pointers to other nodes which comprise various i-lists
 
    node2D* parent; 
    std::vector<node2D*> child;

public:
    // master node only
    // std::vector<std::vector<std::complex<double> > > coeffMpoleList;
    // std::vector<std::complex<double> > centerList;

    node2D(std::vector<std::complex<double> > Z, int maxparts,
        std::vector<int> iz, int lvl, int order, double size, std::complex<double> center);

    ~node2D() {}

    int numNodes();

    int isNeighbour(node2D* other);

    // node2D* findNode(int i);

    void fprintZ(std::vector<std::complex<double> > Z);

    void fprintMpole();

    void evalInode(int i);

    void evalCoeffMpole(std::vector<std::complex<double> > Z, std::vector<double> Q, int p);

    void evalIList();

    void evalCoeffLocalExp(int p);

    /*void node2D::evalCoeffLocalExpSum(); 
    */    

};
