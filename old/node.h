#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

class node2D {
private:
    int lvl; // level of subdivision of this node
    int inode; // index of this node among all nodes in the tree
    int order;
    int isLeaf;

    double size; // length of largest dimension of this node
    std::complex<double> center; // coordinates of center of this node relative to origin 

    std::vector<int> iz; // index of all particles in Z living in this node
    std::vector<std::vector<int> > ilists; // indices to other nodes which comprise various i-lists
    std::vector<std::complex<double> > coeffMpole; // multipole coefficients due to all particles in this node
    std::vector<std::complex<double> > coeffLocalExp; // local expansion coefficients due to all particles in i-list

    node2D *parent; 
    node2D *branch1;
    node2D *branch2;
    node2D *branch3;
    node2D *branch4;

public:
    node2D(std::vector<std::complex<double> > Z, int maxparts, 
        std::vector<int> iz, int lvl, double size, std::complex<double> center);

    ~node2D() {}

    int numNodes();

    node2D* ithNode(int i);

    void evalInode(int i);

    void evalCoeffMpole(std::vector<std::complex<double> > Z, std::vector<double> Q, int p);

    std::vector<std::vector<std::complex<double> > > coeffMpoleAll();

    std::vector<std::complex<double> > centerAll();

    int isNeighbour(node2D* otherNode);

    /* void node2D::evalCoeffLocalExp(
        std::vector<std::vector<std::complex<double> > > Bs, std::vector<std::complex<double> >, int p);

    void node2D::evalCoeffLocalExpSum(); 
    */    

    void fprintZ(std::vector<std::complex<double> > Z);
};
