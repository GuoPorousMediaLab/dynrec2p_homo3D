
#ifndef Simplegeometry_Interface_h
#define Simplegeometry_Interface_h

#include <iostream>
#include <vector>
#include "Node.h"
#include "Cell.h"
using namespace std;

class Layer;
class Cell;

class Interface
{
public:
    Interface ( vector<int> nodeid, vector<int> cellid, double fluxb, double fluxc, double pb, double sc, int interfacetype, int id, Layer* MyLayer) { create_(nodeid, cellid, fluxb, fluxc, pb, sc, interfacetype, id, MyLayer); }
    Interface (){create_();}
    void initialize();// to initialize some of the variables which don't need to be updated
    int get_ID(){return ID_;}
    void set_interfacenodes();
    void set_interfacecells();
    void set_flux();
    void set_Tb();
    void set_Tc();
    Cell* get_cell(int cellid);
    void set_pbrinecoeffiandRHS();
    

    double get_interfacecenterx();
    double get_vectorprojection(){return vectorprojection_; }
    double* get_pbrinecoeffi(){return &(pbrinecoeffi_.front());} //get the pointer which points to the first elment of the coeffi vector
    double* get_pbrineRHS(){return &(pbrineRHS_.front());}//get the pointer which points to the first elment of the RHS vector
    int get_interfacetype(){return interfacetype_;}
    Cell** get_cell(){return &(interfacecells_.front());}
    int get_firstcellid(){return cellid_[0];}
    double* get_fluxc_local(){return &(fluxc_local_.front());}
    double* get_fluxb_local(){return &(fluxb_local_.front());}
    void output_nodeid();
    void output_cellid();
    
private:
    vector<int> nodeid_;
    vector<Node*> interfacenodes_;
    vector<int> cellid_;
    vector<Cell*> interfacecells_;
    vector<double> pbrinecoeffi_; //store the coeffi as a vector, the order is ii(jj)(ii and jj are the same), ij. i is the first cell, j is the second one.
    vector<double> pbrineRHS_; //the order is the same as the interfacecells. 
    double pb_, sc_, fluxb_, fluxc_, kapa_, permeability_, kmobilityb_, kmobilityc_;
    double centerx_, centery_, centerz1_, centerz2_, H_, vectorprojection_, centerdistance_, lengthAB_;
    vector<double> pd_local_, pb_local_, pcap_local_, fluxb_local_, fluxc_local_, mobilityb_local_, mobilityc_local_, permeability_local_;
    int ID_, interfacetype_;  // 1 is dirichlet, -1 is neumann, 0 is interior
    
    Layer* myLayer_;
    
    void create_(vector<int>, vector<int>, double, double, double, double, int, int, Layer*);
    void create_();
    
    
};


#endif
