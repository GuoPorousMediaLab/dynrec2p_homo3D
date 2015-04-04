#ifndef simplegeometry_Layer_h
#define simplegeometry_Layer_h

#include <vector>
#include <sparselib/sparselib.h>
#include <sparselib/typedefs.h>
#include "Node.h"
#include "Cell.h"
#include "Interface.h"
#include "SourceandSink.h"
#include "Const.h"
using namespace std;

class Layer {
    
public:
    Layer(int Nx, int Ny) {create_(Nx, Ny);}
    Layer() {create_();}
    Node* findNodeByID(int);
    Cell* findCellByID(int);
    Interface* findInterfaceByID(int);
    SourceandSink* findSourceandSinkByID(int);
    void initialize();
    void update_coeffi();
    void update_pressure(sparse_lib::vec &b);
    void set_matrixandRHS(sparse_lib::sparse_sym_cds &M, sparse_lib::vec &b);
    
    void update_flux();
    void calculate_saturation();
    
    double get_Nx() const {return Nx_;}
    double get_Ny() const {return Ny_;}
    double get_krcstar(){return krcstar_;}
    double get_miub(){return miub_;}
    double get_miuc(){return miuc_;}
    double get_rhob(){return rhob_;}
    double get_rhoc(){return rhoc_;}
    double get_compressibilityc(){return compressibilityc_;}
    double get_compressibilityb(){return compressibilityb_;}
    double get_compressibilityr(){return compressibilityr_;}
    int get_verRefineNum(){return verRefineNum_;}
    void massBalanceCheck(int n);
    
    int get_i(){return i_;}

    void create_(int, int);
    void create_();
  
private:
    int Nx_;
    int Ny_;
    
    //just for tempory output
    int i_;
    int verRefineNum_;
    
    vector<Node> nodes_;
    vector<Cell> cells_;
    vector<Interface> interfaces_;
    vector<SourceandSink> sourceandsinks_;
   
    double krcstar_, srb_, miub_, miuc_, rhob_, rhoc_, compressibilityc_, compressibilityb_, compressibilityr_, phi_; 

};

#endif
