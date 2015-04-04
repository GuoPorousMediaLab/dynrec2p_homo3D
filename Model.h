#ifndef Simplegeometry_Model_h
#define Simplegeometry_Model_h

#include "sparselib/sparselib.h"
#include "sparselib/typedefs.h"
#include "Node.h"
#include "Interface.h"
#include "Cell.h"
#include "Layer.h"
#include <vector>
using namespace std;

class Model {
    
public:
    Model(int Nx, int Ny):eqNum_((Nx - 1)*(Ny - 1)){ create_(Nx, Ny);}
    Model():eqNum_(0){create_();};
    void cal_Pressure(int Nx, int Ny);
    void cal_Saturation();
    void massBalanceCheck(int n);
private:
    sparse_lib::sparse_sym_cds M_;
    sparse_lib::vec x_;
    sparse_lib::vec b_;
    const int eqNum_;
    void create_(int, int);
    void create_();
    Layer myLayer_;
        
};

#endif
