#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "stdlib.h"
#include "stdio.h"
#include "Node.h"
#include "Cell.h"
#include "Interface.h"
#include "Layer.h"
#include "Model.h"

using namespace std;
class Layer;

void Model::create_(int Nx,int Ny)
{
    // create a layer with 4 x 4 nodes
    myLayer_.create_(Nx, Ny);
    //Initialize the interfaces in the layer.
    myLayer_.initialize();
}

void Model::create_(){
    
}

void Model::cal_Pressure(int Nx, int Ny){
    myLayer_.update_coeffi();
    // initialize the matrix and RHS first
    sparse_lib::sparse_sym_cds M(eqNum_, eqNum_, 3);
    sparse_lib::vec x(eqNum_);
    sparse_lib::vec b(eqNum_);
    cout << "Constructing M and b..." << flush; 
    myLayer_.set_matrixandRHS(M, b);
    cout << "[DONE]" << endl;
    
    if(eqNum_ < 100) {
        cout << "-------------------------------" << endl;
        cout << "M: " << endl;
        cout << scientific << M << endl;
    }
    
    if(eqNum_ < 100) {
        cout << "-------------------------------" << endl;
        cout << "b: " << endl;
        cout << b << endl;
    }

    sparse_lib::ICCPreconditioner<sparse_lib::sparse_sym_cds> ICCP;
    sparse_lib::Iterator<double> it(1e-50, 5000, 1e-50, sparse_lib::ITER_SUMMARY);
    sparse_lib::CG_Solver(M, x, b, ICCP, it, true);

    if(eqNum_ < 100) {
        cout << "-------------------------------" << endl;
        cout << "x: " << endl;
        cout << x << endl;
    }

    
    // ofstream myfile2;
    // myfile2.open ("mypressure.csv",ios::app);
    // if (myfile2.is_open()){
    //     myfile2 << x;
    // }
    // else cout << "failed to open file";
    // myfile2.close();
    
    // update the pressure with the pressure solution
    myLayer_.update_pressure(x);

}

void Model::cal_Saturation(){
    //go through interfaces to update CO2 flux
    myLayer_.update_flux();
    //go through cells to update saturation
    myLayer_.calculate_saturation();
}

void Model::massBalanceCheck(int n){
    myLayer_.massBalanceCheck(n);
}
