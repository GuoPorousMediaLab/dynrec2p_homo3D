#include <iostream>
#include <vector>
#include "Node.h"
#include "Layer.h"
#include "Model.h"
#include <ctime>

int main() {
    std::clock_t start;
    double duration;

    start = std::clock();

    // create a layer with 4 x 4 nodes
    int Nx = 601;
    int Ny = 2;
    Model mymodel(Nx, Ny);
    for (int i = 0; i < 365*500; i++) {
        mymodel.cal_Pressure(Nx,Ny);
        mymodel.cal_Saturation();
        mymodel.massBalanceCheck(i+1);
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"printf: "<< duration <<'\n';

    return (1);
}

