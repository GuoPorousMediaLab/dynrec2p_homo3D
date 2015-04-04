#include <iostream>
#include <algorithm>
#include <string>
#include <sparselib/sparselib.h>
#include <sparselib/typedefs.h>
#include "Node.h"
#include "Cell.h"
#include "Layer.h"
#include "math.h"
#include "Const.h"
#include <iomanip>
#include <fstream>
using namespace std;

void Layer::create_(int Nx, int Ny) {
    Nx_ = Nx;
    Ny_ = Ny;
    
    // Initialize the constants
    krcstar_ = 1.0; miub_ = 3.0E-4, miuc_ = 4.25E-5; rhob_ = 1000.0; rhoc_ = 710.0; compressibilityc_ = 0.0E-9; compressibilityb_ = 0.0E-10; 
    compressibilityr_ = 0.0E-10; verRefineNum_ = 200;
    
    // for output
    i_ = 0;
    
    //create nx by ny nodes
    double centerx = 0.0;
    double centery = 0.0;
    double dx = 1.0*25;
    double dy = 1.0*100;
    double startx = centerx - (Nx_-1) / 2.0 * dx;
    double starty = centery - (Ny_-1) / 2.0 * dy;
    for (int i=1; i<=Ny_; i++){
        for (int j=1; j<=Nx_; j++){
            nodes_.push_back(Node(startx + (j-1)*dx, starty + (i-1)*dy, -1050, -1000.0, (i-1)*Nx_+j, this));
        }
    }

    for (int i = 0; i < Nx_ * Ny_; ++i)
    {
        cout << nodes_[i].get_ID() << endl;
    }
    
    //create cells
    vector<int> cellnodeid;
    vector<int> cellinterfaceid;
    int cellsourceandsinkid = 0;
    for (int i=1; i<Ny_; i++) {
        for (int j=1; j<Nx_; j++){
            cellnodeid.push_back((i-1)*Nx_+j);
            cellnodeid.push_back((i-1)*Nx_+j+1);
            cellnodeid.push_back(i*Nx_+j+1);
            cellnodeid.push_back(i*Nx_+j);  
            cellinterfaceid.push_back((i-1)*(Nx_-1)+j);
            cellinterfaceid.push_back(Ny_*(Nx_-1)+(i-1)*Nx_+j+1);
            cellinterfaceid.push_back(i*(Nx_-1)+j);
            cellinterfaceid.push_back(Ny_*(Nx_-1)+(i-1)*Nx_+j);
            cellsourceandsinkid = 0;
            // if (i == 1 && j == 1){
            //     cellsourceandsinkid = 1;
            // }
            if (j == 1)
            {
                cellsourceandsinkid = 1;
            }
            else if (j == Nx_ - 1)
            {
                cellsourceandsinkid = 2;
            }
            cells_.push_back(Cell(cellnodeid, cellinterfaceid, cellsourceandsinkid, 0.0, 0.0, 0.0, 1.0, 0.0, (i-1)*(Nx_-1)+j, this)); //pb, pc, sc, sb, pcap, id 
            cellnodeid.clear();
            cellinterfaceid.clear();
        }
    }

    //create interfaces
    //horizontal interfaces
    vector<int> interfacecellid;
    vector<int> interfacenodeid;
    int interfacetype;
    for (int i=1; i<=Ny_; i++) {
        for (int j=1; j<Nx_; j++) {
            interfacenodeid.push_back((i-1)*Nx_+j);
            interfacenodeid.push_back((i-1)*Nx_+j+1);
            if (i<Ny_) {
                interfacecellid.push_back((i-1)*(Nx_-1)+j);
            }
            if (i>1) {
                interfacecellid.push_back((i-2)*(Nx_-1)+j);
            }

            interfacetype = 0;
            if (i == Ny_ ){
                //flux boundary
                interfacetype = 1;
            }
            else if (i == 1 ){
                //pressure boundary
                interfacetype = 1;
            }
            interfaces_.push_back(Interface(interfacenodeid, interfacecellid, 0.0, 0.0, 0.0, 0.0, interfacetype, (Nx_-1)*(i-1)+j, this)); 
            // fluxb_, fluxc_, pb_, sc_
            interfacenodeid.clear();
            interfacecellid.clear();
        }
    }

    //vertical interfaces
    for (int i=1; i<Ny_; i++) {
        for (int j=1; j<=Nx_; j++) {
            interfacenodeid.push_back((i-1)*Nx_+j);
            interfacenodeid.push_back(i*Nx_+j);
            if (j<Nx_) {
                interfacecellid.push_back((i-1)*(Nx_-1)+j);
            }
            if (j>1) {
                interfacecellid.push_back((i-1)*(Nx_-1)+j-1);
            }

            interfacetype = 0;
            if (j == Nx_){
                //pressure boundary
                interfacetype = 1;
            }
            else if (j == 1){
                //no flow boundary
                interfacetype = 1;
            }

            interfaces_.push_back(Interface(interfacenodeid, interfacecellid, 0.0, 0.0, 0.0, 0.0, interfacetype, (Nx_-1)*Ny_+(i-1)*Nx_+j, this));
            // fluxb_, fluxc_, pb_, sc_
            interfacenodeid.clear();
            interfacecellid.clear();
        }
    }

 // Create sourceandsinks
    // first source
    sourceandsinks_.push_back(SourceandSink(1, 0.001, 0.0, this)); // id, qc, qb
    // second source 
    sourceandsinks_.push_back(SourceandSink(2, -0.001, 0.0, this));
   
    vector<Interface>::iterator interfaces_iter = interfaces_.begin();
    vector<Interface>::iterator interfaces_end = interfaces_.end();
    while (interfaces_iter != interfaces_end) {
//        cout << "interfacetype" <<(*interfaces_iter).get_interfacetype() << endl;
        (*interfaces_iter).set_interfacenodes();
        (*interfaces_iter).set_interfacecells();
        
        //checking numbering of the interface system
//        (*interfaces_iter).Output_nodeid();
//        (*interfaces_iter).Output_cellid();
        interfaces_iter++;
    }
    
    //setup variable id for cells
    vector<Cell>::iterator cells_iter = cells_.begin();
    vector<Cell>::iterator cells_end = cells_.end();
    int i = 0;
    while (cells_iter!=cells_end) {
        cells_iter->set_cellnodes();
        cells_iter->set_cellinterfaces();
        cells_iter->set_cellsourceandsink();
        cells_iter->set_variableid(i); // set the variableid
        //check numbering of the cell system
        //        cout <<((*cells_iter).Get_ID()) << endl;
                  (*cells_iter).output_nodeid();
        //        (*cells_iter).Output_interfaceid();
        
        i++; cells_iter++;
    } 
    
}

void Layer::create_() {
    Nx_ = 0;
    Ny_ = 0;
}

//find node according id in the layer
Node* Layer::findNodeByID(int nodeID){
    return &(nodes_[nodeID-1]);
}

//find cell according id in the layer
Cell* Layer::findCellByID(int cellID){
    return &(cells_[cellID-1]);
}

//find interface according id in the layer
Interface* Layer::findInterfaceByID(int interfaceID){
    return &(interfaces_[interfaceID-1]);
}

SourceandSink* Layer::findSourceandSinkByID(int sourceandsinkID){
    return &(sourceandsinks_[sourceandsinkID-1]);
}

void Layer::initialize(){
    
    vector<Cell>::iterator cell_iter = cells_.begin();
    vector<Cell>::iterator cell_end = cells_.end();
    while (cell_iter != cell_end) {
        cell_iter->initialize();
        cell_iter++;
    }    
    
    vector<Interface>::iterator interface_iter = interfaces_.begin();
    vector<Interface>::iterator interface_end = interfaces_.end();
    while (interface_iter != interface_end) {
        interface_iter->initialize();
        interface_iter++;
    }

}

//update the brine pressure coefficient
void Layer::update_coeffi(){
    vector<Cell>::iterator cell_iter = cells_.begin();
    vector<Cell>::iterator cell_end = cells_.end();
    while (cell_iter != cell_end) {
        // update coefficient with new mobilities
        cell_iter->set_pbrinecoeffiandRHS();
        cell_iter++;
    }   
    
    vector<Interface>::iterator interface_iter = interfaces_.begin();
    vector<Interface>::iterator interface_end = interfaces_.end();
    while (interface_iter != interface_end) {
        // if (interface_iter->get_interfacetype() == 0) {
        //     if ((*(interface_iter->get_cell()))->get_sc() >1.0E-20 || (*(interface_iter->get_cell()+1))->get_sc() >1.0E-20 || i_ < 10)
        //     {
        //         interface_iter->set_pbrinecoeffiandRHS();
        //     }
        // }
        // else interface_iter->set_pbrinecoeffiandRHS();
        interface_iter->set_pbrinecoeffiandRHS();
        interface_iter++;
    }
}

//update the brine pressure for every cell in this layer after the pressure calculation
void Layer::update_pressure(sparse_lib::vec &x){
    i_++;
    ofstream myfile4;
    //ofstream myfile5;
    vector<Cell>::iterator cell_iter = cells_.begin();
    vector<Cell>::iterator cell_end = cells_.end();
    int i = 0;
    while (cell_iter != cell_end) {
        cell_iter->set_pbrine(x(i));

        if (i_ % (365*50) == 0)
        {
            myfile4.open ("coarse_pbrine.csv",ios::app);
            if (myfile4.is_open()){
                myfile4 << cell_iter->get_ID()<<","<< cell_iter->get_pbrine() << endl;  
            }
            else cout << "failed to open file4" << endl;
            
            myfile4.close();

            // myfile5.open ("finescale_pco2.csv",ios::app);
            // if (myfile5.is_open()){
            //     myfile5 << cell_iter->get_ID() << "," << cell_iter->get_area() << ",";
            //     for (int j=0; j < verRefineNum_; j++){
            //         myfile5 << (*(cell_iter->get_pco2_local() + j)) - 1.1025E07 <<"," ;  
            //     }
            //         myfile5 << endl;
            //     }
            
            // else cout << "failed to open file5" << endl;
            
            // myfile5.close();
        }
        cell_iter++;
        i++;
    }

    if (i_ % (365*50) == 0)
    {
        myfile4.open ("coarse_pbrine.csv",ios::app);
        if (myfile4.is_open())
        {
            myfile4 << endl << "i_ is: "<< i_ << endl;
        }
        myfile4.close();

        // myfile5.open ("finescale_pco2.csv",ios::app);
        // if (myfile5.is_open())
        // {
        //     myfile5 << endl << "i_ is: "<< i_ << endl;
        // }
        // myfile5.close();
    }
    
}

//set the linear system matrix and RHS
void Layer::set_matrixandRHS(sparse_lib::sparse_sym_cds &M, sparse_lib::vec &b){
    vector<Interface>::iterator interface_iter = interfaces_.begin();
    vector<Interface>::iterator interface_end = interfaces_.end();
    while (interface_iter != interface_end) {
       /* cout << "interface id:" << interface_iter->get_ID() <<endl; */
        int i = (*(interface_iter->get_cell()))->get_variableid();
        b(i) += (*(interface_iter->get_pbrineRHS()));
        // cout <<"B"<<"["<<i<<"]"<<"is: "<< B[i]<< endl; 
        if (interface_iter->get_interfacetype()!= 1) {
            M(i,i) += (*(interface_iter->get_pbrinecoeffi()));
            // cout <<"A"<<"["<<i<<"]"<<"["<<i<<"]"<<"is: "<< A[i][i] << endl; 
        }
         
        if (interface_iter->get_interfacetype() == 0) {
            int j = (*(interface_iter->get_cell()+1))->get_variableid();
            // cout <<"cell ID: " << (*(interface_iter->get_cellpointer()+1))->get_ID() << endl; 
            M(j,j) += (*(interface_iter->get_pbrinecoeffi())); 
            if (i < j){
                M(i,j) = (*(interface_iter->get_pbrinecoeffi()+1));
            }
            else M(j,i) = (*(interface_iter->get_pbrinecoeffi()+1));
            b(j) += (*(interface_iter->get_pbrineRHS()+1));
        }
        interface_iter++;
    }
    
    // Go through cells to set coefficient from the time derivative term
    vector<Cell>::iterator cell_iter = cells_.begin();
    vector<Cell>::iterator cell_end = cells_.end();
    while (cell_iter != cell_end) {
        int i = cell_iter->get_variableid();
        M(i,i) += cell_iter->get_pbrinecoeffi();
        /*cout <<"A"<<"["<<i<<"]"<<"["<<i<<"]"<<"is: "<< A[i][i] << endl; */
        b(i) += cell_iter->get_pbrineRHS();
//        cout <<"contribution from time derivative term "<<"B"<<"["<<i<<"]"<<"is: "<< cell_iter->get_pbrineRHS()<< endl;
        cell_iter++;    
    }

    // fixing pressure at the right
    b(Nx_-2) = -9810*M(Nx_-2,Nx_-2);
    b(Nx_-3) -= -9810*M(Nx_-3,Nx_-2);
    // b(2*(Nx_-1)-1) = -9810*M(2*(Nx_-1)-1,2*(Nx_-1)-1);
    // b(2*(Nx_-1)-2) -= -9810*M(2*(Nx_-1)-2,2*(Nx_-1)-1);
    M(Nx_-3,Nx_-2) = 0.0;
    // M(2*(Nx_-1)-2,2*(Nx_-1)-1) = 0.0;
    // M(Nx_-2,2*(Nx_-1)-1) = 0.0;


    
}


void Layer::update_flux(){
    vector<Interface>::iterator interface_iter = interfaces_.begin();
    vector<Interface>::iterator interface_end = interfaces_.end();
    while (interface_iter != interface_end) { 
        // if (interface_iter->get_interfacetype() == 0)
        // {
        //     if ((*(interface_iter->get_cell()))->get_sc() >1.0E-40 || (*(interface_iter->get_cell()+1))->get_sc() >1.0E-40 || i_ < 10)
        //     {
        //         interface_iter->set_flux();
        //     }
        // }
        // else
        // {
        //     if ((*(interface_iter->get_cell()))->get_sc() >1.0E-40 || i_ < 10)
        //     {
        //         interface_iter->set_flux();
        //     }
        // } 
        
        interface_iter->set_flux();
        interface_iter++;
    }
}

void Layer::calculate_saturation(){
    // Calculate the saturation from time derivative term and flux
    cout << i_ << endl;
    //ofstream myfile1;
    // ofstream myfile2;
    ofstream myfile3;
    
    vector<Cell>::iterator cell_iter = cells_.begin();
    vector<Cell>::iterator cell_end = cells_.end();
    while (cell_iter != cell_end) {

        // if (cell_iter->get_neighbourcellsc() > 1.0E-40 || i_ < 10) {
        //     cell_iter->verticalRec();
        // }
        cell_iter->verticalRec();
        if (i_ % (365*50) == 0 && cell_iter->get_sc() > 1.0E-40)
         {
            
    //         myfile1.open ("finescale_sco2.csv",ios::app);
    //         if (myfile1.is_open()){
    //             myfile1 << cell_iter->get_ID() << "," << cell_iter->get_area() << ",";
    //             for (int j=0; j < verRefineNum_; j++){
    //                 myfile1 << (*(cell_iter->get_sc_local() + j)) <<"," ;  
    //             }
    //                 myfile1 << endl;
    //             }
            
    //         else cout << "failed to open file1" << endl;
            
    //         myfile1.close();
            
            
    //         // myfile2.open ("pbrine_local.csv",ios::app);
    //         // if (myfile2.is_open()){
    //         //     for (int j=0; j < verRefineNum_; j++){
    //         //         myfile2 << (*(cell_iter->get_pbrine_local() + j)) <<"," ;  
    //         //     }
    //         //     myfile2 << endl;
    //         // }
            
    //         // else cout << "failed to open file2" << endl;
            
    //         // myfile2.close();
            
            myfile3.open ("coarse_sco2.csv",ios::app);
            if (myfile3.is_open()){
                myfile3 << cell_iter->get_ID() << "," << cell_iter->get_area() << "," << cell_iter->get_sc() << endl;  
            }
            else cout << "failed to open file3" << endl;
            
            myfile3.close();
         }    
         cell_iter++;     
     }

     if (i_ % (365) == 0 )
     {
    //     myfile1.open ("finescale_sco2.csv",ios::app);
    //     if (myfile1.is_open())
    //     {
    //         myfile1 << endl << "i_ is: "<< i_ << endl;
    //     }
    //     myfile1.close();
    //     // myfile2.open ("pbrine_local.csv",ios::app);
    //     // if (myfile2.is_open())
    //     // {
    //     //     myfile2 << endl << "i_ is: "<< i_ << endl;
    //     // }
    //     // myfile2.close();

        myfile3.open ("coarse_sco2.csv",ios::app);
        if (myfile3.is_open())
        {
            myfile3 << endl << "i_ is: "<< i_ << endl;
        }
        myfile3.close();
     }
    
}

void Layer::massBalanceCheck(int n){
    ofstream myfile5;
    int cellSize = cells_.size();
    double co2Mass = 0.0;
    for (int i = 0; i < cellSize; ++i)
    {
        co2Mass += cells_[i].get_co2Mass();
    }
    myfile5.open("massBalanceCheck.csv",ios::app);
    if (myfile5.is_open() && n % (365*50) == 0)
    {
        myfile5 << "Time step is: " << n << endl;
        myfile5 << "co2Mass in the system: " << "," << co2Mass <<","<< " injected amount: " <<","<< 0.001 * DT * n <<","<<" difference: "<<","<< co2Mass - 0.001 * DT * n << endl;
    }
    myfile5.close();
}
