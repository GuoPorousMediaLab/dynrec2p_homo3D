#ifndef simplegeometry_cell_h
#define simplegeometry_cell_h

#include <iostream>
#include "Node.h"
#include "Interface.h"
#include "SourceandSink.h"
#include "Const.h"
#include <vector>
using namespace std;

class Layer;
class Interface;

class Cell {
public:
    Cell(vector <int> nodeid, vector <int> interfaceid, int sourceandsinkid, double pb, double pc, double sc, double sb, double pcap, int ID, Layer* myLayer) 
    {create_(nodeid, interfaceid, sourceandsinkid, pb, pc, sc, sb, pcap, ID, myLayer);}
    Cell(){create_();} 
    void set_variableid(int variableid){variableid_ = variableid;};
    void set_cellnodes();
    void set_cellinterfaces();
    void set_cellsourceandsink();
    
    void initialize();
    void set_H();
    void set_center();
    void set_area(); // assume that the cell is regtangular
    void set_pbrine(double pb){dpbdt_ = (pb - pb_)/DT; pb_ = pb;}//update brine pressure
    void set_pbrinecoeffiandRHS();

    void verticalRec();
    
    int get_ID(){return ID_;}
    int get_variableid(){return variableid_;}
    vector<Interface*>::iterator Get_cellinterfacebegin(){return cellinterfaces_.begin();}
    vector<Interface*>::iterator Get_cellinterfaceend(){return cellinterfaces_.end();}
    SourceandSink* get_sourceandsink(){return cellsourceandsink_;}
    double get_centerx(){return centerx_;}
    double get_centery(){return centery_;}
    double get_centerz1(){return centerz1_;}
    double get_centerz2(){return centerz2_;}
    double get_H(){return H_;}
    double get_phi(){return phi_;}
    double get_area(){return area_;}
    double get_pbrine(){return pb_;}
    double get_dpbdt(){return dpbdt_;}
    double get_pco2(){return pc_;}
    double get_sb(){return sb_;}
    double get_sc(){return sc_;}
    double get_pbrinecoeffi(){return pbrinecoeffi_;}
    double get_pbrineRHS(){return pbrineRHS_;}
    double get_neighbourcellsc();
    double* get_permeability_local(){return &(permeability_local_.front());}
    double* get_mobilityb_local(){return &(mobilityb_local_.front());}
    double* get_mobilityc_local(){return &(mobilityc_local_.front());}
    double* get_fbrine_local(){return &(fbrine_local_.front());}
    double* get_pbrine_local(){return &(pbrine_local_.front());}
    double* get_sc_local(){return &(sc_local_.front());}
    double* get_pco2_local(){return &(pco2_local_.front());}
    double* get_pcap_local(){return &(pcap_local_.front());}
    double get_co2Mass();
    void output_nodeid();
    void output_interfaceid();
    
private:
    vector<Node*> cellnodes_;
    vector<int> nodeid_;
    vector<Interface*> cellinterfaces_;
    vector<int> interfaceid_;
    int sourceandsinkid_;
    SourceandSink* cellsourceandsink_;

    vector<double> pbrine_local_;
    vector<double> pco2_local_;
    vector<double> sb_local_;
    vector<double> sc_local_;
    vector<double> mobilityb_local_;
    vector<double> mobilityc_local_;
    vector<double> interMobilityb_local_;
    vector<double> interMobilityc_local_;
    vector<double> interMobility_local_;
    vector<double> arithMobilityb_local_;
    vector<double> arithMobilityc_local_;
    vector<double> arithMobility_local_;
    vector<double> fbrine_local_;
    vector<double> permeability_local_;
    vector<double> permeabilityz_local_;
    vector<double> pd_local_, phi_local_, pcap_local_, lambda_local_;
    vector<double> dpcapdt_local_, dfbdt_local_, dpbdt_local_, dpcdt_local_;

    double area_, H_, centerx_, centery_, centerz1_, centerz2_, phi_;
    double pb_, pc_, sb_, sc_; 
    double pbrinecoeffi_, pbrineRHS_, dpbdt_, ctc_, ctb_;
    int ID_, variableid_;
    Layer* myLayer_;
    void create_(vector<int>, vector<int>, int, double, double, double, double, double, int, Layer*); 
    void create_();
};

#endif
