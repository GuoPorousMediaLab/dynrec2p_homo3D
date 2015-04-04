#include "Cell.h"
#include "Node.h"
#include "Layer.h"
#include "Const.h"
#include "Material.h"
#include "SourceandSink.h"
#include <math.h>
using namespace std;

void Cell::create_(vector <int> nodeid, vector <int> interfaceid, int sourceandsinkid, double pb, double pc, double sc, double sb, double pcap, int ID, Layer* myLayer) 
{    
    nodeid_ = nodeid;
    interfaceid_ = interfaceid;
    sourceandsinkid_ = sourceandsinkid;
    pb_ = pb;
    pc_ = pc;
    sc_ = sc;
    sb_ = sb;
    ID_ = ID; 
    myLayer_ = myLayer;
    phi_ = 0.15;
    
    for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
    {
        //topdepth = 1000m, pressure gradient = 10.5
        pd_local_.push_back(0.0E4);
        pbrine_local_.push_back( - i * H_/(myLayer_->get_verRefineNum())*myLayer_->get_rhob() * G);
        pco2_local_.push_back( - i * H_/(myLayer_->get_verRefineNum())*myLayer_->get_rhob() * G + pd_local_[i] );
        sb_local_.push_back(1.0);
        sc_local_.push_back(0.0);
        mobilityb_local_.push_back(1.0/myLayer_->get_miub());
        mobilityc_local_.push_back(0.0/myLayer_->get_miuc());
        fbrine_local_.push_back(0.0);
        permeability_local_.push_back(1.0E-14);
        permeabilityz_local_.push_back(1.0E-14);
        phi_local_.push_back(0.15);
        lambda_local_.push_back(2.0);
        pcap_local_.push_back(pd_local_[i]);
        dpcapdt_local_.push_back(0.0);
        dfbdt_local_.push_back(0.0);
        dpbdt_local_.push_back(0.0);
        dpcdt_local_.push_back(0.0);
    }
    for (int i = 0; i < (myLayer_->get_verRefineNum() - 1); ++i)
    {
        interMobilityb_local_.push_back(1.0/myLayer_->get_miub());
        interMobilityc_local_.push_back(0.0/myLayer_->get_miuc());
        interMobility_local_.push_back(0.0);
        arithMobilityb_local_.push_back(1.0/myLayer_->get_miub());
        arithMobilityc_local_.push_back(0.0/myLayer_->get_miuc());
        arithMobility_local_.push_back(0.0);
    }

    for (int i = 0; i < myLayer_->get_verRefineNum()/2; ++i)
    {
        permeability_local_[i] = 1.0E-14;
        permeabilityz_local_[i] = 1.0E-14;
    }
    
}                       

void Cell::create_()
{
    pb_ = 0.0;
    pc_ = 0.0;
    sc_ = 0.0;
    sb_ = 0.0;
    ID_ = 0;
}


void Cell::set_cellnodes()
{
    vector<int>::iterator nodeid_iter = nodeid_.begin();
    vector<int>::iterator nodeid_end = nodeid_.end();
    Node* temp;
    while (nodeid_iter != nodeid_end) {
        temp = (*myLayer_).findNodeByID(*nodeid_iter);
        if (temp == NULL) {
            cout << "can't find node according the given id" << endl;
            break;
        }
        cellnodes_.push_back(temp);
        nodeid_iter++;
    }
    
}

void Cell::set_cellinterfaces()
{
    vector<int>::iterator interfaceid_iter = interfaceid_.begin();
    vector<int>::iterator interfaceid_end = interfaceid_.end();
    Interface* temp;
    while (interfaceid_iter != interfaceid_end) {
        temp = (*myLayer_).findInterfaceByID(*interfaceid_iter);
        if (temp == NULL) {
            cout << "can't find interface according the given id" << endl;
            break;
        }
        cellinterfaces_.push_back(temp);
        
        interfaceid_iter++;
    }
    
}

void Cell::set_cellsourceandsink()
{
    if (sourceandsinkid_ != 0)
    {
        cout <<"sourceandsinkid_: "<<sourceandsinkid_ << endl;
        cellsourceandsink_ = (*myLayer_).findSourceandSinkByID(sourceandsinkid_);
    }
    else {cellsourceandsink_ = NULL;}
}


void Cell::set_area(){
    double side1, side2;
    side1 = sqrt(pow(cellnodes_[0]->get_x() - cellnodes_[1]->get_x(), 2) + pow(cellnodes_[0]->get_y() - cellnodes_[1]->get_y(), 2));
    side2 = sqrt(pow(cellnodes_[1]->get_x() - cellnodes_[2]->get_x(), 2) + pow(cellnodes_[1]->get_y() - cellnodes_[2]->get_y(), 2));
    area_ = side1*side2;
}

void Cell::set_center(){
    double temp1 = 0; double temp2 = 0; double temp3 = 0; double temp4 = 0;
    for (int i = 0; i< cellnodes_.size(); i++){
        temp1 += cellnodes_[i]->get_x();
        temp2 += cellnodes_[i]->get_y();
        temp3 += cellnodes_[i]->get_z1();
        temp4 += cellnodes_[i]->get_z2();
    }
    centerx_ = temp1/cellnodes_.size();
    centery_ = temp2/cellnodes_.size();
    centerz1_ = temp3/cellnodes_.size();
    centerz2_ = temp4/cellnodes_.size();
}

void Cell::set_H(){
    double z1 = 0; 
    double z2 = 0; 
    for (int i = 0; i < cellnodes_.size(); i++) {
        z1 += cellnodes_[i]->get_z1();
        z2 += cellnodes_[i]->get_z2();
    }
    H_ = (z2-z1)/cellnodes_.size();
}

void Cell::initialize(){
    double phi = 0;
    for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
    {
        phi += phi_local_[i];
    }
    
    phi = phi/(myLayer_->get_verRefineNum());
    ctc_ = phi*myLayer_->get_compressibilityc() + myLayer_->get_compressibilityr();
    ctb_ = phi*myLayer_->get_compressibilityb() + myLayer_->get_compressibilityr();
    dpbdt_ = 0.0;
    set_H();
    set_center();
    set_area();

    //double S1 = 0.01, S2 = 0.99;
    int m = myLayer_->get_verRefineNum();
    for (int j = 0; j < m; ++j)
    {
        //sc_local_[j] = S2 + (S1 - S2)*0.5*(1+tanh((j + 0.5)/m - 0.5))*0.5*(1-tanh((j + 0.5)/m-10000.0));
        sc_local_[j] = 0.0;
        sb_local_[j] = 1.0 - sc_local_[j];
        pcap_local_[j] = pcap(sc_local_[j], pd_local_[j], lambda_local_[j]);
        mobilityb_local_[j] = mobilitybrine(sb_local_[j], lambda_local_[j]);
        mobilityc_local_[j] = mobilityco2(sc_local_[j], lambda_local_[j]);
    }

}

void Cell::set_pbrinecoeffiandRHS(){
    pbrinecoeffi_ = ( sc_*ctc_ + sb_*ctb_)/DT*H_*area_;
    // time dirivative term and source and sink
    double tempdpcapdt = 0.0;
    double tempdfbdt = 0.0;
    for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
    {
        tempdfbdt += ( ctc_ * sc_local_[i] + ctb_ * sb_local_[i] ) * dfbdt_local_[i]* H_/myLayer_->get_verRefineNum();
        tempdpcapdt += ctc_ * sc_local_[i] * dpcapdt_local_[i] * H_/myLayer_->get_verRefineNum();
    }
    pbrineRHS_ = pb_*pbrinecoeffi_ - area_ * (tempdfbdt + tempdpcapdt);
    if (sourceandsinkid_ != 0)
    {
        /*
         cout << "reached here 3" << endl;
         */
        //  if (myLayer_->get_i() > 365*50)
        // {
        //     pbrineRHS_ += 0.0;
        // }
        // else pbrineRHS_ += cellsourceandsink_->get_co2source() + cellsourceandsink_->get_brinesource();
        pbrineRHS_ += cellsourceandsink_->get_co2source() + cellsourceandsink_->get_brinesource();
    }
}

double Cell::get_neighbourcellsc(){
    double temp = 0.0;
    for (int i=0; i < cellinterfaces_.size(); i++) {
        if (cellinterfaces_[i]->get_interfacetype() == 0) {
            temp += ((*(cellinterfaces_[i]->get_cell()))->get_sc() + (*(cellinterfaces_[i]->get_cell()+1))->get_sc());
        }
        else if (cellinterfaces_[i]->get_interfacetype() == -1)
        {
            temp += (*cellinterfaces_[i]->get_cell())->get_sc();
        }
    }
    return temp + sc_;
}

void Cell::verticalRec(){
    int m = myLayer_->get_verRefineNum();
    double dz = H_/m, deltasat0;
    double deltarho = myLayer_->get_rhob() - myLayer_->get_rhoc();
    double fluxc[m-1], fluxb[m-1];
    double interPermeability[m-1];
    
    for (int j = 0; j < m - 1; ++j)
    {
        fluxc[j] = 0.0;
        fluxb[j] = 0.0;
        interPermeability[j] = 2.0/(1.0/permeabilityz_local_[j] + 1.0/permeabilityz_local_[j+1]);
    } 
    
    double horiFlux[m], horiFluxc[m], sourceandsink[m], sourceandsinkc[m], totalflux[m+1];
    for (int i = 0; i < m; ++i)
    {
        horiFlux[i] = 0.0;
        horiFluxc[i] = 0.0;
        if (sourceandsinkid_ == 1)
        {
            sourceandsink[i] = (cellsourceandsink_->get_co2source() + cellsourceandsink_->get_brinesource())/m;
            sourceandsinkc[i] = cellsourceandsink_->get_co2source()/m;
            // sourceandsink[0] = (cellsourceandsink_->get_co2source() + cellsourceandsink_->get_brinesource());
            // sourceandsinkc[0] = cellsourceandsink_->get_co2source();
            // sourceandsink[0] = 0.0;
            // sourceandsinkc[0] = 0.0;
            // sourceandsink[1] = 0.0;
            // sourceandsinkc[1] = 0.0;
            // sourceandsink[2] = 0.0;
            // sourceandsinkc[2] = 0.0;
        }
        else if (sourceandsinkid_ == 2) { 
            sourceandsink[i] = (cellsourceandsink_->get_co2source() + cellsourceandsink_->get_brinesource())/m;
            sourceandsinkc[i] = cellsourceandsink_->get_co2source()/m * sc_local_[i];
            //cout << "sourceandsink " << sourceandsink[i] << endl;
            // sourceandsink[m-1] = 0.0;
            // sourceandsinkc[m-1] = 0.0;
            // sourceandsink[m-2] = 0.0;
            // sourceandsinkc[m-2] = 0.0;
            // sourceandsink[m-3] = 0.0;
            // sourceandsinkc[m-3] = 0.0;
        }
        else
        {
            sourceandsink[i] = 0.0;
            sourceandsinkc[i] = 0.0;
        }
    }
    vector<Interface*>::iterator interface_iter = cellinterfaces_.begin();
    vector<Interface*>::iterator interface_end = cellinterfaces_.end();
    while (interface_iter != interface_end) {
        // internal interface, flux is from cell[0] to cell[1]
        if ((*interface_iter)->get_interfacetype() == 0) {
            if ((*interface_iter)->get_firstcellid() != ID_ ) {
                for (int i = 0; i < m; ++i)
                {
                    horiFlux[i] += - (*((*interface_iter)->get_fluxc_local() + i)) - (*((*interface_iter)->get_fluxb_local() + i));
                    horiFluxc[i] += - (*((*interface_iter)->get_fluxc_local() + i));
                    
                    /*
                     cout << horiFlux[i] << endl;
                     cout << horiFluxc[i] << endl;
                     */
                }
                
            }
            else {
                for (int i = 0; i < m; ++i)
                {
                    horiFlux[i] += (*((*interface_iter)->get_fluxc_local() + i)) + (*((*interface_iter)->get_fluxb_local() + i));
                    horiFluxc[i] += (*((*interface_iter)->get_fluxc_local() + i));
                    
                    /*
                     cout << horiFlux[i] << endl;
                     cout << horiFluxc[i] << endl;
                     */
                }
            }
            
        }
        // pressure boundary, flux is from cell to boundary
        else if ((*interface_iter)->get_interfacetype() == -1){
            for (int i = 0; i < m; ++i)
            {
                horiFlux[i] += (*((*interface_iter)->get_fluxc_local() + i))+ (*((*interface_iter)->get_fluxb_local() + i));
                horiFluxc[i] += (*((*interface_iter)->get_fluxc_local() + i));
                /*
                 cout << horiFlux[i] << endl;
                 cout << horiFluxc[i] << endl;
                 */
            }
        }
        // flux bounday, flux is going into the cell
        else if ((*interface_iter)->get_interfacetype() == 1){
            for (int i = 0; i < m; ++i)
            {
                horiFlux[i] += - (*((*interface_iter)->get_fluxc_local() + i)) - (*((*interface_iter)->get_fluxb_local() + i));
                horiFluxc[i] += - (*((*interface_iter)->get_fluxc_local() + i));
                /*
                 cout << horiFlux[i] << endl;
                 cout << horiFluxc[i] << endl;
                 */
            }
        } 
        else cout << "interfacetype is out of scope" << endl;
        interface_iter++;
    }
    /*
     for (int i = 0; i < m; ++i)
     {
     cout << horiFlux[i] << endl;
     }
     */
    
    // Calculate total fluxes on each horizontal interfaces of the fine cells
    for (int i = 0; i < m+1; ++i)
    {
        totalflux[i] = 0.0;
    }
    for (int i = 1; i < m+1 ; ++i)
    {
        totalflux[i] = totalflux[i-1] + ( sourceandsink[i-1] - horiFlux[i-1] * dz - (ctb_*sb_local_[i-1]*(dpbdt_ + dfbdt_local_[i-1]) + ctc_*sc_local_[i-1]*(dpbdt_ + dfbdt_local_[i-1] + dpcapdt_local_[i-1])) * area_ * dz)/area_; 
    }

    double verticalmassbalance = totalflux[m-1] + ( sourceandsink[m-1] - horiFlux[m-1] * dz - (ctb_*sb_local_[m-1]*(dpbdt_ + dfbdt_local_[m-1]) + ctc_*sc_local_[m-1]*(dpbdt_ + dfbdt_local_[m-1] + dpcapdt_local_[m-1])) * area_ * dz)/area_;;
    double horizontalmassbalance = 0.0;

    for (int i = 0; i < m; ++i)
    {
        horizontalmassbalance += ( sourceandsink[i] - horiFlux[i] * dz)/area_;
    }
    
    // Check mass balance
    //if (abs(totalflux[m-1] * area_ + sourceandsink[m-1] - horiFlux[m-1] * dz - (ctb_*sb_local_[i-1]*dpbdt_ + ctc_*sc_local_[i-1]*dpbdt_) * area_ * dz) > 1.0E-6 * area_)
    //{
        
         //cout <<"away from mass balance: "<< tempdistri << endl;
         
        //cout << "There is a mass balance problem in vertical collumn" << endl;
    //}
    //for (int i = 0; i < m-1; ++i)
    //{
        //totalflux[i] += tempdistri/m;
    //}
    //totalflux[m-1] -= tempdistri/m*(m-1);

    // Calculate co2 saturation
    for (int j = 0; j < m; j++) {
        // cout << sco2_local_[j] << " ";
        double tempdeltafluxc = 0.0;
        if (j == 0) {
            tempdeltafluxc = interMobilityc_local_[j]/(interMobilityb_local_[j] + interMobilityc_local_[j]) * totalflux[j+1] 
            + interPermeability[j] * (interMobility_local_[j] * deltarho*G - arithMobility_local_[j] * (pcap_local_[j+1] - pcap_local_[j])/dz);
            sc_local_[j] = (- tempdeltafluxc * area_ + sourceandsinkc[j] - horiFluxc[j] * dz - 0.0 * ctc_*sc_local_[j]*(dpbdt_ + dpcapdt_local_[j]) * area_ * dz)*DT/(area_ * dz * phi_local_[j]) 
                            + sc_local_[j];
        }
        else if ( j == m-1){
            tempdeltafluxc = - interMobilityc_local_[j-1]/(interMobilityb_local_[j-1] + interMobilityc_local_[j-1]) * totalflux[j] 
            - interPermeability[j-1] * (interMobility_local_[j-1] * deltarho*G - arithMobility_local_[j-1] * (pcap_local_[j] - pcap_local_[j-1])/dz);
            sc_local_[j] = (- tempdeltafluxc * area_ + sourceandsinkc[j] - horiFluxc[j] * dz - 0.0 * ctc_*sc_local_[j]*(dpbdt_ + dpcapdt_local_[j]) * area_ * dz)*DT/(area_ * dz * phi_local_[j]) 
                            + sc_local_[j];
        }
        else {
            tempdeltafluxc = interMobilityc_local_[j]/(interMobilityb_local_[j] + interMobilityc_local_[j]) * totalflux[j+1] 
            - interMobilityc_local_[j-1]/(interMobilityb_local_[j-1] + interMobilityc_local_[j-1]) * totalflux[j]
            + interPermeability[j] * (interMobility_local_[j] * deltarho*G - arithMobility_local_[j] * (pcap_local_[j+1] - pcap_local_[j])/dz) 
            - interPermeability[j-1] * (interMobility_local_[j-1] * deltarho*G - arithMobility_local_[j-1] * (pcap_local_[j] - pcap_local_[j-1])/dz);
            sc_local_[j] = (- tempdeltafluxc * area_ + sourceandsinkc[j] - horiFluxc[j] * dz - 0.0 * ctc_*sc_local_[j]*(dpbdt_ + dpcapdt_local_[j]) * area_ * dz)*DT/(area_ * dz * phi_local_[j])                      
                            + sc_local_[j];
        }
        
    }
    // jump the transport calculation
    // for (int j = 0; j < m; j++)
    // {
    //     sc_local_[j] = 0;
    // }


    for (int j = 0; j < m-1; j++) {
        fluxc[j] = interMobilityc_local_[j]/(interMobilityb_local_[j] + interMobilityc_local_[j]) * totalflux[j+1] + interPermeability[j] * interMobility_local_[j] * (deltarho*G-(pcap_local_[j+1] - pcap_local_[j])/dz);
        fluxb[j] = totalflux[j+1] - fluxc[j];
    }

    //update brine saturaton, capillary pressure, mobility
    for (int j = 0; j < m; j++) 
    {   
        sb_local_[j] = 1.0 - sc_local_[j];
        double temppcap = pcap_local_[j];
        pcap_local_[j] = pcap(sc_local_[j], pd_local_[j], lambda_local_[j]);
        dpcapdt_local_[j] = (pcap_local_[j] - temppcap)/DT;
        mobilityb_local_[j] = mobilitybrine(sb_local_[j], lambda_local_[j]);
        mobilityc_local_[j] = mobilityco2(sc_local_[j], lambda_local_[j]);
    }

    //upstream weighting scheme
    for (int j = 0; j < m-1; ++j)
    {
        double GB = - interPermeability[j]*(myLayer_->get_rhoc() - myLayer_->get_rhob())*G;
        //cout << GB << endl;
        // if GB and ut have the same sign
        if ( totalflux[j+1] >= 0.0)
        {
            interMobilityc_local_[j] = mobilityc_local_[j];
            double thetab = totalflux[j+1] - interMobilityc_local_[j]*GB;
            if (thetab >= 0.0)
            {
                 interMobilityb_local_[j] = mobilityb_local_[j];
            } 
            else interMobilityb_local_[j] = mobilityb_local_[j+1];
        }
        else {
            interMobilityb_local_[j] = mobilityb_local_[j+1];
            double thetac = totalflux[j+1] + interMobilityb_local_[j]*GB;
            if (thetac >= 0.0)
            {
                interMobilityc_local_[j] = mobilityc_local_[j];
            }
            else interMobilityc_local_[j] = mobilityc_local_[j+1];
        }  
        interMobility_local_[j] = (interMobilityb_local_[j] * interMobilityc_local_[j])/(interMobilityb_local_[j] + interMobilityc_local_[j]);
    }

    for (int j = 0; j < m-1; ++j)
    {
        arithMobilityb_local_[j] = (mobilityb_local_[j] + mobilityb_local_[j+1])/2.0;
        arithMobilityc_local_[j] = (mobilityc_local_[j] + mobilityc_local_[j+1])/2.0;
        arithMobility_local_[j] = (arithMobilityb_local_[j] * arithMobilityc_local_[j])/(arithMobilityb_local_[j] + arithMobilityc_local_[j]);
    }

    //Update local interface mobility, by upstream scheme
    // for (int j = 0; j < m-1; j++) {
    //     if (fluxb[j] <= 0.0) {
    //         interMobilityb_local_[j] = mobilityb_local_[j+1];
    //     }
    //     else {
    //         interMobilityb_local_[j] = mobilityb_local_[j];
    //     }
        
    //     if (fluxc[j] <= 0.0) {
    //         interMobilityc_local_[j] = mobilityc_local_[j+1];
    //     }
        
    //     else {
    //         interMobilityc_local_[j] = mobilityc_local_[j];
    //     }
        
    //     interMobility_local_[j] = (interMobilityc_local_[j]*interMobilityb_local_[j]/(interMobilityc_local_[j] + interMobilityb_local_[j]));
    // }      
    // update brine and co2 flux, calculate interface mobility
    //for (int j = 0; j < m-1; j++) {
        // if ((deltarho*G - (pcap_local_[j+1] - pcap_local_[j])/dz) >= 0.0) {
        //     interMobilityc[j] = mobilityc_local_[j];
        //     interMobilityb[j] = mobilityb_local_[j+1];
        // }
        // else {
        //     interMobilityc[j] = mobilityc_local_[j+1];
        //     interMobilityb[j] = mobilityb_local_[j];
        // }
        // interMobility[j] = (interMobilityc[j]*interMobilityb[j]/(interMobilityc[j] + interMobilityb[j]));
        
        //fluxc[j] = interMobilityc[j]/(interMobilityb[j] + interMobilityc[j]) * totalflux[j+1] + interPermeability[j] * interMobility[j] * (deltarho*G-(pcap_local_[j+1] - pcap_local_[j])/dz);
        //fluxb[j] = totalflux[j+1] - fluxc[j];
        //cout << "totalflux" << totalflux[j+1] << endl;
        //cout << "fluxb" << fluxb[j] << endl;
        //cout <<"fluxb is: "<< fluxb[j] << endl;
        // if (fluxb[j] > 0.0)
        // {
        //     interMobilityb[j] = mobilityb_local_[j];
        // }
        // else {
        //     interMobilityb[j] = mobilityb_local_[j+1];
        // }
        // cout << interPermeability[i] << endl;
    //}
    
    // for (int i = 0; i < m+1; ++i)
    // {
    //     totalflux[i] = 0.0;
    // }
    
    // reconstruct local brine and co2 pressure  
    // dpbdt_local_[0] = (pb_ - pbrine_local_[0])/DT;  
    //pbrine_local_[0] = pb_;
    // for (int j = 0; j < m-1 ; j++) {
        // double temppb = pbrine_local_[j+1];
        //if (abs(interMobilityb[j]) > 0.0)
        //{
            //cout << interMobilityb[j] << endl;
          //  pbrine_local_[j+1] = pbrine_local_[j] - ( totalflux[j+1]/ interPermeability[j]
            //                                        + G*(interMobilityb_local_[j]*myLayer_->get_rhob() + interMobilityc_local_[j]*myLayer_->get_rhoc()))/(interMobilityb_local_[j] + interMobilityc_local_[j])*dz
              //                                      + arithMobilityc_local_[j]*(pcap_local_[j+1]-pcap_local_[j])/(arithMobilityb_local_[j] + arithMobilityc_local_[j]);

        //cout << totalflux[j+1]/ interPermeability[j] << endl;

        //}
        //else pbrine_local_[j+1] = pbrine_local_[j] - myLayer_->get_rhob() * G * dz;
        //dpbdt_local_[j+1] = dpbdt_local_[j];
        // dpbdt_local_[j+1] = (pbrine_local_[j+1] - temppb)/DT;
    // }
    // for (int i = 0; i < m; ++i)
    // {
        // dpcdt_local_[i] = dpbdt_local_[i] + dpcapdt_local_[i];
    // }
    
    // Calculate coarse scale saturation
    double temp1 = 0;
    double temp2 = 0;
    double temp3 = 0;
    for (int i = 0; i < m; ++i)
    {
        temp1 += phi_local_[i] * sc_local_[i] * dz;
        temp2 += phi_local_[i] * sb_local_[i] * dz;
        temp3 += phi_local_[i] * dz;
    }
    sc_ = temp1/temp3;
    sb_ = temp2/temp3;
    if (abs(sc_ + sb_ - 1.0) > 1.0E-6)
    {
        cout << "ERROR: sc_ + sb != 1.0 after reconstruction" << endl;
    }
    
    // double tempi = 0;
    fbrine_local_[0] = 0;
    // dfbdt_local_[0] = 0;
    // pcap_ = pcap_local_[m/2];
    // pbrine_local_[m/2] = pb_;
    pbrine_local_[0] = pb_;
    // fbrine_local_[m/2] = 0;
    // fco2_local_[m/2] = 0;
    // pco2_local_[m/2] = pb_ + pcap_;
    for (int i = 1; i < m; ++i)
    {
    //     //double tempfb = fbrine_local_[i];
        fbrine_local_[i] = fbrine_local_[i-1] - (myLayer_->get_rhob() * sb_local_[i] + myLayer_->get_rhoc()*sc_local_[i]) * G * dz - sc_local_[i] * (pcap_local_[i] - pcap_local_[i-1]);
        // fbrine_local_[i] = 0.0;
        pbrine_local_[i] = pb_ + fbrine_local_[i];


    //     //fbrine_local_[i] = fbrine_local_[i-1] + sc_local_[i] * myLayer_->get_rhoc() * G * dz + sb_local_[i] * myLayer_->get_rhob() * G * dz;
    //     // if (sc_local_[i] < 1.0E-40)
    //     // {
    //     //     tempi++;
    //     //     fbrine_local_[i] = deltarho * G * tempi * dz;
    //     // }
    //     // else {
    //     //     fbrine_local_[i] = pcap_local_[tempi] + deltarho * G * i * dz - pcap_local_[i];
    //     //     //fbrine_local_[i] = deltarho * G * (tempi - 1.0) * dz - pcap_local_[i] + pd_local_[i];
    //     // }

    //     // if (sc_local_[i] < 1.0E-2)
    //     // {
    //     //     fbrine_local_[i] = -myLayer_->get_rhob() * G * dz * tempi;
    //     //     tempi++;
    //     // }
    //     // else fbrine_local_[i] = -myLayer_->get_rhob() * G * dz * tempi - myLayer_->get_rhoc() * G * (i - tempi) * dz;

    //     //fbrine_local_[i] = 0.0;
    //     //fbrine_local_[i] = - pcap_local_[i] + deltarho * G * i * dz + pcap_local_[0];
    //     //fbrine_local_[i] = pbrine_local_[i] - ( pb_ - i * myLayer_->get_rhob()*dz*G );
    //     //dfbdt_local_[i] = (fbrine_local_[i] - tempfb)/DT;
    }
    // for (int i = m/2 - 1; i >= 0; --i)
    // {
    //     fbrine_local_[i] = fbrine_local_[i+1] + (myLayer_->get_rhob() * sb_local_[i] + myLayer_->get_rhoc()*sc_local_[i]) * G * dz + sc_local_[i] * (pcap_local_[i+1] - pcap_local_[i]);
    // }
    // for (int i = m/2 + 1; i < m; ++i)
    // {
    //     fbrine_local_[i] = fbrine_local_[i-1] - (myLayer_->get_rhob() * sb_local_[i] + myLayer_->get_rhoc()*sc_local_[i]) * G * dz - sc_local_[i] * (pcap_local_[i] - pcap_local_[i-1]);
    // }

    if (myLayer_->get_i() % (365*100) == 0 && sc_ > 1.0E-40)
    {
        ofstream myfile6;
        myfile6.open ("fine_totalflux.csv",ios::app);
        if (myfile6.is_open()){
            myfile6 << myLayer_->get_i() << ","<< ID_;
            for (int i = 0; i < m+1; ++i)
            {
                myfile6 << ","<< totalflux[i];  
            }         
            myfile6 << endl;
        }
        else cout << "failed to open file6" << endl;
        myfile6.close();

        ofstream myfile7;
        myfile7.open ("fbrine.csv",ios::app);
        if (myfile7.is_open()){
            myfile7 << myLayer_->get_i() << ","<< ID_;
            for (int i = 0; i < m; ++i)
            {
                myfile7 << ","<< fbrine_local_[i];  
            }         
            myfile7 << endl;
        }
        else cout << "failed to open file7" << endl;
        myfile7.close();

        ofstream myfile8;
        myfile8.open ("massdivergence.csv",ios::app);
        if (myfile8.is_open()){
            myfile8 << myLayer_->get_i() << ","<< ID_;
            myfile8 <<","<< verticalmassbalance << "," << horizontalmassbalance << endl;        
        }
        else cout << "failed to open file8" << endl;
        myfile8.close();

        // ofstream myfile9;
        // myfile9.open ("verticalfluxb.csv",ios::app);
        // if (myfile9.is_open()){
        //     myfile9 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m-1; ++i)
        //     {
        //         myfile9 <<","<< 1.0E13*fluxb[i]*100;  
        //     }         
        //     myfile9 << endl;
        // }
        // else cout << "failed to open file9" << endl;
        // myfile9.close();

        // ofstream myfile10;
        // myfile10.open("verticalfluxc.csv",ios::app);
        // if (myfile10.is_open()){
        //     myfile10 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m-1; ++i)
        //     {
        //         myfile10 <<","<< 1.0E13*fluxc[i]*100;  
        //     }         
        //     myfile10 << endl;
        // }
        // else cout << "failed to open file10" << endl;
        // myfile10.close();

        ofstream myfile11;
        myfile11.open("fine_pbrine.csv",ios::app);
        if (myfile11.is_open()){
            myfile11 << myLayer_->get_i() << ","<< ID_;
            for (int i = 0; i < m; ++i)
            {
                myfile11 <<","<< pbrine_local_[i];  
            }         
            myfile11 << endl;
        }
        else cout << "failed to open file11" << endl;
        myfile11.close();

        // ofstream myfile12;
        // myfile12.open("fine_pco2.csv",ios::app);
        // if (myfile12.is_open()){
        //     myfile12 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m; ++i)
        //     {
        //         myfile12 <<","<< pco2_local_[i];  
        //     }         
        //     myfile12 << endl;
        // }
        // else cout << "failed to open file12" << endl;
        // myfile12.close();

        ofstream myfile13;
        myfile13.open("fine_sco2.csv",ios::app);
        if (myfile13.is_open()){
            myfile13 << myLayer_->get_i() << ","<< ID_;
            for (int i = 0; i < m; ++i)
            {
                myfile13 <<","<< sc_local_[i];  
            }         
            myfile13 << endl;
        }
        else cout << "failed to open file13" << endl;
        myfile13.close();

        // ofstream myfile14;
        // myfile14.open("fine_pcap.csv",ios::app);
        // if (myfile14.is_open()){
        //     myfile14 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m; ++i)
        //     {
        //         myfile14 <<","<< pcap_local_[i];  
        //     }         
        //     myfile14 << endl;
        // }
        // else cout << "failed to open file14" << endl;
        // myfile14.close();

        // ofstream myfile15;
        // myfile15.open("fine_interMobilityb.csv",ios::app);
        // if (myfile15.is_open()){
        //     myfile15 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m-1; ++i)
        //     {
        //         myfile15 <<","<< interMobilityb_local_[i];  
        //     }         
        //     myfile15 << endl;
        // }
        // else cout << "failed to open file15" << endl;
        // myfile15.close();

        // ofstream myfile16;
        // myfile16.open("fine_interMobilityc.csv",ios::app);
        // if (myfile16.is_open()){
        //     myfile16 << myLayer_->get_i() << ","<< ID_;
        //     for (int i = 0; i < m-1; ++i)
        //     {
        //         myfile16 <<","<< interMobilityc_local_[i];  
        //     }         
        //     myfile16 << endl;
        // }
        // else cout << "failed to open file16" << endl;
        // myfile16.close();


    }
  
    
            
}

double Cell::get_co2Mass(){
    return area_ * H_ * sc_ * phi_;
}

void Cell::output_nodeid()
{
    vector<Node*>::iterator cellnodes_iter = cellnodes_.begin();
    vector<Node*>::iterator cellnodes_end = cellnodes_.end();
    while (cellnodes_iter!=cellnodes_end) {
        //        cout << (*cellnodes_iter)->Get_ID()<< endl;
        cellnodes_iter++;
    }    
}

void Cell::output_interfaceid()
{
    vector<int>::iterator interfaceid_iter = interfaceid_.begin();
    vector<int>::iterator interfaceid_end = interfaceid_.end();
    while (interfaceid_iter!=interfaceid_end) {
        cout << (*interfaceid_iter) << endl;
        interfaceid_iter++;
    }  
}


