#include "Interface.h"
#include "Layer.h"
#include "Const.h"
#include "math.h"
#include "fstream"
using namespace std;

void Interface::create_( vector<int> nodeid, vector<int>cellid, double fluxb, double fluxc, double pb, double sc, int interfacetype, int id, Layer* myLayer){
    nodeid_ = nodeid;
    cellid_ = cellid;
    fluxb_ = fluxb;
    fluxc_ = fluxc;
    pb_ = pb;
    sc_ = sc; 
    ID_ = id;
    interfacetype_ = interfacetype;
    myLayer_ = myLayer;

    for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
    {
        pd_local_.push_back(0.0E4);
        mobilityb_local_.push_back(0.0);
        mobilityc_local_.push_back(0.0);
        fluxb_local_.push_back(0.0);
        fluxc_local_.push_back(0.0);
    }

}

void Interface::create_(){
    pb_ = 0.0;
    sc_ = 0.0;
    interfacetype_ = 0;
}

void Interface::initialize(){
    double xC, xD, yC, yD, zC1, zD1, zC2, zD2, zC3, zD3;
    double xA = interfacenodes_[0]->get_x(); double xB = interfacenodes_[1]->get_x(); 
    double yA = interfacenodes_[0]->get_y(); double yB = interfacenodes_[1]->get_y();
    double z1A = interfacenodes_[0]->get_z1(); double z1B = interfacenodes_[1]->get_z1();
    double z2A = interfacenodes_[0]->get_z2(); double z2B = interfacenodes_[1]->get_z2();
    lengthAB_ = sqrt(pow((xA-xB),2)+pow((yA-yB), 2)); // length of the interface
    
    
    
    // Initialize the interface center
    centerx_ = (xA + xB)/2;
    centery_ = (yA + yB)/2;
    centerz1_ = (z1A + z1B)/2;
    centerz2_ = (z2A + z2B)/2;
    
    H_ = centerz2_ - centerz1_;
    
    
    if (interfacetype_ == 0) {
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            // Initialize local permeability, using harmonic average of the two cells.
            permeability_local_.push_back(2/(1/(*(interfacecells_[0]->get_permeability_local() + i))+1/(*(interfacecells_[1]->get_permeability_local() + i))));
        }
        
        xC = interfacecells_[0]->get_centerx(); yC = interfacecells_[0]->get_centery();
        xD = interfacecells_[1]->get_centerx(); yD = interfacecells_[1]->get_centery();
        
        zC1 = interfacecells_[0]->get_centerz2(); zD1 = interfacecells_[1]->get_centerz2();
        zC3 = interfacecells_[0]->get_centerz1(); zD3 = interfacecells_[1]->get_centerz1();
        zC2 = (zC1 + zC3)/2;; zD2 = (zD1 + zD3)/2;
        
        centerdistance_ = pow( pow((xC - xD),2) + pow((yC - yD),2) + pow((zC2 - zD2),2),0.5);
        
        
    }
    else if (interfacetype_ == -1){
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            permeability_local_.push_back(*(interfacecells_[0]->get_permeability_local() + i));
        }

        xC = (xA + xB)/2; yC = (yA + yB)/2; 
        xD = interfacecells_[0]->get_centerx(); yD = interfacecells_[0]->get_centery();
        
        zC1 = (z2A + z2B)/2; zD1 = interfacecells_[0]->get_centerz2();
        zC3 = (z1A + z1B)/2; zD3 = interfacecells_[0]->get_centerz1();
        zC2 = (zC1 + zC3)/2; zD2 = (zD1 + zD3)/2;
        
        centerdistance_ = pow( pow((xC - xD),2) + pow((yC - yD),2) + pow((zC2 - zD2),2),0.5);

        // Give pressure boundary
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            pb_local_.push_back( pb_ - i * H_/(myLayer_->get_verRefineNum())*myLayer_->get_rhob() * G);
            pcap_local_.push_back( pd_local_[i] );
        }
        
    }

    else if (interfacetype_ == 1){     
        xC = (xA + xB)/2; yC = (yA + yB)/2; 
        xD = interfacecells_[0]->get_centerx(); yD = interfacecells_[0]->get_centery();
        
        zC1 = (z2A + z2B)/2; zD1 = interfacecells_[0]->get_centerz2();
        zC3 = (z1A + z1B)/2; zD3 = interfacecells_[0]->get_centerz1();
        zC2 = (zC1 + zC3)/2; zD2 = (zD1 + zD3)/2;
        
        centerdistance_ = pow( pow((xC - xD),2) + pow((yC - yD),2) + pow((zC2 - zD2),2),0.5);

        // Give no flow bounddary
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i){
            fluxb_local_[i]=fluxb_;
            fluxc_local_[i]=fluxc_;
        }
    }
    else  {
        cout << "ERROR: interfacetype is out of scope" << endl;
    }
    
    // CD*n*|AB|/|CD|^2   vectorprojection
    double temp = -((xD-xC)*(yA-yB)+(yD-yC)*(xA-xB))/(centerdistance_*pow( pow((xC - xD),2) + pow((yC - yD),2),0.5));
    if (temp>0){
        vectorprojection_ = temp;
    }
    else vectorprojection_ = -temp;  
    
}

 

void Interface::set_interfacenodes(){
    Node* temp;
    vector<int>::iterator nodeid_iter = nodeid_.begin();
    vector<int>::iterator nodeid_end = nodeid_.end();
    while (nodeid_iter != nodeid_end) {
        temp = (*myLayer_).findNodeByID(*nodeid_iter);
        interfacenodes_.push_back(temp);
        nodeid_iter++;}
}

void Interface::set_interfacecells(){
    vector<int>::iterator cellid_iter = cellid_.begin();
    vector<int>::iterator cellid_end = cellid_.end();
    while (cellid_iter != cellid_end){
        interfacecells_.push_back((*myLayer_).findCellByID(*cellid_iter));
        cellid_iter++;}
    //check interfacecells
//    vector<Cell*>::iterator interfacecells_iter = interfacecells_.begin();
//    vector<Cell*>::iterator interfacecells_end = interfacecells_.end();
//    while (interfacecells_iter!=interfacecells_end) {
//        cout << "interfacecellid" << (*interfacecells_iter)->Get_ID() << endl;
//        interfacecells_iter++;
//    }
}

Cell* Interface::get_cell(int cellid){
    vector<Cell*>::iterator interfacecells_iter = interfacecells_.begin();
    vector<Cell*>::iterator interfacecells_end = interfacecells_.end();
    while (interfacecells_iter != interfacecells_end) {
        int temp = (*interfacecells_iter)->get_ID();
//        cout <<"interfacecellid" <<temp <<endl;
        if ( temp!= cellid) {
            return *interfacecells_iter; 
        }
        interfacecells_iter++;
    }
    return NULL;
}

void Interface::set_flux(){
    if (interfacetype_ == 0) {
        //fluxc is from cell[0] to cell[1]
        
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {      
            //tempc and tempb are for checking the flux direction, then we can define interface mobility based on that.
            double tempc = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine()
                + (*(interfacecells_[1]->get_pcap_local() + i)) - (*(interfacecells_[0]->get_pcap_local() + i))
                + myLayer_->get_rhoc()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1()) 
                + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            double tempb = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine()
                    + myLayer_->get_rhob()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1())
                    + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
        
            mobilityb_local_[i] = (*(interfacecells_[0]->get_mobilityb_local() + i))*(tempb>=0) + (*(interfacecells_[1]->get_mobilityb_local() + i))*(tempb<0);
            mobilityc_local_[i] = (*(interfacecells_[0]->get_mobilityc_local() + i))*(tempc>=0) + (*(interfacecells_[1]->get_mobilityc_local() + i))*(tempc<0);        
            
            double tempc1 = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine() 
                + (*(interfacecells_[1]->get_pcap_local() + i)) - (*(interfacecells_[0]->get_pcap_local() + i))
                + myLayer_->get_rhoc()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1()) 
                + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            double tempb1 = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine() 
                    + myLayer_->get_rhob()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1())
                    + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;


            fluxc_local_[i] = permeability_local_[i] * mobilityc_local_[i] * tempc1;
            fluxb_local_[i] = permeability_local_[i] * mobilityb_local_[i] * tempb1;
            /*
            cout << tempc << endl;
            cout << tempb << endl;
            cout << "fluxc:" << fluxc_local_[i] << "fluxb:" << fluxb_local_[i] << endl;
            */
            
        }

        // if (myLayer_->get_i() % (365) == 0 && ID_ > (myLayer_->get_Nx()-1)*myLayer_->get_Ny())
        // {
        //     ofstream myfile1;
        //     myfile1.open ("horizontalfluxb.csv",ios::app);
        //     if (myfile1.is_open()){
        //         myfile1 << myLayer_->get_i() << ","<< ID_;
        //         for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        //         {
        //             if (interfacecells_[0]->get_centerx() < interfacecells_[1]->get_centerx() )
        //             {
        //                 myfile1 << ","<< 1.0E13*fluxb_local_[i]/100.; 
        //             }
        //             else myfile1 << ","<< - 1.0E13*fluxb_local_[i]/100.;
        //         }         
        //         myfile1 << endl;
        //     }
        //     else cout << "failed to open myfile1 in interface" << endl;
        //     myfile1.close();

        //     ofstream myfile2;
        //     myfile2.open ("horizontalfluxc.csv",ios::app);
        //     if (myfile2.is_open()){
        //         myfile2 << myLayer_->get_i() << ","<< ID_;
        //         for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        //         {
        //             if (interfacecells_[0]->get_centerx() < interfacecells_[1]->get_centerx() )
        //             {
        //                 myfile2 << ","<< 1.0E13*fluxc_local_[i]/100.; 
        //             }
        //             else myfile2 << ","<< - 1.0E13*fluxc_local_[i]/100.;
        //         }         
        //         myfile2 << endl;
        //     }
        //     else cout << "failed to open myfile2 in interface" << endl;
        //     myfile2.close();

        //     ofstream myfile3;
        //     myfile3.open ("mobilityb_local.csv",ios::app);
        //     if (myfile3.is_open()){
        //         myfile3 << myLayer_->get_i() << ","<< ID_;
        //         for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        //         {

        //             myfile3 << ","<< mobilityb_local_[i]; 

        //         }         
        //         myfile3 << endl;
        //     }
        //     else cout << "failed to open myfile3 in interface" << endl;
        //     myfile3.close();

        //     ofstream myfile4;
        //     myfile4.open ("mobilityc_local.csv",ios::app);
        //     if (myfile4.is_open()){
        //         myfile4 << myLayer_->get_i() << ","<< ID_;
        //         for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        //         {

        //             myfile4 << ","<< mobilityc_local_[i]; 

        //         }         
        //         myfile4 << endl;
        //     }
        //     else cout << "failed to open myfile4 in interface" << endl;
        //     myfile4.close();
        // }
    }
    else if (interfacetype_ == -1) {
        // fluxc is from cell center to the boundary
        
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            mobilityb_local_[i] = *(interfacecells_[0]->get_mobilityb_local() + i);
            mobilityc_local_[i] = *(interfacecells_[0]->get_mobilityc_local() + i);
        
            fluxc_local_[i] = - permeability_local_[i] * mobilityc_local_[i] * (pb_ - interfacecells_[0]->get_pbrine() 
                + 0*pcap_local_[i] - 0*(*(interfacecells_[0]->get_pcap_local() + i))
                + myLayer_->get_rhob()*G*(centerz1_ - interfacecells_[0]->get_centerz1())
                + 0.0 - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            /*
            cout << pb_ << endl;
            cout << interfacecells_[0]->get_pbrine() << endl;
            cout << pcap_local_[i] << endl;
            cout << (*(interfacecells_[0]->get_pcap_local() + i)) << endl;
            cout << (*(interfacecells_[0]->get_fbrine_local() + i)) << endl;
            */
            
            fluxb_local_[i] = - permeability_local_[i] * mobilityb_local_[i] * (pb_ - interfacecells_[0]->get_pbrine() 
                + myLayer_->get_rhob()*G*(centerz1_ - interfacecells_[0]->get_centerz1())
                + 0.0 - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            /*
            cout << pb_ << endl;
            cout << interfacecells_[0]->get_pbrine() << endl;
            cout << centerz1_ << endl;
            cout << interfacecells_[0]->get_centerz1() << endl;
            cout << (*(interfacecells_[0]->get_fbrine_local() + i)) << endl;
            cout << pb_ - interfacecells_[0]->get_pbrine() 
                + myLayer_->get_rhob()*G*(centerz1_ - interfacecells_[0]->get_centerz1())
                + 0.0 - (*(interfacecells_[0]->get_fbrine_local() + i)) << endl;

            cout << "fluxc:" << fluxc_local_[i] << "fluxb:" << fluxb_local_[i] << endl;
            */
        }
        
    }
    else if (interfacetype_ == 1)
    {
       for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
       {
            fluxb_local_[i] = fluxb_;
            fluxc_local_[i] = fluxc_;
       }
    }
    else {
        cout << "ERROR: Interfacetype is out of scope" << endl;
    }
    
}

void Interface::set_pbrinecoeffiandRHS(){
    pbrinecoeffi_.clear();
    pbrineRHS_.clear();
    //for interior interface(interfacetype = 0), there are two coefficients, A[i][i] and A[i][j]
    if (interfacetype_ == 0) {
        
        // tempmobilityb and tempmobilityc are for store temp integrated mobilities 
        double tempmobilityb = 0.0;
        double tempmobilityc = 0.0;
        // vertically-integrated non-equilibrium fluxes
        double tempdeviationb = 0.0;
        double temppcap = 0.0;
        double dz = H_/myLayer_->get_verRefineNum();
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {      
            //tempc and tempb are for checking the flux direction, then we can define interface mobility based on that.
            double tempc = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine() 
                + (*(interfacecells_[1]->get_pcap_local() + i)) - (*(interfacecells_[0]->get_pcap_local() + i))
                + myLayer_->get_rhoc()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1()) 
                + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            double tempb = - (interfacecells_[1]->get_pbrine() - interfacecells_[0]->get_pbrine() 
                    + myLayer_->get_rhob()*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1())
                    + (*(interfacecells_[1]->get_fbrine_local() + i)) - (*(interfacecells_[0]->get_fbrine_local() + i))
                )*vectorprojection_;
            
            // fine scale mobilities upstream weighted from cells
            mobilityb_local_[i] = (*(interfacecells_[0]->get_mobilityb_local() + i))*(tempb>=0) + (*(interfacecells_[1]->get_mobilityb_local() + i))*(tempb<0);
            mobilityc_local_[i] = (*(interfacecells_[0]->get_mobilityc_local() + i))*(tempc>=0) + (*(interfacecells_[1]->get_mobilityc_local() + i))*(tempc<0); 

            tempmobilityb += permeability_local_[i] * mobilityb_local_[i] * dz;
            tempmobilityc += permeability_local_[i] * mobilityc_local_[i] * dz;

            tempdeviationb += permeability_local_[i] * (mobilityb_local_[i] + mobilityc_local_[i]) * (*(interfacecells_[1]->get_fbrine_local() + i) - (*(interfacecells_[0]->get_fbrine_local()+i))) * dz; 
            temppcap += permeability_local_[i] * mobilityc_local_[i] * (*(interfacecells_[1]->get_pcap_local() + i) - (*(interfacecells_[0]->get_pcap_local() + i))) * dz;
        }
        
        kmobilityb_ = tempmobilityb;
        kmobilityc_ = tempmobilityc;
        pbrinecoeffi_.push_back((kmobilityb_ + kmobilityc_) * vectorprojection_);
        pbrinecoeffi_.push_back(-(kmobilityb_ + kmobilityc_) * vectorprojection_); // first element is A[i][i] and A[j][j], second element is A[i][j] and A[j][i]
        
        // temp is the pbrineRHS
        double temp = ( - (kmobilityb_ * myLayer_->get_rhob() + kmobilityc_ * myLayer_->get_rhoc())*G*(interfacecells_[1]->get_centerz1() - interfacecells_[0]->get_centerz1()) 
                            + tempdeviationb + temppcap
                        )*vectorprojection_;
        pbrineRHS_.push_back(temp);
        pbrineRHS_.push_back(-temp); // first element is for i, second is for j

        // if (myLayer_->get_i() % (365) == 0 && ID_ > (myLayer_->get_Nx()-1)*myLayer_->get_Ny())
        // {

        //     ofstream myfile3;
        //     myfile3.open ("tempmobilityb.csv",ios::app);
        //     if (myfile3.is_open()){
        //         myfile3 << myLayer_->get_i() << ","<< ID_;
        //         myfile3 << ","<< tempmobilityb*1.0E14/2 << endl;        
        //     }
        //     else cout << "failed to open myfile3 in interface" << endl;
        //     myfile3.close();

        //     ofstream myfile4;
        //     myfile4.open ("tempmobilityc.csv",ios::app);
        //     if (myfile4.is_open()){
        //         myfile4 << myLayer_->get_i() << ","<< ID_;
        //         myfile4 << ","<< tempmobilityc*1.0E14/2 << endl;        
        //     }
        //     else cout << "failed to open myfile4 in interface" << endl;
        //     myfile4.close();
        // }




    }
    //for pressure boundary(interfacetype = -1), there is only one coefficient, A[i][i]
    else if (interfacetype_ == -1){
        // tempmobilityb and tempmobilityc are for store temp integrated mobilities
        double tempmobilityb = 0.0;
        double tempmobilityc = 0.0;
        // vertically-integrated non-equilibrium fluxes
        double tempdeviationb = 0.0;
        double temppcap = 0.0;
        double dz = H_/myLayer_->get_verRefineNum();
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            mobilityb_local_[i] = (*(interfacecells_[0]->get_mobilityb_local() + i));
            mobilityc_local_[i] = (*(interfacecells_[0]->get_mobilityc_local() + i));

            tempmobilityb += permeability_local_[i] * mobilityb_local_[i] * dz;
            tempmobilityc += permeability_local_[i] * mobilityc_local_[i] * dz;

            tempdeviationb += permeability_local_[i] * (mobilityb_local_[i] + mobilityc_local_[i]) * ( 0.0 - *(interfacecells_[0]->get_fbrine_local() + i) );
            temppcap += permeability_local_[i] * mobilityc_local_[i] * (pcap_local_[i]- (*(interfacecells_[0]->get_pcap_local() + i)) );
        }

        kmobilityb_ = tempmobilityb;
        kmobilityc_ = tempmobilityc;

        pbrinecoeffi_.push_back((kmobilityb_ + kmobilityc_)*vectorprojection_); //this coeffi is for A[j][j]
          
        // temp is the pbrineRHS
        double temp = ((kmobilityb_ + kmobilityc_)*pb_ 
                        - (kmobilityb_ * myLayer_->get_rhob() + kmobilityc_ * myLayer_->get_rhoc())*G*(centerz1_ - interfacecells_[0]->get_centerz1()) 
                        + tempdeviationb + temppcap
                        )*vectorprojection_;
        pbrineRHS_.push_back(temp);

    }
    //for flux boundary(interfacetype = 1), there is no matrix coeffient, all the contribution goes into the RHS
    else if (interfacetype_ == 1){
        double temp = 0.0;
        for (int i = 0; i < myLayer_->get_verRefineNum(); ++i)
        {
            temp += (fluxb_local_[i]+fluxc_local_[i]);
        }
        pbrineRHS_.push_back(temp * vectorprojection_ * centerdistance_);
    }

    else {
        cout << "ERROR: interfacetype is out of scope" << endl;
    }
}


void Interface::output_nodeid(){
    vector<int>::iterator nodeid_iter = nodeid_.begin();
    vector<int>::iterator nodeid_end = nodeid_.end();
    while (nodeid_iter != nodeid_end) {
        cout << (*nodeid_iter) << endl;
        nodeid_iter++;
    }
}

void Interface::output_cellid(){
    vector<Cell*>::iterator interfacecells_iter = interfacecells_.begin();
    vector<Cell*>::iterator interfacecells_end = interfacecells_.end();
    while (interfacecells_iter != interfacecells_end) {
        cout << (*interfacecells_iter)->get_ID() << endl;
        interfacecells_iter++;
    }
}
