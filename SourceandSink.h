//
//  SourceandSink.h
//  Simplegeometry
//
//  Created by Bo Guo on 8/14/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Simplegeometry_SourceandSink_h
#define Simplegeometry_SourceandSink_h
class Layer;

class SourceandSink
{
public:
    SourceandSink(int ID, double qc, double qb, Layer* myLayer){create_(ID, qc, qb, myLayer);}
    SourceandSink(){create_();}
    int get_ID(){return ID_;}
    double get_co2source(){return qc_;}
    double get_brinesource(){return qb_;}
private:
    Layer* myLayer_;
    double qc_, qb_;
    int ID_;
    void create_(int, double, double, Layer*);
    void create_();
};


#endif
