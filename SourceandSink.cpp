//
//  SourceandSink.cpp
//  Simplegeometry
//
//  Created by Bo Guo on 8/14/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "SourceandSink.h"
#include "Layer.h"
void SourceandSink::create_(int ID, double qc, double qb, Layer* myLayer){
    ID_ = ID;
    qc_ = qc;
    qb_ = qb;
    myLayer_ = myLayer;
}

void SourceandSink::create_(){
    ID_ = 0;
    qc_ = 0;
    qb_ = 0;
    myLayer_ = 0;
}