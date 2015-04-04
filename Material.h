//
//  function.h
//  FractionalFlow_VR
//
//  Created by Bo Guo on 7/31/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef FractionalFlow_VR_function_h
#define FractionalFlow_VR_function_h
#include <fstream>
#include <math.h>
#include <cmath>

double pcap ( double Sc, double Pd, double lambda );//the function of capillary pressure
// double dpcapdsc ( double Sc, double Pd, double lambda ); //the function of dpcap/ds
double mobilityco2 ( double Sc, double lambda);  //kapa function for CO2
double mobilitybrine ( double Sb, double lambda);  //kapa function for Brine
// double dlambdacdsc ( double Sc, double lambda); // dlambdac/dsc, evaluate by upstream
// double dlambdabdsc ( double Sc, double lambda); // dlambdab/dsc, evaluate by upstream


#endif
