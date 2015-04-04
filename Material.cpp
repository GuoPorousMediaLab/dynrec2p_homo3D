//
//  function.cpp
//  FractionalFlow_VR
//
//  Created by Bo Guo on 7/31/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Material.h"
using namespace std;

double pcap ( double Sc, double Pd, double lambda ) //the function of capillary pressure
{
    //    Brooks Corey
    double ans, Sbr = 0.0, S1;
    S1 = (1 - Sc - Sbr)/(1.0 - Sbr);
    // if (Pd*pow(S1, - 1.0/lambda) <= 1.0E6)
    // {
    //     ans = Pd*pow(S1, - 1.0/lambda);
    // }
    // else ans = 1.0E6;
    ans = Pd*pow(S1, - 1.0/lambda);
    return 0.0;

    // V-G, zhou
    // double ans, Sbr = 0.3, S1, m = 0.46;
    // S1 = (1 - Sc - Sbr)/(1.0 - Sbr);
    // ans = Pd*pow(pow(S1, -1.0/m)-1, 1-m);
    //return ans;
    
}

// double dpcapdsc ( double Sc, double Pd, double lambda) //the function of dcap/ds
// {
//     //   Brooks Corey
//     double ans, Sbr = 0.2, S1;
//     S1 = (1 - Sc - Sbr)/(1.0 - Sbr);
//     ans = Pd/lambda/(1.0 - Sbr)*pow(S1, - 1/lambda - 1);
//     return ans;
// }


double mobilityco2 ( double Sc, double lambda)  //kapa function for CO2
{
    
    //    Brooks Corey
    double ans, Se, krc, Sbr = 0.0, Scr = 0.0, miuc = 4.25E-5;
    Se = (1 - Sc - Sbr)/(1-Sbr-Scr);
    krc = pow((1-Se),2)*(1-pow(Se,(2+lambda)/lambda));
    ans = krc/miuc;
    //ans = (1-Se)/miuc;
    //ans = pow((1-Se),2)/miuc;
    //ans = pow(1 - Se,3)/miuc;
    if (Sc < 0.0 )
    {
       ans = 0.0;
    }
    else if (Sc > 1.0)
    {
        ans = 1.0/miuc;
    }
    return ans;

    // V-G, zhou
    // double ans, S1, krc, Sbr = 0.3, Scr = 0.0, miuc = 4.25E-5, m = 0.46;
    // S1 = (1 - Sc - Sbr)/(1-Sbr-Scr);
    // krc = pow((1-S1),0.5)*pow(1-pow(S1,1.0/m),2*m);
    // ans = krc/miuc;
    // return ans;
    
}

double mobilitybrine ( double Sb, double lambda)  //kapa function for Brine
{
    //    Brooks Corey
    double ans, Se, krb, Sbr = 0.0, Scr = 0.0, miub = 3.0E-4;
    Se = (Sb - Sbr)/(1-Sbr-Scr);
    krb = pow(Se, (2.0+3.0*lambda)/lambda);
    ans = krb/miub;
    //ans = Se/miub;
    //ans = pow(Se,2)/miub;
    //ans = pow(Se,2)/miub;
    if (Sb < 0.0)
    {
       ans = 0.0;
    }
    else if (Sb > 1.0)
    {
        ans = 1.0/miub;
    }
    return ans;

    // V-G, zhou
    // double ans, S1, krb, Sbr = 0.3, Scr = 0.0, miub = 3.0E-4, m = 0.46;
    // S1 = (Sb - Sbr)/(1-Sbr-Scr);
    // krb = pow(S1, 0.5)*pow(1 - pow(1 - pow(S1, 1.0/m), m), 2); 
    // ans = krb/miub;
    // return ans;
}

// double dlambdacdsc ( double Sc, double lambda) // dlambdac/dsc
// {
//     // Brooks Corey
//     double ans, Se, Sbr = 0.2, Scr = 0.0, miuc = 4.25E-5;
//     Se = (1 - Sc - Sbr)/(1-Sbr-Scr);
//     ans = (1 - Se)/(miuc*(1 - Sbr - Scr))*(2*(1 - pow(Se, (2+lambda)/lambda))  + (1-Se)*(2+lambda)/lambda*pow(Se, 2/lambda));
//     return ans;
// }

// double dlambdabdsc ( double Sc, double lambda) // dlambdab/dsc
// {
//     // Brooks Corey
//     double ans, Se, Sbr = 0.2, Scr = 0.0, miub = 3.0E-4;
//     Se = (1 - Sc - Sbr)/(1-Sbr-Scr);
//     ans = -1/(miub*(1 - Sbr - Scr))*(2+3*lambda)/lambda*pow(Se, 2*(1+lambda)/lambda);
//     return ans;
// }
