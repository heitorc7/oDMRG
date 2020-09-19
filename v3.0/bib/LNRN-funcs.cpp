#include <iostream>
#include <fstream>
#include "math.h"
#include <complex.h>
#include <vector>

using namespace itensor;

void MakeIVEC_LNRN(itensor::MPS &Ivec, int L)    
// Overwrites the MPS entry with a vectorized Identity with dimension L ready to be used in calculations
{

    auto A =Ivec.A(1);
    auto indi = A.inds();
    auto b2 = Index(2);

    //First and last sites are set by hand
    auto B = ITensor(indi[0],b2);
    B.set(indi[1],b2(1),1);
    /*
    B.set(indi[0](1),b2(2),0);
    B.set(indi[0](2),b2(1),0);
    B.set(indi[0](2),b2(2),1);

    Ivec.setA(1,B);
    for(int i=1;i<L;i++){
        auto A =Ivec.A(2*i);
        auto indi = A.inds();

        auto A0 = Ivec.A(2*i-1);
        auto b1 = (A0.inds())[(i==1)? 1:2];
        auto b2 = Index(1);

        auto B = ITensor(b1,indi[1],b2);

        B.set(b1(1),indi[1](1),b2(1),1);
        B.set(b1(2),indi[1](1),b2(1),0);
        B.set(b1(1),indi[1](2),b2(1),0);
        B.set(b1(2),indi[1](2),b2(1),1);

        Ivec.setA(2*i,B);


        A =Ivec.A(2*i+1);
        indi = A.inds();

        b1 = b2;
        b2 = Index(2);

        B = ITensor(b1,indi[1],b2);

        B.set(b1(1),indi[1](1),b2(1),1);
        B.set(b1(1),indi[1](1),b2(2),0);
        B.set(b1(1),indi[1](2),b2(1),0);
        B.set(b1(1),indi[1](2),b2(2),1);

        Ivec.setA(2*i+1,B);
    }
    A =Ivec.A(2*L);
    indi = A.inds();

    auto A0 = Ivec.A(2*L-1);
    auto b1 = (A0.inds())[2];

    B = ITensor(b1,indi[1]);
    B.set(b1(1),indi[1](1),1);
    B.set(b1(2),indi[1](1),0);
    B.set(b1(1),indi[1](2),0);
    B.set(b1(2),indi[1](2),1);

    Ivec.setA(2*L,B);
    */
    int a = 1;
} //MakeIVEC

AutoMPO LiouvConstruct_LNRN(double gamma, double f1, double fL, double h, double Delta, SiteSet sites)       
// Constructs a Liouvillian given the relevant parameters and a SiteSet object
// to be used in pair with LiouvDConstruct_LNRN(...)
{                                                                                                       
    int L = length(sites)/2;                // relevant for looping through
    AutoMPO ampo = AutoMPO(sites);

    /* Construct D_1_m */
    ampo +=  gamma*f1,     "S-", 1, "S-", 2;
    ampo += -gamma*f1*0.5, "S+", 1, "S-", 1;
    ampo += -gamma*f1*0.5, "S+", 2, "S-", 2;

    /* Construct D_1_p */
    ampo +=  gamma*(1.0-f1),     "S+", 1, "S+", 2;
    ampo += -gamma*(1.0-f1)*0.5, "S-", 1, "S+", 1;
    ampo += -gamma*(1.0-f1)*0.5, "S-", 2, "S+", 2;

    /* Construct D_L_m */
    ampo +=  gamma*fL,     "S-", 2*L-1, "S-", 2*L;
    ampo += -gamma*fL*0.5, "S+", 2*L-1, "S-", 2*L-1;
    ampo += -gamma*fL*0.5, "S+", 2*L, "S-", 2*L;

    /* Construct D_L_p */
    ampo +=  gamma*(1.0-fL),     "S+", 2*L-1, "S+", 2*L;
    ampo += -gamma*(1.0-fL)*0.5, "S-", 2*L-1, "S+", 2*L-1;
    ampo += -gamma*(1.0-fL)*0.5, "S-", 2*L, "S+", 2*L;

    /* Construct local Hamiltonian */
    for(int j = 1; j <= L; ++j)
    {
        ampo += -2*1_i*h, "Sz", 2*j-1;
        ampo +=  2*1_i*h, "Sz", 2*j;
    }

    /* Construct nearest-neighbor interactions */
    for(int j = 1; j < L; ++j)
    {
        ampo += -4*1_i, "Sx", 2*j-1, "Sx", 2*j+1;
        ampo +=  4*1_i, "Sx", 2*j,   "Sx", 2*j+2;

        ampo += -4*1_i, "Sy", 2*j-1, "Sy", 2*j+1;
        ampo +=  4*1_i, "Sy", 2*j,   "Sy", 2*j+2;

        ampo += -4*1_i*Delta, "Sz", 2*j-1, "Sz", 2*j+1;
        ampo +=  4*1_i*Delta, "Sz", 2*j,   "Sz", 2*j+2;

    }

    return ampo;
} //LiouvConstruct



AutoMPO LiouvDConstruct_LNRN(double gamma, double f1, double fL, double h, double Delta, SiteSet sites)          
// Constructs a Liouvillian (dagger) given the relevant parameters and a site set
// To be used in pair with LiouvConstruct_LNRN
{                                                                                                           
    int L = length(sites)/2;                // relevant for looping through
    AutoMPO ampoD = AutoMPO(sites);
    /* Construct D_1_m */
    ampoD +=  gamma*f1,     "S+", 1, "S+", 2;
    ampoD += -gamma*f1*0.5, "S+", 1, "S-", 1;
    ampoD += -gamma*f1*0.5, "S+", 2, "S-", 2;

    /* Construct D_1_p */
    ampoD +=  gamma*(1.0-f1),     "S-", 1, "S-", 2;
    ampoD += -gamma*(1.0-f1)*0.5, "S-", 1, "S+", 1;
    ampoD += -gamma*(1.0-f1)*0.5, "S-", 2, "S+", 2;

    /* Construct D_L_m */
    ampoD +=  gamma*fL,     "S+", 2*L-1, "S+", 2*L;
    ampoD += -gamma*fL*0.5, "S+", 2*L-1, "S-", 2*L-1;
    ampoD += -gamma*fL*0.5, "S+", 2*L, "S-", 2*L;

    /* Construct D_L_p */
    ampoD +=  gamma*(1.0-fL),     "S-", 2*L-1, "S-", 2*L;
    ampoD += -gamma*(1.0-fL)*0.5, "S-", 2*L-1, "S+", 2*L-1;
    ampoD += -gamma*(1.0-fL)*0.5, "S-", 2*L, "S+", 2*L;

    /* Construct local Hamiltonian */
    for(int j = 1; j <= L; ++j)
    {
        ampoD +=  2*1_i*h, "Sz", 2*j-1;
        ampoD += -2*1_i*h, "Sz", 2*j;
    }

    /* Construct nearest-neighbor interactions */
    for(int j = 1; j < L; ++j)
    {
        ampoD +=  4*1_i, "Sx", 2*j-1, "Sx", 2*j+1;
        ampoD += -4*1_i, "Sx", 2*j,   "Sx", 2*j+2;

        ampoD +=  4*1_i, "Sy", 2*j-1, "Sy", 2*j+1;
        ampoD += -4*1_i, "Sy", 2*j,   "Sy", 2*j+2;

        ampoD +=  4*1_i*Delta, "Sz", 2*j-1, "Sz", 2*j+1;
        ampoD += -4*1_i*Delta, "Sz", 2*j,   "Sz", 2*j+2;

    }

    return ampoD;
} //LiouvDConstruct



// Experimental stuff.
//using Gate = itensor::BondGate<itensor::ITensor>;
//std::vector<Gate> TrotterConstruct_LNRN(double gamma, double f1, double fL, double h, double Delta, SiteSet sites, Real tstep)       //Constructs a TrotterLiouv illian given the relevant parameters and a SiteSet object
//{

//    int L = length(sites)/2;                // relevant for looping through
//    auto gates = std::vector<Gate>();

 //    /* Construct D_1_m */
 //    auto TrotterLiouv  = gamma*f1*sites.op("S-", 1)*sites.op("S-", 2);
 //    //TrotterLiouv += -gamma*f1*0.5*sites.op("S+", 1)*sites.op("S-", 1);
 //    TrotterLiouv  += -gamma*f1*0.5*sites.op("Sz", 1)*sites.op("Id", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    //TrotterLiouv  += -gamma*f1*0.5*sites.op("S+", 2)*sites.op("S-", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*sites.op("Id", 1)*sites.op("Sz", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);

	// /* Construct D_1_p */
 //    TrotterLiouv  +=  gamma*(1.0-f1)*sites.op("S+", 1)*sites.op("S+", 2);
 //    //TrotterLiouv  += -gamma*(1.0-f1)*0.5*sites.op("S-", 1)*sites.op("S+", 1);
 //    TrotterLiouv  += -gamma*(1.0-f1)*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    TrotterLiouv  -= -gamma*(1.0-f1)*0.5*sites.op("Sz", 1)*sites.op("Id", 2);
 //    //TrotterLiouv  += -gamma*(1.0-f1)*0.5*sites.op("S-", 2)*sites.op("S+", 2);
 //    TrotterLiouv  += -gamma*(1.0-f1)*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    TrotterLiouv  -= -gamma*(1.0-f1)*0.5*sites.op("Id", 1)*sites.op("Sz", 2);

 //    /* Local Hamiltonian for the first-site case */
 //    TrotterLiouv  += -2*1_i*h*sites.op("Sz", 1)*sites.op("Id",2);
 //    TrotterLiouv  +=  2*1_i*h*sites.op("Id",1)*sites.op("Sz", 2);

 //    /* Nearest-neighbor interactions for the first-site case */
 //    // TrotterLiouv  += -4*1_i*sites.op("Sx", 1)*sites.op("Sx", 3);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sx", 2)*sites.op("Sx", 4);

 //    // TrotterLiouv  += -4*1_i*sites.op("Sy", 1)*sites.op("Sy", 3);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sy", 2)*sites.op("Sy", 4);

 //    // TrotterLiouv  += -4*1_i*sites.op("Sz", 1)*sites.op("Sz", 3);
 //    // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 2)*sites.op("Sz", 4);

 //    auto g = Gate(sites,1,2,Gate::tReal,tstep/2.,TrotterLiouv );
 //    //gates.push_back(g);



 //    /* Construct local Hamiltonian + First Neighbor interactions for the interior of the chain */
 //    for(int j = 1; j <= 2*L-1; ++j)
 //    {
 //        TrotterLiouv  = -2*1_i*h*sites.op("Sz", j)*sites.op("Id", j+1);
 //        TrotterLiouv  +=  2*1_i*h*sites.op("Id", j)*sites.op("Sz", j+1);

 //        // TrotterLiouv  += -4*1_i*sites.op("Sx", 2*j-1)*sites.op("Sx", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*sites.op("Sx", 2*j)*sites.op("Sx", 2*j+2);

 //        // TrotterLiouv  += -4*1_i*sites.op("Sy", 2*j-1)*sites.op("Sy", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*sites.op("Sy", 2*j)*sites.op("Sy", 2*j+2);

 //        // TrotterLiouv  += -4*1_i*Delta*sites.op("Sz", 2*j-1)*sites.op("Sz", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 2*j)*sites.op("Sz", 2*j+2);

 //        g = Gate(sites,j,j+1,Gate::tReal,tstep/2.,TrotterLiouv );
 //        gates.push_back(g);
 //    }




 //    /* Construct D_L_m */
 //    TrotterLiouv  =  gamma*fL*sites.op("S-", 2*L-1)*sites.op("S-", 2*L);
 //    //TrotterLiouv  += -gamma*fL*0.5*sites.op("S+", 2*L-1)*sites.op("S-", 2*L-1);
	// TrotterLiouv  += -gamma*fL*0.5*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
	// TrotterLiouv  += -gamma*fL*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
 //    //TrotterLiouv  += -gamma*fL*0.5*sites.op("S+", 2*L)*sites.op("S-", 2*L);
	// TrotterLiouv  += -gamma*fL*0.5*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);
	// TrotterLiouv  += -gamma*fL*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);

 //    /* Construct D_L_p */
 //    TrotterLiouv  +=  gamma*(1.0-fL)*sites.op("S+", 2*L-1)*sites.op("S+", 2*L);
 //    //TrotterLiouv  += -gamma*(1.0-fL)*0.5*sites.op("S-", 2*L-1)*sites.op("S+", 2*L-1);
 //    TrotterLiouv  += -gamma*(1.0-fL)*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  -= -gamma*(1.0-fL)*0.5*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
 //    //TrotterLiouv  += -gamma*(1.0-fL)*0.5*sites.op("S-", 2*L)*sites.op("S+", 2*L);
 //    TrotterLiouv  += -gamma*(1.0-fL)*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  -= -gamma*(1.0-fL)*0.5*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);

 //    /* Local Hamiltonian for the last-site case */
 //    TrotterLiouv  += -2*1_i*h*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  +=  2*1_i*h*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);

 //    /* Nearest-neighbor interactions for the last-site case */
 //    // TrotterLiouv  += -4*1_i*sites.op("Sx", 2*L-1)*sites.op("Sx", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sx", 2*L)*sites.op("Sx", 2*L+2);

 //    // TrotterLiouv  += -4*1_i*sites.op("Sy", 2*L-1)*sites.op("Sy", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sy", 2*L)*sites.op("Sy", 2*L+2);

 //    // TrotterLiouv  += -4*1_i*Delta*sites.op("Sz", 2*L-1)*sites.op("Sz", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 2*L)*sites.op("Sz", 2*L+2);

 //    g = Gate(sites,(2*L-1),(2*L),Gate::tReal,tstep/2.,TrotterLiouv );
 //    //gates.push_back(g);

 // //    //FIRST PART ENDS HERE

 // //    /* Construct D_L_m */
 //    TrotterLiouv  =  gamma*fL*sites.op("S-", 2*L-1)*sites.op("S-", 2*L);
 //    //TrotterLiouv  += -gamma*fL*0.5*sites.op("S+", 2*L-1)*sites.op("S-", 2*L-1);
 //    TrotterLiouv  += -gamma*fL*0.5*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  += -gamma*fL*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
 //    //TrotterLiouv  += -gamma*fL*0.5*sites.op("S+", 2*L)*sites.op("S-", 2*L);
 //    TrotterLiouv  += -gamma*fL*0.5*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);
 //    TrotterLiouv  += -gamma*fL*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);

 //    /* Construct D_L_p */
 //    TrotterLiouv  +=  gamma*(1.0-fL)*sites.op("S+", 2*L-1)*sites.op("S+", 2*L);
 //    //TrotterLiouv  += -gamma*(1.0-fL)*0.5*sites.op("S-", 2*L-1)*sites.op("S+", 2*L-1);
 //    TrotterLiouv  += -gamma*(1.0-fL)*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  -= -gamma*(1.0-fL)*0.5*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
 //    //TrotterLiouv  += -gamma*(1.0-fL)*0.5*sites.op("S-", 2*L)*sites.op("S+", 2*L);
 //    TrotterLiouv  += -gamma*(1.0-fL)*0.5*(1/2)*sites.op("Id", 2*L-1)*sites.op("Id", 2*L);
	// TrotterLiouv  -= -gamma*(1.0-fL)*0.5*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);


 //    /* Local Hamiltonian for the last-site case */
 //    TrotterLiouv  += -2*1_i*h*sites.op("Sz", 2*L-1)*sites.op("Id", 2*L);
 //    TrotterLiouv  +=  2*1_i*h*sites.op("Id", 2*L-1)*sites.op("Sz", 2*L);

 //    /* Nearest-neighbor interactions for the last-site case */
 //    // TrotterLiouv  += -4*1_i*sites.op("Sx", 2*L-1)*sites.op("Sx", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sx", 2*L)*sites.op("Sx", 2*L+2);

 //    // TrotterLiouv  += -4*1_i*sites.op("Sy", 2*L-1)*sites.op("Sy", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*sites.op("Sy", 2*L)*sites.op("Sy", 2*L+2);

 //    // TrotterLiouv  += -4*1_i*Delta*sites.op("Sz", 2*L-1)*sites.op("Sz", 2*L+1);
 //    // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 2*L)*sites.op("Sz", 2*L+2);

 //    g = Gate(sites,(2*L-1),2*L,Gate::tReal,tstep/2.,TrotterLiouv );
 //    //gates.push_back(g);

 //        /* Construct local Hamiltonian + First Neighbor interactions for the interior of the chain */
 //        for(int j = 2*L-1; j >= 1; --j)
 //    {
 //        TrotterLiouv  = -2*1_i*h*sites.op("Sz", j)*sites.op("Id", j+1);
 //        TrotterLiouv  +=  2*1_i*h*sites.op("Id", j)*sites.op("Sz", j+1);

 //        // TrotterLiouv  += -4*1_i*sites.op("Sx", 2*j-1)*sites.op("Sx", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*sites.op("Sx", 2*j)*sites.op("Sx", 2*j+2);

 //        // TrotterLiouv  += -4*1_i*sites.op("Sy", 2*j-1)*sites.op("Sy", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*sites.op("Sy", 2*j)*sites.op("Sy", 2*j+2);

 //        // TrotterLiouv  += -4*1_i*Delta*sites.op("Sz", 2*j-1)*sites.op("Sz", 2*j+1);
 //        // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 2*j)*sites.op("Sz", 2*j+2);

 //        g = Gate(sites,j,j+1,Gate::tReal,tstep/2.,TrotterLiouv );
 //        gates.push_back(g);
 //    }


 //    /* Construct D_1_m */
 //    TrotterLiouv  = gamma*f1*sites.op("S-", 1)*sites.op("S-", 2);
 //    //TrotterLiouv  += -gamma*f1*0.5*sites.op("S+", 1)*sites.op("S-", 1);
 //    TrotterLiouv  += -gamma*f1*0.5*sites.op("Sz", 1)*sites.op("Id", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    //TrotterLiouv  += -gamma*f1*0.5*sites.op("S+", 2)*sites.op("S-", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*sites.op("Id", 1)*sites.op("Sz", 2);
 //    TrotterLiouv  += -gamma*f1*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);

	// /* Construct D_1_p */
 //    TrotterLiouv  +=  gamma*(1.0-f1)*sites.op("S+", 1)*sites.op("S+", 2);
 //    //TrotterLiouv  += -gamma*(1.0-f1)*0.5*sites.op("S-", 1)*sites.op("S+", 1);
 //    TrotterLiouv  += -gamma*(1.0-f1)*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    TrotterLiouv  -= -gamma*(1.0-f1)*0.5*sites.op("Sz", 1)*sites.op("Id", 2);
 //    //TrotterLiouv  += -gamma*(1.0-f1)*0.5*sites.op("S-", 2)*sites.op("S+", 2);
 //    TrotterLiouv  += -gamma*(1.0-f1)*0.5*(1/2)*sites.op("Id", 1)*sites.op("Id", 2);
 //    TrotterLiouv  -= -gamma*(1.0-f1)*0.5*sites.op("Id", 1)*sites.op("Sz", 2);

 //    /* Local Hamiltonian for the first-site case */
 //    TrotterLiouv  += -2*1_i*h*sites.op("Sz", 1)*sites.op("Id", 2);
 //    TrotterLiouv  +=  2*1_i*h*sites.op("Id", 1)*sites.op("Sz", 2);

 // //    /* Nearest-neighbor interactions for the first-site case */
 // //    // TrotterLiouv  += -4*1_i*sites.op("Sx", 1)*sites.op("Sx", 2);
 // //    // TrotterLiouv  +=  4*1_i*sites.op("Sx", 1)*sites.op("Sx", 2);

 // //    // TrotterLiouv  += -4*1_i*sites.op("Sy", 1)*sites.op("Sy", 2);
 // //    // TrotterLiouv  +=  4*1_i*sites.op("Sy", 1)*sites.op("Sy", 2);

 // //    // TrotterLiouv  += -4*1_i*sites.op("Sz", 1)*sites.op("Sz", 2);
 // //    // TrotterLiouv  +=  4*1_i*Delta*sites.op("Sz", 1)*sites.op("Sz", 2);

 //    g = Gate(sites,1,2,Gate::tReal,tstep/2.,TrotterLiouv );
 //    //gates.push_back(g);

//	return gates;
//} //TrotterLiouv Construct
