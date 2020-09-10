#include "math.h"
#include <complex.h>
#include <vector> 

using namespace itensor;

// Constructs a Liouvillian given the relevant parameters and a SiteSet object
// Returns an AutoMPO object
AutoMPO 
LiouvXXZConstruct(
	double Delta, 
	SiteSet &sites, 
	std::vector<int> DissipatorsLoc, 
	std::vector<double> f, 
	std::vector<double> gamma, 
	std::vector<double> h)       																							
{                                                                                                       
    int L = sites.N();
    AutoMPO ampo = AutoMPO(sites);

    if (DissipatorsLoc.size() != f.size())
	{throw std::invalid_argument( "ERROR: Misleading dissipators temperature vector while constructing Liouvillian, please check your input. Aborting." );}
    if (DissipatorsLoc.size() != gamma.size())
	{throw std::invalid_argument( "ERROR: Misleading dissipators to gamma vector while constructing Liouvillian, please check your input. Aborting." );}
    if (L != h.size())
	{throw std::invalid_argument( "ERROR: Misleading magnetic field vector while constructing Liouvillian, please check your input. Aborting." );}

	std::cout << "\nConstructing Liouvillian for a spin chain of size " << L << " using " << DissipatorsLoc.size() <<" dissipators\n ";
	for (int j = 0; j < DissipatorsLoc.size(); ++j)
    {
		/* Constructs set of D_m interactions */
		ampo +=  gamma[j]*f[j],     "SmL", DissipatorsLoc[j], "SpR", DissipatorsLoc[j];
		ampo += -gamma[j]*f[j]*0.5, "SpSmL", DissipatorsLoc[j];
		ampo += -gamma[j]*f[j]*0.5, "SpSmR", DissipatorsLoc[j];

		/* Constructs set of D_p interactions */
		ampo +=  gamma[j]*(1.0-f[j]),     "SpL", DissipatorsLoc[j], "SmR", DissipatorsLoc[j];
		ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpL", DissipatorsLoc[j];
		ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpR", DissipatorsLoc[j];
	}

    /* Construct local Hamiltonian */
    for(int j = 1; j <= L; ++j)
    {
        ampo += -2*1_i*h[j-1], "SzL", j;
        ampo +=  2*1_i*h[j-1], "SzR", j;
    }

    /* Construct nearest-neighbor interactions */
    for(int j = 1; j < L; ++j)
    {
        ampo += -4*1_i, "SxL", j, "SxL", j+1;
        ampo +=  4*1_i, "SxR", j, "SxR", j+1;

        ampo += -4*1_i, "SyL", j, "SyL", j+1;
        ampo +=  4*1_i, "SyR", j, "SyR", j+1;

        ampo += -4*1_i*Delta, "SzL", j, "SzL", j+1;
        ampo +=  4*1_i*Delta, "SzR", j, "SzR", j+1;
    }

    return ampo;
} //LiouvConstruct

// Constructs a Liouvillian dagger given the relevant parameters and a SiteSet object (XXZ chain)
// Returns an AutoMPO object
AutoMPO 
LiouvXXZDConstruct(
	double Delta, 
	SiteSet &sites, 
	std::vector<int> DissipatorsLoc, 
	std::vector<double> f, 
	std::vector<double> gamma, 
	std::vector<double> h)          
{                                                                                                           							
    int L = sites.N();
    AutoMPO ampoD = AutoMPO(sites);

    if (DissipatorsLoc.size() != f.size())
	{throw std::invalid_argument( "ERROR: Misleading dissipators temperature vector while constructing Liouvillian dagger, please check your input. Aborting." );}
    if (DissipatorsLoc.size() != gamma.size())
	{throw std::invalid_argument( "ERROR: Misleading dissipators to gamma vector while constructing Liouvillian dagger, please check your input. Aborting." );}
    if (L != h.size())
	{throw std::invalid_argument( "ERROR: Misleading magnetic field vector while constructing Liouvillian dagger, please check your input. Aborting." );}

	std::cout << "\nConstructing Liouvillian dagger for a spin chain of size " << L << " using " << DissipatorsLoc.size() <<" dissipators\n ";
    for (int j = 0; j < DissipatorsLoc.size(); ++j)
    {
	    /* Constructs set of D_m interactions */
	   	ampoD +=  gamma[j]*f[j], "SpL", DissipatorsLoc[j], "SmR", DissipatorsLoc[j];
	    ampoD += -gamma[j]*f[j]*0.5, "SpSmL", DissipatorsLoc[j];
	    ampoD += -gamma[j]*f[j]*0.5, "SpSmR", DissipatorsLoc[j];

	    /* Constructs set of D_p interactions */
	    ampoD +=  gamma[j]*(1.0-f[j]), "SmL", DissipatorsLoc[j], "SpR", DissipatorsLoc[j];
	    ampoD += -gamma[j]*(1.0-f[j])*0.5, "SmSpL", DissipatorsLoc[j];
	    ampoD += -gamma[j]*(1.0-f[j])*0.5, "SmSpR", DissipatorsLoc[j];
    }

    /* Construct local Hamiltonian */
    for(int j = 1; j <= L; ++j)
    {
        ampoD +=  2*1_i*h[j-1], "SzL", j;
        ampoD += -2*1_i*h[j-1], "SzR", j;
    }

    /* Construct nearest-neighbor interactions */
    for(int j = 1; j < L; ++j)
    {
        ampoD +=  4*1_i, "SxL", j, "SxL", j+1;
        ampoD += -4*1_i, "SxR", j, "SxR", j+1;
        
        ampoD +=  4*1_i, "SyL", j, "SyL", j+1;
        ampoD += -4*1_i, "SyR", j, "SyR", j+1;
        
        ampoD +=  4*1_i*Delta, "SzL", j, "SzL", j+1;
        ampoD += -4*1_i*Delta, "SzR", j, "SzR", j+1;
    }

    return ampoD;
} //LiouvDConstruct


itensor::MPO 
LdLXXZConstruct(
    SiteSet &sites, 
	double Delta, 
	std::vector<int> dissipatorsVec, 
	std::vector<double> dissipatorsTempValues, 
	std::vector<double> gammaVec, 
	std::vector<double> hVec)
{
	int L = sites.N();

	itensor::AutoMPO LiouvAMPO = LiouvXXZConstruct(Delta, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec);
    itensor::AutoMPO LiouvDAMPO = LiouvXXZDConstruct(Delta, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec);
    
    itensor::MPO Liouv = MPO(LiouvAMPO);
    itensor::MPO LiouvD = MPO(LiouvDAMPO);

    //print("Liouv = ", Liouv);
    //print("LiouvD = ", LiouvD);


    itensor::MPO LdL;
    itensor::nmultMPO(Liouv, LiouvD, LdL, {"Maxm", 50000,"Cutoff",1E-8}); // constructs Ld*L
    //print("LdL = ", LdL);
   	return LdL;
}