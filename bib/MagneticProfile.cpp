#include "math.h"
#include <complex.h>
#include <vector> 

using namespace itensor;

// Simple function to create different magnetic profiles to your chain. 
// Use ProfileTag to choose if you want a linear profile ('l') or constant ('c')
std::vector<double> 
MagneticProfile(
	itensor::SiteSet &sites, 
	double hIni, 
	double hFin , 
	char ProfileTag)
{
	int L = sites.N();
	std::vector<double> outVec; 

	// Constant magnetic field profile
	if (ProfileTag == 'c')
	{
		for (int cont = 0; cont < L; ++cont)
    	{
       		outVec.push_back(hIni);
    	}
	}

	// Linear magnetic field profile
	if (ProfileTag == 'l')
	{
		for (int cont = 0; cont < L; ++cont)
    	{
        	outVec.push_back(hIni + cont*(hFin - hIni)/(L-1));
    	}
	}

	// Parabolic magnetic field profile
	if (ProfileTag == 'p')
	{
		for (int cont = 1; cont <= L; ++cont)
    	{
        	outVec.push_back( 
        		( 4*(hIni - hFin)/(-L*L + (L+1)*2*L + 4 - (L+1)*4) )*cont*cont + ( (4*(L+1)*(hFin-hIni))/(-L*L + (L+1)*2*L + 4 - (L+1)*4) )*cont + (hIni - ( 4*(hIni - hFin)/(-L*L + (L+1)*2*L + 4 - (L+1)*4) ) - ( (4*(L+1)*(hFin-hIni))/(-L*L + (L+1)*2*L + 4 - (L+1)*4) ) )
        	);
    	}
	}

	// Staggered magnetic field profile
	if (ProfileTag == 's')
	{
		for (int cont = 0; cont < L; ++cont)
    	{
        	outVec.push_back( hIni * pow(-1, cont) );
    	}
	}

	// +h for half of the chain, -h for the other half magnetic profile
	if (ProfileTag == 'h')
	{
		for (int cont = 1; cont <= L; ++cont)
    	{
    		if (cont <= L/2)
    		{
    			outVec.push_back( hIni );
    		}
    		else
    		{
    			outVec.push_back( -hIni );
    		}
        	
    	}
	}

	//else{
	//		printf("I don't recognize this tag. Sorry. Possible tags are (lower or uppercase): 'c' for constant, 'l' for linear, 'p' for parabolic, 's' for staggered, 'h' for half'n'half\n");
	//	}
			return outVec;
}

void
PrintMagneticProfile(
	std::vector<double> hVec)
{
	std::cout << "hProfile: ";
    for (int i = 0; i < hVec.size(); ++i)
    {
    	std::cout << hVec[i] <<", ";
    }
    std::cout << "\n\n ";

}