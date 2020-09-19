#include "MakeIvec.cpp"
#include "Liouv.cpp"
#include "WarmUp.cpp"
#include "itensor/all.h"
#include "MagneticProfile.cpp"
#include "FileManagement.cpp"
//#include "OpenTrotterGates.cpp"

void 
MakeIVEC(
		   itensor::MPS &Ivec, 
		   int L);

itensor::AutoMPO 
LiouvXXZConstruct(
    				double Delta, 
    				itensor::SiteSet &sites,
    				std::vector<int> DissipatorsLoc, 
    				std::vector<double> f, 
    				std::vector<double> gamma, 
    				std::vector<double> h);

itensor::AutoMPO 
LiouvXXZDConstruct(
    				double Delta, 
    				itensor::SiteSet &sites, 
    				std::vector<int> DissipatorsLoc, 
    				std::vector<double> f, 
    				std::vector<double> gamma, 
    				std::vector<double> h);

itensor::MPO 
LdLXXZConstruct(
    	itensor::SiteSet &sites, 
    	double Delta, 
    	std::vector<int> dissipatorsVec, 
    	std::vector<double> dissipatorsTempValues, 
    	std::vector<double> gammaVec, 
    	std::vector<double> hVec);

itensor::MPO 
LdLXXZConstruct(
        itensor::SiteSet &sites, 
        double Delta, 
        double f1,
        double fL, 
        double gamma,
        double h)
    {
        std::vector<int> dissipatorsVec;
        std::vector<double> dissipatorsTempValues;
        std::vector<double> gammaVec;
        std::vector<double> hVec;

        dissipatorsVec.push_back(1);
        dissipatorsVec.push_back(length(sites));

        dissipatorsTempValues.push_back(f1);
        dissipatorsTempValues.push_back(fL);

        gammaVec.push_back(gamma);
        gammaVec.push_back(gamma);

        for (int i = 0; i < length(sites); ++i)
        {
            hVec.push_back(h);
        }
        return LdLXXZConstruct(sites, Delta, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec);
    };

itensor::MPO 
LdLXXZConstruct(
        itensor::SiteSet &sites, 
        double Delta, 
        double f1,
        double fL, 
        double gamma,
        std::vector<double> hVec)
    {
        std::vector<int> dissipatorsVec;
        std::vector<double> dissipatorsTempValues;
        std::vector<double> gammaVec;

        dissipatorsVec.push_back(1);
        dissipatorsVec.push_back(length(sites));

        dissipatorsTempValues.push_back(f1);
        dissipatorsTempValues.push_back(fL);

        gammaVec.push_back(gamma);
        gammaVec.push_back(gamma);

        return LdLXXZConstruct(sites, Delta, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec);
    };

itensor::MPO 
LdLXXZConstruct(
        itensor::SiteSet &sites, 
        double Delta, 
        double f1,
        double fN, 
        double gamma1,
        double gammaN,
        double h)
    {
        std::vector<int> dissipatorsVec;
        std::vector<double> dissipatorsTempValues;
        std::vector<double> gammaVec;
        std::vector<double> hVec;

        dissipatorsVec.push_back(1);
        dissipatorsVec.push_back(length(sites));

        dissipatorsTempValues.push_back(f1);
        dissipatorsTempValues.push_back(fN);

        gammaVec.push_back(gamma1);
        gammaVec.push_back(gammaN);

        for (int i = 0; i < length(sites); ++i)
        {
            hVec.push_back(h);
        }

        return LdLXXZConstruct(sites, Delta, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec);
    };
    
void 
WarmUp(
    	itensor::MPS &rho_RLn, 
    	itensor::MPO &M2, 
    	double convergenceThreshold, 
    	const Args& args, 
    	const Args& argsDMRG);

std::vector<double> 
MagneticProfile(
    			itensor::SiteSet &sites, 
    			double hIni, 
    			double hFin, 
    			char ProfileTag);

void
PrintMagneticProfile(
    std::vector<double> hVec);                                                                                                          

void SaveMPS(
    	itensor::SiteSet &sites, 
    	itensor::MPS &chainMPS, 
    	double gamma, 
    	double f1, 
    	double fL, 
    	double h, 
    	double Delta, 
    	int sweepNum, 
    	int maxBond);

itensor::MPS LoadMPS(
    				itensor::SiteSet &sites, 
    				double gamma, 
    				double f1, 
    				double fL, 
    				double h, 
    				double Delta, 
    				int sweepNum, 
    				int maxBond);

/*std::vector<itensor::BondGate<itensor::ITensor>> OpenTrotterXXZConstruct(
                                                            double gamma, 
                                                            double f1, 
                                                            double fL, 
                                                            double h, 
                                                            double Delta, 
                                                            SiteSet sites, 
                                                            Real tstep);
*/

std::vector<double> 
ExtractValuesFromFile(
    				  std::string filename);

int 
smartBDchange(
    double sweepCont,
    int profile);