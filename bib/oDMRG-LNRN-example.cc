#include "math.h"
#include <complex.h>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <string>
#include <climits>

#include "itensor/all.h"
#include "all.h"

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double gamma = atof(argv[2]);
	double f1 = atof(argv[3]);
	double fN = atof(argv[4]);
	double h = atof(argv[5]);
	double Delta = atof(argv[6]);
    double hIni = h;
    double hFin = h;
    int DissipatorNum = 2;

  //char magProf = 'c';
  // 'c' for constant, 'l' for linear, 'p' for parabolic, 's' for staggered, 'h' for half h, half -h

  std::vector<double> energy;
  std::vector<itensor::MPO> obsCurrVec;
  std::vector<itensor::MPO> obsMagVec;

////////////////////////////////////////////////
// SIMULATION INTERNAL PARAMETERS
    int MaxBond = -1;
    int sweepNumber = 1;
    double sweepBDChangeThresholdValue = 0.001; // If the difference between the energy before & after "sweepNumber" sweeps falls
    int BDinc = 1;                             // under this value, the bond dimension is increased by BDinc

    int bdIni = 1;                             // Choose another initial bond dimension  if you wish

    int loadedSweepNum = 0;
    int loadedMaxBD = 0;
    bool saveMPS = false;
    bool loadMPS = false;
    bool importInputFromFile = false;
    int sweepIni = 1;
    double timeElapsed;
////////////////////////////////////////////////

    PrintSysInfo(N, gamma, f1, fN, hIni, hFin, Delta, MaxBond, saveMPS, loadMPS, DissipatorNum, loadedSweepNum, loadedMaxBD);

    SiteSet sites = SpinHalf(2*N);			// loading sites in (LR)N formalism
    MPS rho = MPS(sites);       // creating MPS object for density matrix rho
    MakeIVEC_LNRN(rho, N);				// start from Id

    MPS Ivec = MPS(sites);          // creates Ivec as MPS object related to sites to be input in the MakeIVEC function
    MakeIVEC_LNRN(Ivec, N);              // sets Ivec as desired for calculations


    /*************************************************************************************
    *            Pre-loading spin flux and magnetization operators in vectors            *
    **************************************************************************************/
    //  Current measurement loop start
    for(int i=1; i<N; i++)
    {
        auto aobs  = AutoMPO(sites);
        aobs +=  4.0,"Sx",2*i-1,"Sy",2*i+1;
        aobs += -4.0,"Sy",2*i-1,"Sx",2*i+1;
        obsCurrVec.push_back(MPO(aobs));
    }   // Current measurement loop end

    //  Magnetization measurement loop start
    for(int i=1; i<=2*N; i++)
    {
        auto aobs  = AutoMPO(sites);
        aobs += 2.0,"Sz",i;
        obsMagVec.push_back(MPO(aobs));
    }   // Magnetization measurement loop end

    //std::vector<double> hVec = MagneticProfile(sites, hIni, hFin, magProf);		//constructing magnetic field profile
    //PrintMagneticProfile(hVec);


    AutoMPO LiouvAMPO = LiouvConstruct_LNRN(gamma, f1, fN, h, Delta, sites);
    AutoMPO LiouvDAMPO = LiouvDConstruct_LNRN(gamma, f1, fN, h, Delta, sites);
    MPO Liouv = MPO(LiouvAMPO);
    MPO LiouvD = MPO(LiouvDAMPO);
    MPO LdL;
    nmultMPO(Liouv, LiouvD, LdL, {"Maxm", 50000,"Cutoff",1E-8}); // constructs Ld*L


	/************************************************
    *       		 	 WARM UP		            *
    *************************************************/
    std::clock_t begin = clock();

    WarmUp(rho, LdL, 0.001, {"Quiet", true});

    timeElapsed = double(clock() - begin) / CLOCKS_PER_SEC;
    std::cout << "timeElapsedInWarmUpRound " << timeElapsed;


    /*************************************************
    *      		   SWEEPS PARAMETERS                 *
    **************************************************/
    int bd = bdIni;
    auto sweeps = Sweeps(sweepNumber);
    sweeps.maxm() = bd;
    sweeps.cutoff() = 1E-8;
    double energyIni; double energyFin = LONG_MAX;


    /********************************************
    *					   						*
    *             DMRG + MEASUREMENTS           *
    *					   						*
    *********************************************/
    bool stopCond = true;
    for (int i = sweepIni; stopCond; i += sweepNumber)
    {

        /**********************************
        *              DMRG               *
        ***********************************/

            energyIni = energyFin;
            energyFin = dmrg(rho,LdL,sweeps,{"Quiet",true});


        /*********************************
        *          MEASUREMENTS          *
        **********************************/
        std::cout << "\n\n" << i << "\n\n";

        // Sz measurement (Magnon Density) loop for each site start
        for(int j=1;j<=2*N;j++)
        {
            auto ev = overlapC(Ivec,obsMagVec[j-1],rho)/overlapC(Ivec,rho);
            std::cout << "\nMagnonDensityAtSite " << j << " is " << ev;
        }   // Sz measurement (Magnon Density) loop for each site end
        std::cout << "\n\n";

        //  Current measurement loop start
        for(int j=1; j<N; j++)
        {
            auto current = overlapC(Ivec,obsCurrVec[j-1],rho)/overlapC(Ivec,rho);
            std::cout << "\nCurrentAtSite " << j << "-" << j+1 << " is " << current;
        }   // Current measurement loop end

        std::cout << "\n\n\nEnergyAfterSweeps " << i << " is " << energyFin << "\n";
        std::cout << "\nBondAfterSweeps " << i << " is " << bd << "\n";

        if (MaxBond > 0 && bd > MaxBond)
        {
            stopCond = false;
        }

        /* Checking for BD increase */
        if ((energyIni-energyFin)/energyIni < sweepBDChangeThresholdValue)
        {
            bd += BDinc;
            sweeps.maxm() = bd;
        }

    }//LOOP FOR SWEEPS - END
    return 0;
}   // main end
