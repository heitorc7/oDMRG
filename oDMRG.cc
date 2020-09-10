#include "math.h"
#include <complex.h>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <string>
#include <climits>

#include "itensor/all.h"
#include "bib/all.h"

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double gamma = atof(argv[2]);
	double f1 = atof(argv[3]);
	double fN = atof(argv[4]);
	double h = atof(argv[5]);
	double Delta = atof(argv[6]);
    double hIni = h;
    double hFin = -h;
    int DissipatorNum = 2;


  char magProf = 'c';
  // 'c' for constant,  'l' for linear, 'p' or for parabolic, 
  // 's' for staggered, 'h' for first half h, last half -h

  //std::vector<double> energy;
  std::vector<itensor::MPO> obsCurrVec;
  std::vector<itensor::MPO> obsMagVec;

////////////////////////////////////////////////
// SIMULATION INTERNAL PARAMETERS
    int MaxBond = 3;
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

    auto sites = LRN(N);			// loading sites in (LR)N formalism
    MPS rho = MPS(sites);           // creating MPS object for density matrix rho
    MakeIVEC(rho, N);				// start from Id

    MPS Ivec = MPS(sites);          // creates Ivec as MPS object related to sites to be input in the MakeIVEC function
    MakeIVEC(Ivec, N);              // sets Ivec as desired for calculations


    /*************************************************************************************
    *            Pre-loading spin flux and magnetization operators in vectors            *
    **************************************************************************************/
    //  Current measurement loop start
    for(int i=1; i<N; i++)
    {
        auto aobs  = AutoMPO(sites);
        aobs +=  4.0,"SxL",i,"SyL",i+1;
        aobs += -4.0,"SyL",i,"SxL",i+1;
        obsCurrVec.push_back(MPO(aobs));
        //print(aobs);
    }   // Current measurement loop end

    //  Magnetization measurement loop start
    for(int i=1; i<=N; i++)
    {
        auto aobs  = AutoMPO(sites);
        aobs += 2.0,"SzL",i;

        obsMagVec.push_back(MPO(aobs));
        //print(aobs);
    }   // Magnetization measurement loop end


    std::vector<double> hVec = MagneticProfile(sites, hIni, hFin, magProf);		//constructing magnetic field profile
    PrintMagneticProfile(hVec);

    MPO LdL = LdLXXZConstruct(sites, Delta, f1, fN, gamma, hVec);		//constructing LdL by  L(dag)*L

	/************************************************
    *       		 	 WARM UP		            *
    *************************************************/
    std::clock_t begin = clock();

    WarmUp(rho, LdL, 0.001, {"Quiet", false});

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
    *					   					    *
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
        for(int j=1;j<=N;j++)
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

        /* Stopping condition */
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
