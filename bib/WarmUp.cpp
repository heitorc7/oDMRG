void WarmUp(
	itensor::MPS &rho_L, 
	itensor::MPO &M2, 
	double convergenceThreshold = 0.001, 
	const itensor::Args& args = {"Quiet", false}, 
	const itensor::Args& argsDMRG = {"Quiet", true})
{	
	const bool quiet = args.getBool("Quiet");
	std::clock_t begin = clock(); 
    auto sweeps = Sweeps(1);
    sweeps.maxm() = 1;
    sweeps.cutoff() = 1E-8;
    bool tag = true;
    double energyIni;
    double energyFin = 100000;

    if (!quiet){
        std::cout << "\n WarmingUp! \n";}

    for (int count = 1; tag; count += 1)
    {
        energyIni = energyFin;
        energyFin = dmrg(rho_L,M2,sweeps,argsDMRG);
        std::cout << count;

        if (energyIni-energyFin < convergenceThreshold)
        {
        	tag = false;
        	if (!quiet)
        	{
        		std::cout << "\n\nEndedWarmUpAfter " << count << " sweeps \n\n ";
    			std::cout << "timeElapsedInWarmUpRound " << double(clock() - begin) / CLOCKS_PER_SEC;
    		}//writing out after ending warm-up
        }//warmUp ending condition
    }//warmUp looping

} //warmUp

int smartBDchange(
    double sweepCont,
    int profile)
{   
    if (profile == 1)
    {
        if (sweepCont > 50)
        {
            return 2;
        }

        if (sweepCont > 150)
        {
            return 5;
        }

        if (sweepCont > 500)
        {
            return 10;
        }

        if (sweepCont > 1000)
        {
            return 25;
        }
    }
    return 1;

} //smart is an overstatement