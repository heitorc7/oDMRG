void 
PrintSysInfo(
    int N,
    double gamma,
    double f1,
    double fN,
    double hIni,
    double hFin,
    double Delta,
    int MaxBond,
    bool saveMPS,
    bool loadMPS,
    int DissipatorNum,
    int loadedSweepNum,
    int loadedMaxBD
    )
{
    std::cout << "N:           " << N << "\n";
    std::cout << "gamma:       " << gamma << "\n";
    std::cout << "f1:          " << f1 << "\n";
    std::cout << "fL:          " << fN << "\n";
    std::cout << "h1:          " << hIni << "\n";
    std::cout << "hL:          " << hFin << "\n";
    std::cout << "Delta:       " << Delta<< "\n";
    std::cout << "bondMax:     " << MaxBond << "\n";
    std::cout << "SaveMPS:     " << saveMPS << "\n";
    std::cout << "LoadMPS:     " << loadMPS << "\n";
    std::cout << "DissipatorsNumber: " << DissipatorNum << "\n";

    if (loadMPS)
     {
     std::cout << "LoadedPreviousMPS \n";
     std::cout << "LoadedSweepNum: " << loadedSweepNum << "\n";
     std::cout << "LoadedMaxBD: " << loadedMaxBD << "\n";
     }

 }

void 
SaveMPS(
	itensor::SiteSet &sites, 
	itensor::MPS &chainMPS, 
	double gamma, 
	double f1, 
	double fL, 
	double h, 
	double Delta, 
	int sweepNum, 
	int maxBond)
{
	std::cout << "\n\n Saving MPS \n";

	int L = length(sites);

    std::stringstream sitesString; 
    std::stringstream psiString; 

    sitesString << "Sites-FileL_L-" << L << "_gamma-" << gamma << "_f1-" << f1 << "_fL-" << fL << "_h-" << h << "_Delta-" << Delta << "_SweepsDone-" << sweepNum << "_MaxBondDimension-" << maxBond;
    psiString   << "Psi-FileL_L-" << L << "_gamma-" << gamma << "_f1-" << f1 << "_fL-" << fL << "_h-" << h << "_Delta-" << Delta << "_SweepsDone-" << sweepNum << "_MaxBondDimension-" << maxBond;
    const std::string sitesTmp = sitesString.str(); 
    const char* sitesName = sitesTmp.c_str();

    const std::string psiTmp = psiString.str(); 
    const char* psiName = psiTmp.c_str();

    writeToFile(sitesName, sites);
    writeToFile(psiName,chainMPS); 
	std::cout << "\n MPS sucessfully saved \n\n";
}

/* NOT WORKING DO NOT USE IT */
itensor::MPS 
LoadMPS(itensor::SiteSet &sites, 
	double gamma, 
	double f1, 
	double fL, 
	double h, 
	double Delta, 
	int sweepNum, 
	int maxBond)
{
	int L = length(sites); 

	std::cout << "\n\n Loading MPS from \n" 
	<< "Psi-FileL_L-" << L << "_gamma-" << gamma << "_f1-" << f1 << "_fL-" << fL << "_h-" << h << "_Delta-" << Delta 
	<< "_SweepsDone-" << sweepNum << "_MaxBondDimension-" << maxBond <<"\n\n";

	std::stringstream sitesString; 
    std::stringstream psiString; 

    sitesString << "Sites-FileL_L-" << L << "_gamma-" << gamma << "_f1-" << f1 << "_fL-" << fL << "_h-" << h << "_Delta-" << Delta << "_SweepsDone-" << sweepNum << "_MaxBondDimension-" << maxBond;
    psiString   << "Psi-FileL_L-" << L << "_gamma-" << gamma << "_f1-" << f1 << "_fL-" << fL << "_h-" << h << "_Delta-" << Delta << "_SweepsDone-" << sweepNum << "_MaxBondDimension-" << maxBond;
    const std::string sitesTmp = sitesString.str(); const char* sitesName = sitesTmp.c_str();
    const std::string psiTmp = psiString.str(); const char* psiName = psiTmp.c_str();


    readFromFile(sitesName, sites);
    MPS chainMPS = MPS(sites);
    readFromFile(psiName, chainMPS);
    return chainMPS;
    std::cout << "\n MPS sucessfully Loaded. Previously, "<< sweepNum <<" were performed, and the final bond dimension was " << maxBond <<"\n\n";
}

// Returns a vector with information extracted from file under filename
std::vector<double> 
ExtractValuesFromFile(
	std::string filename)
{
	std::vector<double> outVec;

    std::string STRING;
    std::ifstream infile;

    infile.open (filename);
    int cont = 0;

    std::vector<double> values;
    std::vector<int> dissipatorsVec;
    std::vector<double> dissipatorsTempValues;
    std::vector<double> gammaVec;
    std::vector<double> hVec;

    while(!infile.eof()) // To get you all the lines.
    {
        cont += 1;
        getline(infile,STRING); // Saves the line in STRING.
        std::string::size_type spaceVal = STRING.find(" " , 0)+1;
        if (STRING == "end")
        {
            break;
        }
        values.push_back(stod(STRING.substr(spaceVal)));

        if (cont > 11 && cont < 12 + values[10])
        {
            dissipatorsVec.push_back(stod(STRING.substr(1)));
            dissipatorsTempValues.push_back(stod(STRING.substr(spaceVal)));
            gammaVec.push_back(values[1]);
        }

    }
    	outVec.push_back(/*L = */values[0]);
    	outVec.push_back(/*gamma = */values[1]);
    	outVec.push_back(/*f1 = */values[2]);
    	outVec.push_back(/*fL = */values[3]);    
    	outVec.push_back(/*hIni = */values[4]);
    	outVec.push_back(/*hFin = */values[5]);
    	outVec.push_back(/*Delta = */values[6]);
    	outVec.push_back(/*MaxBond =  */values[7]);
    	outVec.push_back(/*saveMPS = */values[8]);
    	outVec.push_back(/*loadMPS = */values[9]);
    	outVec.push_back(/*DissipatorNum = */values[10]);

    if (values[9])
    {
        outVec.push_back(/*loadedSweepNum = */values[11 + values[10]]);
        outVec.push_back(/*loadedMaxBD = */values[12 + values[10]]);
	}
	else
	{
		outVec.push_back(0);
		outVec.push_back(0);
	}

	outVec.insert(outVec.end(), dissipatorsVec.begin(), dissipatorsVec.end());
	outVec.insert(outVec.end(), dissipatorsTempValues.begin(), dissipatorsTempValues.end());
	outVec.insert(outVec.end(), gammaVec.begin(), gammaVec.end());

	return outVec;
}