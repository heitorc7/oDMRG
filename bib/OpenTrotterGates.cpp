using Gate = itensor::BondGate<itensor::ITensor>;
std::vector<Gate> 
OpenTrotterXXZConstruct(
	double gamma, 
	double f1, 
	double fL, 
	double h, 
	double Delta, 
	SiteSet sites, 
	Real tstep)
{                                                                                                       				  
    int L = sites.N();                // relevant for looping through
    auto gates = std::vector<Gate>();
    
     //From 1 to L
    /* Construct D_1_m */
    auto ampo =  gamma*f1*sites.op("SmLSpR", 1)*sites.op("Id", 2);
    ampo += -gamma*f1*0.5*sites.op("SpSmL", 1)*sites.op("Id", 2);
    ampo += -gamma*f1*0.5*sites.op("SpSmR", 1)*sites.op("Id", 2);

    /* Construct D_1_p */
    ampo +=  gamma*(1.0-f1)*sites.op("SpLSmR", 1)*sites.op("Id", 2);
    ampo += -gamma*(1.0-f1)*sites.op("SmSpL", 1)*sites.op("Id", 2);
    ampo += -gamma*(1.0-f1)*sites.op("SmSpR", 1)*sites.op("Id", 2);

    /* Hamiltonian */
    ampo = -2*1_i*h*sites.op("SzL", 1)*sites.op("Id", 2);
    ampo +=  2*1_i*h*sites.op("SzR", 1)*sites.op("Id", 2);
    
    /* Nearest Neighbor interaction */
    ampo += -4*1_i*sites.op("SxL", 1)*sites.op("SxL", 2);
    ampo +=  4*1_i*sites.op("SxR", 1)*sites.op("SxR", 2);
    ampo += -4*1_i*sites.op("SyL", 1)*sites.op("SyL", 2);
    ampo +=  4*1_i*sites.op("SyR", 1)*sites.op("SyR", 2);
    ampo += -4*1_i*Delta*sites.op("SzL", 1)*sites.op("SzL", 2);
    ampo +=  4*1_i*Delta*sites.op("SzR", 1)*sites.op("SzR", 2);

    auto g = Gate(sites,1,2,Gate::tReal,tstep/2.,ampo);
    gates.push_back(g);

    /* Construct local Hamiltonian + Neighbor interaction */
    for(int j = 2; j < L-1; ++j)
    {
        ampo = -2*1_i*h*sites.op("SzL", j)*sites.op("Id", j+1);
        ampo +=  2*1_i*h*sites.op("SzR", j)*sites.op("Id", j+1);

        ampo += -4*1_i*sites.op("SxL", j)*sites.op("SxL", j+1);
        ampo +=  4*1_i*sites.op("SxR", j)*sites.op("SxR", j+1);
        ampo += -4*1_i*sites.op("SyL", j)*sites.op("SyL", j+1);
        ampo +=  4*1_i*sites.op("SyR", j)*sites.op("SyR", j+1);
        ampo += -4*1_i*Delta*sites.op("SzL", j)*sites.op("SzL", j+1);
        ampo +=  4*1_i*Delta*sites.op("SzR", j)*sites.op("SzR", j+1);

        auto g = Gate(sites,j,j+1,Gate::tReal,tstep/2.,ampo);
        gates.push_back(g);
    }

    /* Construct D_L_m */
    ampo =  gamma*fL*sites.op("Id", L-1)*sites.op("SmLSpR", L);
    ampo += -gamma*fL*0.5*sites.op("Id", L-1)*sites.op("SpSmL", L);
    ampo += -gamma*fL*0.5*sites.op("Id", L-1)*sites.op("SpSmR", L);

    /* Construct D_L_p */
    ampo +=  gamma*(1.0-fL)*sites.op("Id", L-1)*sites.op("SpLSmR", L);
    ampo += -gamma*(1.0-fL)*0.5*sites.op("Id", L-1)*sites.op("SmSpL", L);
    ampo += -gamma*(1.0-fL)*0.5*sites.op("Id", L-1)*sites.op("SmSpR", L);

    /* Hamiltonian */
    ampo += -2*1_i*h*sites.op("SzL", L-1)*sites.op("Id", L);
    ampo +=  2*1_i*h*sites.op("SzR", L-1)*sites.op("Id", L);

    /* Neighbor interaction */
    ampo += -4*1_i*sites.op("SxL", L-1)*sites.op("SxL", L);
    ampo +=  4*1_i*sites.op("SxR", L-1)*sites.op("SxR", L);
    ampo += -4*1_i*sites.op("SyL", L-1)*sites.op("SyL", L);
    ampo +=  4*1_i*sites.op("SyR", L-1)*sites.op("SyR", L);
    ampo += -4*1_i*Delta*sites.op("SzL", L-1)*sites.op("SzL", L);
    ampo +=  4*1_i*Delta*sites.op("SzR", L-1)*sites.op("SzR", L);

    g = Gate(sites,L-1,L,Gate::tReal,tstep/2.,ampo);
    gates.push_back(g);

    //From L to 1

    /* Construct D_L_m */
    ampo =  gamma*fL*sites.op("Id", L-1)*sites.op("SmLSpR", L);
    ampo += -gamma*fL*0.5*sites.op("Id", L-1)*sites.op("SpSmL", L);
    ampo += -gamma*fL*0.5*sites.op("Id", L-1)*sites.op("SpSmR", L);

    /* Construct D_L_p */
    ampo +=  gamma*(1.0-fL)*sites.op("Id", L-1)*sites.op("SpLSmR", L);
    ampo += -gamma*(1.0-fL)*0.5*sites.op("Id", L-1)*sites.op("SmSpL", L);
    ampo += -gamma*(1.0-fL)*0.5*sites.op("Id", L-1)*sites.op("SmSpR", L);

    ampo += -2*1_i*h*sites.op("SzL", L-1)*sites.op("Id", L);
    ampo +=  2*1_i*h*sites.op("SzR", L-1)*sites.op("Id", L);

    ampo += -4*1_i*sites.op("SxL", L-1)*sites.op("SxL", L);
    ampo +=  4*1_i*sites.op("SxR", L-1)*sites.op("SxR", L);
    ampo += -4*1_i*sites.op("SyL", L-1)*sites.op("SyL", L);
    ampo +=  4*1_i*sites.op("SyR", L-1)*sites.op("SyR", L);
    ampo += -4*1_i*Delta*sites.op("SzL", L-1)*sites.op("SzL", L);
    ampo +=  4*1_i*Delta*sites.op("SzR", L-1)*sites.op("SzR", L);

    g = Gate(sites,L-1,L,Gate::tReal,tstep/2.,ampo);
    gates.push_back(g);

    /* Construct local Hamiltonian in inverted order */
    for(int j = L-2; j >= 2; --j)
    {
        ampo = -2*1_i*h*sites.op("SzL", j)*sites.op("Id", j+1);
        ampo +=  2*1_i*h*sites.op("SzR", j)*sites.op("Id", j+1);

        ampo += -4*1_i*sites.op("SxL", j)*sites.op("SxL", j+1);
        ampo +=  4*1_i*sites.op("SxR", j)*sites.op("SxR", j+1);
        ampo += -4*1_i*sites.op("SyL", j)*sites.op("SyL", j+1);
        ampo +=  4*1_i*sites.op("SyR", j)*sites.op("SyR", j+1);
        ampo += -4*1_i*Delta*sites.op("SzL", j)*sites.op("SzL", j+1);
        ampo +=  4*1_i*Delta*sites.op("SzR", j)*sites.op("SzR", j+1);

        auto g = Gate(sites,j,j+1,Gate::tReal,tstep/2.,ampo);
        gates.push_back(g);
    }

    //From L to 1

    /* Construct D_1_m */
    ampo =  gamma*f1*sites.op("SmLSpR", 1)*sites.op("Id", 2);
    ampo += -gamma*f1*0.5*sites.op("SpSmL", 1)*sites.op("Id", 2);
    ampo += -gamma*f1*0.5*sites.op("SpSmR", 1)*sites.op("Id", 2);

    /* Construct D_1_p */
    ampo +=  gamma*(1.0-f1)*sites.op("SpLSmR", 1)*sites.op("Id", 2);
    ampo += -gamma*(1.0-f1)*sites.op("SmSpL", 1)*sites.op("Id", 2);
    ampo += -gamma*(1.0-f1)*sites.op("SmSpR", 1)*sites.op("Id", 2);

    /* Hamiltonian */
    ampo += -2*1_i*h*sites.op("SzL", 1)*sites.op("Id", 2);
    ampo +=  2*1_i*h*sites.op("SzR", 1)*sites.op("Id", 2);
    
    /* Nearest Neighbor interaction */
    ampo += -4*1_i*sites.op("SxL", 1)*sites.op("SxL", 2);
    ampo +=  4*1_i*sites.op("SxR", 1)*sites.op("SxR", 2);
    ampo += -4*1_i*sites.op("SyL", 1)*sites.op("SyL", 2);
    ampo +=  4*1_i*sites.op("SyR", 1)*sites.op("SyR", 2);
    ampo += -4*1_i*Delta*sites.op("SzL", 1)*sites.op("SzL", 2);
    ampo +=  4*1_i*Delta*sites.op("SzR", 1)*sites.op("SzR", 2);

    g = Gate(sites,1,2,Gate::tReal,tstep/2.,ampo);
    gates.push_back(g);

    return gates;
} //TrotterLiouv Construct