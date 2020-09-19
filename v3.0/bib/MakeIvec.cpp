#include "math.h"
#include <complex.h>
#include <vector> 

using namespace itensor;

void 
MakeIVEC(
	itensor::MPS &Ivec, int L)                        
{
    auto A =Ivec.A(1);
    auto indi = A.inds();
    
    //First and last sites are set by hand
    auto B = ITensor(indi[0],indi[1]);
    B.set(indi[0](1),indi[1](1),1);
    B.set(indi[0](2),indi[1](1),0);
    B.set(indi[0](3),indi[1](1),0);
    B.set(indi[0](4),indi[1](1),1);
    
    Ivec.setA(1,B);
    
    for(int i=2;i<L;i++){
        auto A =Ivec.A(i);
        auto indi = A.inds();
        
        auto B = ITensor(indi[0],indi[1],indi[2]);
        
        B.set(indi[0](1),indi[1](1),indi[2](1),1);
        B.set(indi[0](1),indi[1](2),indi[2](1),0);
        B.set(indi[0](1),indi[1](3),indi[2](1),0);
        B.set(indi[0](1),indi[1](4),indi[2](1),1);
        
        Ivec.setA(i,B);
    }
    
    A =Ivec.A(L);
    indi = A.inds();
    B = ITensor(indi[0],indi[1]);
    
    B.set(indi[0](1),indi[1](1),1);
    B.set(indi[0](1),indi[1](2),0);
    B.set(indi[0](1),indi[1](3),0);
    B.set(indi[0](1),indi[1](4),1);
    
    Ivec.setA(L,B);
} //MakeIVEC