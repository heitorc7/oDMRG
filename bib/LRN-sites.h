#ifndef __ITENSOR_TWOSPINHALF_H
#define __ITENSOR_TWOSPINHALF_H
#include "itensor/mps/siteset.h"
#include "math.h"

namespace itensor {

class TwoSpinHalfSite;
using LRN = BasicSiteSet<TwoSpinHalfSite>;

class TwoSpinHalfSite
    {
    IQIndex s;
    public:

    TwoSpinHalfSite() { }

    TwoSpinHalfSite(IQIndex _I) : s(_I) { }

    TwoSpinHalfSite(int n, Args const& args = Args::global())
        {
        s = IQIndex{nameint("S=1/2 ",n),
               Index(nameint("UpUp ",n),1,Site),QN(),
               Index(nameint("UpDn ",n),1,Site),QN(),
               Index(nameint("DnUp ",n),1,Site),QN(),
               Index(nameint("DnDn ",n),1,Site),QN()};
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "UpUp")
            {
            return s(1);
            }
        else 
        if(state == "UpDn")
            {
            return s(2);
            }
        else
        if(state == "DnUp")
        {
            return s(3);
        }
        else
            if(state == "DnDn")
            {
                return s(4);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
            auto sP = prime(s);

            auto UpUp = s(1);
            auto UpDn = s(2);
            auto DnUp = s(3);
            auto DnDn = s(4);

            auto UpUpP = sP(1);
            auto UpDnP = sP(2);
            auto DnUpP = sP(3);
            auto DnDnP = sP(4);

            auto Op = IQTensor(dag(s),sP);

            if(opname == "SzL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,UpUpP,0.5);
                Op.set(UpDn,UpDnP,-0.5);
                Op.set(DnUp,DnUpP,0.5);
                Op.set(DnDn,DnDnP,-0.5);
            }
            else
            if(opname == "SzR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,UpUpP,0.5);
                Op.set(UpDn,UpDnP,0.5);
                Op.set(DnUp,DnUpP,-0.5);
                Op.set(DnDn,DnDnP,-0.5);
            }

            else
            if(opname == "SxL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpDn,UpUpP,0.5);
                Op.set(UpUp,UpDnP,0.5);
                Op.set(DnDn,DnUpP,0.5);
                Op.set(DnUp,DnDnP,0.5);
            }

            else
            if(opname == "SxR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(DnUp,UpUpP,0.5);
                Op.set(DnDn,UpDnP,0.5);
                Op.set(UpUp,DnUpP,0.5);
                Op.set(UpDn,DnDnP,0.5);
            }

            else
            if(opname == "SyL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpDn,UpUpP,-0.5*Cplx_i);
                Op.set(UpUp,UpDnP,0.5*Cplx_i);
                Op.set(DnDn,DnUpP,-0.5*Cplx_i);
                Op.set(DnUp,DnDnP,0.5*Cplx_i);
            }

            else
            if(opname == "SyR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(DnUp,UpUpP,0.5*Cplx_i);
                Op.set(DnDn,UpDnP,0.5*Cplx_i);
                Op.set(UpUp,DnUpP,-0.5*Cplx_i);
                Op.set(UpDn,DnDnP,-0.5*Cplx_i);
            }

            else
            if(opname == "SmL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,UpDnP,1.);
                Op.set(DnUp,DnDnP,1.);
            }

            else
            if(opname == "SmR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(DnUp,UpUpP,1.);
                Op.set(DnDn,UpDnP,1.);
            }           

            else
            if(opname == "SpL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpDn,UpUpP,1.);
                Op.set(DnDn,DnUpP,1.);
            }

            else
            if(opname == "SpR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,DnUpP,1.);
                Op.set(UpDn,DnDnP,1.);
            } 

            else
            if(opname == "SmLSpR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,DnDnP,1.);
            }

            else
            if(opname == "SpLSmR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(DnDn,UpUpP,1.);
            }            

            else
            if(opname == "SpSmL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,UpUpP,1.);
                Op.set(DnUp,DnUpP,1.);
            }

            else
            if(opname == "SpSmR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpUp,UpUpP,1.);
                Op.set(UpDn,UpDnP,1.);
            }

            else
            if(opname == "SmSpL")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(UpDn,UpDnP,1.);
                Op.set(DnDn,DnDnP,1.);
            }

            else
            if(opname == "SmSpR")
            {
//                Op = mixedIQTensor(dag(s),sP);
                Op.set(DnUp,DnUpP,1.);
                Op.set(DnDn,DnDnP,1.);
            }

            return Op;
        }
    };

} //namespace itensor

#endif
