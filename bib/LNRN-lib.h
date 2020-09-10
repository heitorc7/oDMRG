// Here are indexed functions in bib.cpp. If you want to add some, don't forget to include them here!
// For information on what they do check LNRN-funcs.cpp!
#ifndef LNRN
#define LNRN
#include "itensor/all.h"
#include "LNRN-funcs.cpp"

    void MakeIVEC_LNRN(itensor::MPS &Ivec, int L);

    itensor::AutoMPO LiouvConstruct_LNRN(double gamma, double f1, double fL, double h, double Delta, itensor::SiteSet sites);
    itensor::AutoMPO LiouvDConstruct_LNRN(double gamma, double f1, double fL, double h, double Delta, itensor::SiteSet sites);

#endif
