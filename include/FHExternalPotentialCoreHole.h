#pragma once

#include "FHExternalPotential.h"



class FHExternalPotentialCoreHole : public FHExternalPotential
{
public:

    FHExternalPotentialCoreHole(FHLattice* lat, uint iSite, double pot);

    std::vector<double>         getParameters();
    Json::Value                 writePotentialJson();   


    ~FHExternalPotentialCoreHole();

protected :


    uint    iSite;
    double  holePotential;


};