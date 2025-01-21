#pragma once

#include "FHExternalPotential.h"



class FHExternalPotentialConstant : public FHExternalPotential
{
public:

    FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot);
    FHExternalPotentialConstant(FHLattice* lat, uint iSite, double intensityA0);


    ~FHExternalPotentialConstant();

protected :


};