#pragma once

#include "FHExternalPotential.h"



class FHExternalPotentialConstant : public FHExternalPotential
{
public:

    FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot);

    std::vector<double>         getParameters();
    Json::Value         writePotentialJson();   


    ~FHExternalPotentialConstant();

protected :

    


};