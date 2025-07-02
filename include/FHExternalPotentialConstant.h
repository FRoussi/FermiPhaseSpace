#pragma once

#include "FHExternalPotential.h"



class FHExternalPotentialConstant : public FHExternalPotential
{
public:

    FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot);

    std::vector<double>         getParameters();
    void                        addSite(double pot);

    Json::Value                 writePotentialJson();   



    ~FHExternalPotentialConstant();

protected :

    


};