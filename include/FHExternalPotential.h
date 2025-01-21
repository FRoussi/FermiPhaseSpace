#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
//#include "FHSystem.h"


enum FHExternalPotentialType{
    noneEP,
    null,
    constant,
    coreHole,
    pulse
};

class FHExternalPotential
{
public:


    FHExternalPotential();
    FHExternalPotential(FHLattice* lat);

    bool                        getExternalPotential(double time, arma::cx_dmat& potential);
    arma::cx_dmat*              getExternalPotentialPointer();
    arma::cx_double             getExternalPotential(uint iSite, uint jSite);
    arma::cx_double             getExternalPotential(double t, uint iSite, uint jSite);
    arma::cx_dmat               getExternalPotential(double t);
    FHExternalPotentialType     getExternalPotentialType();

    virtual                     std::vector<double>         getParameters();
    virtual Json::Value         writePotentialJson();                 

    void                        updatePotential(double t);

    bool                        isThePotentialNull();
    bool                        isThePotentialTimeDependant();

protected:  

    bool                        nullPotential;

    bool                        timeDependantPotential;

    double                      currentPotentialTime;

    FHLattice*            fhLattice;


    arma::cx_dmat               externalPotential;

    FHExternalPotentialType     potentialType;    
    

};