#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
#include "FHCoreHole.h"
#include "FHTunnelingInteraction.h"
//#include "FHSystem.h"


enum FHExternalPotentialType{
    EPnull,
    EPconstant,
    EPtimeDependant
};

class FHExternalPotential
{
public:


    FHExternalPotential();
    FHExternalPotential(FHLattice* lat);
    FHExternalPotential(FHLattice* lat, arma::dvec Pot);
    FHExternalPotential(FHLattice* lat, uint iSite, double pot);


    bool                        getExternalPotential(double time, arma::cx_dmat& potential);
    arma::cx_dmat*              getExternalPotentialPointer();
    arma::cx_double             getExternalPotential(uint iSite, uint jSite);
    arma::cx_double             getExternalPotential(double t, uint iSite, uint jSite);
    arma::cx_dmat               getExternalPotential(double t);
    FHExternalPotentialType     getExternalPotentialType();

    std::vector<double>         getParameters();
    Json::Value                 writePotentialJson();                 

    void                        updatePotential(double t);

    bool                        isThePotentialNull();
    bool                        isThePotentialTimeDependant();

    void                        addSite(double pot);

protected:  



    double                      currentPotentialTime;

    FHLattice*                  fhLattice;


    arma::cx_dmat               externalPotential;

    FHExternalPotentialType     potentialType;    
    

};