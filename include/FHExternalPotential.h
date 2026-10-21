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
    const arma::cx_dmat&        getExternalPotentialReference() const;
    arma::cx_double             getExternalPotential(uint iSite, uint jSite) const;
    arma::cx_double             getExternalPotential(double t, uint iSite, uint jSite);
    arma::cx_dmat               getExternalPotential(double t) const;
    FHExternalPotentialType     getExternalPotentialType() const;

    std::vector<double>         getParameters() const;
    Json::Value                 writePotentialJson();                 

    void                        updatePotential(double t);

    bool                        isThePotentialNull() const;
    bool                        isThePotentialTimeDependant() const;

    void                        addSite(double pot);

protected:  



    double                      currentPotentialTime{0};

    FHLattice*                  fhLattice;


    arma::cx_dmat               externalPotential;

    FHExternalPotentialType     potentialType{FHExternalPotentialType::EPnull};    
    

};