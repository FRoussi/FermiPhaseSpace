#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "FHSystem.h"




class FHExternalPotential
{
public:

    FHExternalPotential();
    FHExternalPotential(const arma::dmat* sitePos);

    bool                    getExternalPotential(double time, arma::cx_dmat& potential);
    arma::cx_dmat*          getExternalPotentialPointer();
    arma::cx_double         getExternalPotential(uint iSite, uint jSite);
    arma::cx_double         getExternalPotential(double t, uint iSite, uint jSite);

protected:

    bool                    isThereAPotential;
    
    double                  currentPotentialTime;

    const arma::dmat*       sitePositions;



    arma::cx_dmat           externalPotential;

    void                    updatePotential(double t);

};