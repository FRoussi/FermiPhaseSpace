#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
//#include "FHSystem.h"


enum FHTunnelingInteractionType{
    TEnull,
    TEnearestNeigbors,
    TEdistanceInteraction,
    TEgiven
};

class FHTunnelingInteraction
{
public:


    FHTunnelingInteraction();
    FHTunnelingInteraction(FHLattice* lat);
    FHTunnelingInteraction(FHLattice* lat, Json::Value tunnJson);
    FHTunnelingInteraction(FHLattice* lat, double neighbourDistance, double tunnelingEnergy);
    FHTunnelingInteraction(FHLattice* lat, std::vector<double> tunnelingDistances, std::vector<double> tunnelingEnergies);

    arma::cx_dmat*              getTunnelingInteractionPointer();
    arma::cx_double             getTunnelingInteraction(uint iSite, uint jSite);
    FHTunnelingInteractionType       getTunnelingInteractionType();

    virtual                     std::vector<double>         getParameters(uint set);
    virtual Json::Value         writeTunnelingInteractionJson();        


    ~FHTunnelingInteraction();

private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               tunnelingInteraction;

    FHTunnelingInteractionType       tunnelingInteractionType;    
    

    std::vector<double>         parameterSet1;      //<- Energies
    std::vector<double>         parameterSet2;      //<- Distances


    void                        addTunneling(double distance, double energy);

};