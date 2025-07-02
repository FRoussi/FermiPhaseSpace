#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
//#include "FHSystem.h"


enum FHCoulombInteractionType{
    CEnull,
    CEonSite,
    CEdistanceInteraction,
    CEgiven
};

class FHCoulombInteraction
{
public:


    FHCoulombInteraction();
    FHCoulombInteraction(FHLattice* lat);
    FHCoulombInteraction(FHLattice* lat, Json::Value coulJson);
    FHCoulombInteraction(FHLattice* lat, double onSiteInteraction);
    FHCoulombInteraction(FHLattice* lat, std::vector<double> coulombDistances, std::vector<double> coulombEnergies);

    arma::cx_dmat*                  getCoulombUInteractionPointer();
    arma::cx_double                 getCoulombUInteraction(uint iSite, uint jSite);
    arma::cx_dmat*                  getCoulombVInteractionPointer();
    arma::cx_double                 getCoulombVInteraction(uint iSite, uint jSite);
    FHCoulombInteractionType        getCoulombInteractionType();

    virtual                     std::vector<double>         getParameters(uint set);
    virtual Json::Value         writeCoulombInteractionJson();                 


    ~FHCoulombInteraction();

private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               U;      //<- Opposite spin Interaction
    arma::cx_dmat               V;      //<- Same spin Interaction

    FHCoulombInteractionType       CoulombInteractionType;    
    

    std::vector<double>         parameterSet1;   //<- Energies
    std::vector<double>         parameterSet2;  //<- Distances


    void                        addCoulomb(double distance, double Interaction);
    void                        addCoulombinU(double distance, double Interaction);
    void                        addCoulombinV(double distance, double Interaction, bool fillDiagonal = false);


};