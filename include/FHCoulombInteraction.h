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

    const arma::cx_dmat&                  getCoulombUInteractionReference() const;
    arma::cx_double                 getCoulombUInteraction(uint iSite, uint jSite) const;
    const arma::cx_dmat&                  getCoulombVInteractionReference() const;
    arma::cx_double                 getCoulombVInteraction(uint iSite, uint jSite) const;
    FHCoulombInteractionType        getCoulombInteractionType() const;

    virtual                         std::vector<double>         getParameters(uint set) const;
    virtual Json::Value             writeCoulombInteractionJson();                 

    void                            addSite(arma::cx_dvec interactions);

    ~FHCoulombInteraction();

private:  


    FHLattice*                      fhLattice;

    arma::cx_dmat                   U;      //<- Opposite spin Interaction
    arma::cx_dmat                   V;      //<- Same spin Interaction

    FHCoulombInteractionType        CoulombInteractionType{FHCoulombInteractionType::CEnull};    
    

    std::vector<double>             parameterSet1{std::vector<double>()};   //<- Energies
    std::vector<double>             parameterSet2{std::vector<double>()};  //<- Distances


    void                            addCoulomb(double distance, double Interaction);
    void                            addCoulombinU(double distance, double Interaction);
    void                            addCoulombinV(double distance, double Interaction, bool fillDiagonal = false);


};