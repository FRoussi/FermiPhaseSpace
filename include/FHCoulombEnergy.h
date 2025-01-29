#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
//#include "FHSystem.h"


enum FHCoulombEnergyType{
    CEnull,
    CEonSite,
    CEdistanceEnergy,
    CEgiven
};

class FHCoulombEnergy
{
public:


    FHCoulombEnergy();
    FHCoulombEnergy(FHLattice* lat);
    FHCoulombEnergy(FHLattice* lat, double onSiteEnergy);
    FHCoulombEnergy(FHLattice* lat, std::vector<double> coulombDistances, std::vector<double> coulombEnergies);

    arma::cx_dmat*              getCoulombUEnergyPointer();
    arma::cx_double             getCoulombUEnergy(uint iSite, uint jSite);
    arma::cx_dmat*              getCoulombVEnergyPointer();
    arma::cx_double             getCoulombVEnergy(uint iSite, uint jSite);
    FHCoulombEnergyType         getCoulombEnergyType();

    virtual                     std::vector<double>         getParameters(uint set);
    virtual Json::Value         writeTunnelinEnergyJson();                 


private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               U;      //<- Opposite spin Interaction
    arma::cx_dmat               V;      //<- Same spin Interaction

    FHCoulombEnergyType       CoulombEnergyType;    
    

    std::vector<double>         parameterSet1;
    std::vector<double>         parameterSet2;


    void                        addCoulomb(double distance, double energy);

};