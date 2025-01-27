#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
//#include "FHSystem.h"


enum FHTunnelingEnergyType{
    TEnull,
    TEnearestNeigbors,
    TEdistanceEnergy,
    TEgiven
};

class FHTunnelingEnergy
{
public:


    FHTunnelingEnergy();
    FHTunnelingEnergy(FHLattice* lat);
    FHTunnelingEnergy(FHLattice* lat, double neighbourDistance, double tunnelingEnergy);
    FHTunnelingEnergy(FHLattice* lat, std::vector<double> tunnelingDistances, std::vector<double> tunnelingEnergies);

    arma::cx_dmat*              getTunnelingEnergyPointer();
    arma::cx_double             getTunnelingEnergy(uint iSite, uint jSite);
    FHTunnelingEnergyType       getTunnelingEnergyType();

    virtual                     std::vector<double>         getParameters(uint set);
    virtual Json::Value         writeTunnelinEnergyJson();                 


private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               tunnelingEnergy;

    FHTunnelingEnergyType       tunnelingEnergyType;    
    

    std::vector<double>         parameterSet1;
    std::vector<double>         parameterSet2;


    void                        addTunneling(double distance, double energy);

};