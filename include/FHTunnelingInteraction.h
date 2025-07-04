#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
#include "FHCoreHole.h"
//#include "FHSystem.h"


enum FHTunnelingInteractionType{
    TEnull,
    TEnearestNeigbors,
    TEdistanceInteraction0,
    TEdistanceInteraction1,
    TEdistanceInteraction2,
    TEgiven,
    TEsum
};

enum FHTunnelingInteractionRelation{
    TERall,
    TERsameSpecies,
    TERAB,
};



class FHTunnelingInteraction
{
public:


    FHTunnelingInteraction();
    FHTunnelingInteraction(FHLattice* lat);
    FHTunnelingInteraction(FHLattice* lat, Json::Value tunnJson);
    FHTunnelingInteraction(FHLattice* lat, FHCoreHole* coreHole, Json::Value tunnJson);
    FHTunnelingInteraction(FHLattice* lat, double neighbourDistance, double tunnelingEnergy);
    FHTunnelingInteraction(FHLattice* lat, std::vector<double> tunnelingDistances, std::vector<double> tunnelingEnergies);

    arma::cx_dmat*              getTunnelingInteractionPointer();
    arma::cx_double             getTunnelingInteraction(uint iSite, uint jSite);
    FHTunnelingInteractionType       getTunnelingInteractionType();

    virtual                     std::vector<double>         getParameters(uint set);
    virtual Json::Value         writeTunnelingInteractionJson();        

    void                        addSite(arma::cx_dvec interactions);


    ~FHTunnelingInteraction();

private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               tunnelingInteraction;

    FHTunnelingInteractionType       tunnelingInteractionType;    





    std::vector<double>         parameterSet1;      //<- Energies
    std::vector<double>         parameterSet2;      //<- Distances

    FHTunnelingInteractionRelation  tunnelingInteractionRelation;
    std::string speciesA;
    std::string speciesB;


    void                        addTunneling(double distance, double energy);
    void                        addTunnelingOrder1();
    double                      getTunnelingEnergyOrder1(double distance);
    void                        addTunnelingOrder2();
    double                      getTunnelingEnergyOrder2(double distance);


    std::map<FHTunnelingInteractionType, std::string>           tunnelingInteractionTypeNames;
    std::map<std::string, FHTunnelingInteractionType>           tunnelingInteractionTypeNamesReverse;
    std::map<FHTunnelingInteractionRelation, std::string>       tunnelingInteractionRelationNames;
    std::map<std::string, FHTunnelingInteractionRelation>       tunnelingInteractionRelationNamesReverse;


    void                                    init_tunnelingMaps();

};