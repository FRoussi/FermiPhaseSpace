#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>

#include "json/json.h"
#include "FHLattice.h"
#include "FHCoreHole.h"
#include "MapDoubleEntry.h"


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

    const arma::cx_dmat&                  getTunnelingInteractionReference() const;
    arma::cx_double                 getTunnelingInteraction(uint iSite, uint jSite) const;
    FHTunnelingInteractionType      getTunnelingInteractionType() const;

    virtual std::vector<double>     getParameters(uint set) const;
    virtual Json::Value             writeTunnelingInteractionJson();        

    void                            addSite(arma::cx_dvec interactions);


    ~FHTunnelingInteraction();

private:  


    FHLattice*                  fhLattice;

    arma::cx_dmat               tunnelingInteraction;

    FHTunnelingInteractionType       tunnelingInteractionType{FHTunnelingInteractionType::TEnull};    


    std::vector<double>         parameterSet1{std::vector<double>()};      //<- Energies
    std::vector<double>         parameterSet2{std::vector<double>()};      //<- Distances

    FHTunnelingInteractionRelation  tunnelingInteractionRelation{FHTunnelingInteractionRelation::TERall};
    std::string speciesA;
    std::string speciesB;


    void                        addTunneling(double distance, double energy);
    void                        addTunnelingOrder1();
    double                      getTunnelingEnergyOrder1(double distance);
    void                        addTunnelingOrder2();
    double                      getTunnelingEnergyOrder2(double distance);


    MapDoubleEntry<FHTunnelingInteractionType, std::string>                   tunnelingInteractionTypeNames;
    MapDoubleEntry<FHTunnelingInteractionRelation, std::string>               tunnelingInteractionRelationNames;


    void                                    init_tunnelingMaps();

};