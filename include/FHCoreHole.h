#pragma once

#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

#include "json/json.h"



class FHCoreHole
{
public:

    FHCoreHole();
    FHCoreHole(Json::Value coreHoleJson);


    uint    getSiteNumber();
    double  getHolePotential();
    double  getHoleTunneling();
    double  getHoleInteraction();
    double  getHoleSiteInteraction();
    double  getUpSpinOccupation();
    double  getDownSpinOccupation();

    Json::Value writeFHCoreHoleJson();

    ~FHCoreHole();


private:

   // bool    presence;
    uint    siteNumber;
    double  holePotential;
    double  holeTunneling;
    double  holeInteraction;
    double  holeSiteInteraction;
    double  upSpinOccupation;
    double  downSpinOccupation;       


};