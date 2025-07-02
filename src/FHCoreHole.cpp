#include "../include/FHCoreHole.h"



FHCoreHole::FHCoreHole(Json::Value coreHoleJson) 
{
    siteNumber = coreHoleJson.get("site", 0).asUInt();
    holePotential = coreHoleJson.get("eps", -2).asDouble();
    holeTunneling = coreHoleJson.get("X", -1).asDouble();
    holeInteraction = coreHoleJson.get("G", 1).asDouble();


}


uint    FHCoreHole::getSiteNumber() {return siteNumber;}
double  FHCoreHole::getHolePotential() {return holePotential;}
double  FHCoreHole::getHoleTunneling() {return holeTunneling;}
double  FHCoreHole::getHoleInteraction() {return holeInteraction;}


FHCoreHole::~FHCoreHole()
{

}
