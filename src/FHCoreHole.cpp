#include "../include/FHCoreHole.h"


FHCoreHole::FHCoreHole()
{
    
}


FHCoreHole::FHCoreHole(Json::Value coreHoleJson) 
{
    siteNumber = coreHoleJson.get("site", 0).asUInt();
    holePotential = coreHoleJson.get("eps", -2).asDouble();
    holeTunneling = coreHoleJson.get("X", -1).asDouble();
    holeInteraction = coreHoleJson.get("G", 1).asDouble();

    downSpinOccupation = 0;
    upSpinOccupation = 1;

    Json::Value occ = coreHoleJson["IC"];
    if (occ.size()>0)         upSpinOccupation = occ[0].asUInt();
    if (occ.size()>1)         downSpinOccupation = occ[1].asUInt();
           

}


uint    FHCoreHole::getSiteNumber()         { return siteNumber;}
double  FHCoreHole::getHolePotential()      { return holePotential;}
double  FHCoreHole::getHoleTunneling()      { return holeTunneling;}
double  FHCoreHole::getHoleInteraction()    { return holeInteraction;}
uint    FHCoreHole::getUpSpinOccupation()   { return upSpinOccupation;}
uint    FHCoreHole::getDownSpinOccupation() { return downSpinOccupation;}


Json::Value FHCoreHole::writeFHCoreHoleJson()
{
    Json::Value jsonCH;
    jsonCH["presence"] = true;
    jsonCH["site"] = siteNumber;
    jsonCH["eps"] = holePotential;
    jsonCH["X"] = holeTunneling;
    jsonCH["G"] = holeInteraction;

    Json::Value occ(Json::arrayValue);
    occ.append(upSpinOccupation);
    occ.append(downSpinOccupation);
    jsonCH["IC"] = occ;

    return jsonCH;
}


FHCoreHole::~FHCoreHole()
{

}
