#include "../include/FHExternalPotentialCoreHole.h"





FHExternalPotentialCoreHole::FHExternalPotentialCoreHole(FHLattice* lat, uint iSite, double pot): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    externalPotential(iSite, iSite) = std::complex<double>(pot);
    nullPotential = false;
    potentialType = FHExternalPotentialType::coreHole;
}


std::vector<double>         FHExternalPotentialCoreHole::getParameters()
{
    std::vector<double> v = {(double)iSite, holePotential};

    
    return v;
}


Json::Value        FHExternalPotentialCoreHole::writePotentialJson()
{
   Json::Value jpot;

   jpot["type"] = "coreHole";
   jpot["site"] = iSite;
   jpot["holePotential"] = holePotential;
    
    return jpot;
}

FHExternalPotentialCoreHole::~FHExternalPotentialCoreHole()
{

}
