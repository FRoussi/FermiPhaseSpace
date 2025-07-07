#include "../include/FHExternalPotentialConstant.h"



FHExternalPotentialConstant::FHExternalPotentialConstant(FHLattice* lat): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    potentialType = FHExternalPotentialType::EPconstant;
}


FHExternalPotentialConstant::FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    externalPotential.diag() = arma::cx_dmat(Pot, arma::dvec(Pot.size(), arma::fill::zeros));
    potentialType = FHExternalPotentialType::EPconstant;
}




std::vector<double>         FHExternalPotentialConstant::getParameters()
{

    std::vector<double> v = std::vector<double>();
    for (uint i = 0; i<externalPotential.n_cols ; i++)
    {
        v.push_back(externalPotential(i,i).real());
    }

    return v;
}


Json::Value        FHExternalPotentialConstant::writePotentialJson()
{
   Json::Value jpot;

   jpot["type"] = "constant";
    Json::Value jsonPotVec(Json::arrayValue);
    for (uint i = 0; i < externalPotential.n_cols; i++)
    {
        jsonPotVec.append(externalPotential(i,i).real());
    }  
   jpot["sitePotentials"] = jsonPotVec;
    
    return jpot;
}

void    FHExternalPotentialConstant::addSite(double pot)
{
   // std::cout << "fhEstConst" << std::endl;
    uint nSites = externalPotential.n_cols;

    arma::cx_dmat newPot = arma::cx_dmat(nSites+1, nSites+1, arma::fill::zeros);


    newPot.submat(0, 0, nSites-1, nSites - 1) = externalPotential;
    newPot.at(nSites, nSites) = pot;

}

FHExternalPotentialConstant::~FHExternalPotentialConstant()
{

}
