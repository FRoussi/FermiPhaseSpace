#include "../include/FHExternalPotentialConstant.h"


FHExternalPotentialConstant::FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    externalPotential.diag() = arma::cx_dmat(Pot, arma::dvec(Pot.size(), arma::fill::zeros));
    nullPotential = false;
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

FHExternalPotentialConstant::~FHExternalPotentialConstant()
{

}
