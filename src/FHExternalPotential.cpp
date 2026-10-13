#include "../include/FHExternalPotential.h"


FHExternalPotential::FHExternalPotential()
{
    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPnull;
}

FHExternalPotential::FHExternalPotential(FHLattice* lat)
{

    //uint nSites = lat->getnSites();    
    //externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    
    externalPotential = arma::cx_dmat();

    fhLattice = lat;

    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPnull;
}


FHExternalPotential::FHExternalPotential(FHLattice* lat, arma::dvec Pot)
{
    fhLattice = lat;
    
    uint nSites = fhLattice->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    arma::dvec diagPot = arma::dvec(nSites);
    diagPot.subvec(0,arma::size(Pot)) = Pot; 
    arma::cx_dvec diagPotC = arma::cx_dvec(diagPot, arma::dvec(nSites));

    externalPotential = arma::diagmat(diagPotC);


    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPconstant;
}

FHExternalPotential::FHExternalPotential(FHLattice* lat, uint iSite, double pot)
{

    fhLattice = lat;
    
    uint nSites = fhLattice->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    
    if (iSite < nSites)
    {
        externalPotential.at(iSite, iSite) = pot;
    }


    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPconstant;
}



bool                    FHExternalPotential::isThePotentialNull() const
{
    return (potentialType == FHExternalPotentialType::EPnull);
}

bool                    FHExternalPotential::isThePotentialTimeDependant() const
{
    return (potentialType == FHExternalPotentialType::EPtimeDependant);
}

bool                    FHExternalPotential::getExternalPotential(double t, arma::cx_dmat& potential)
{
    if (potentialType != FHExternalPotentialType::EPnull)
    {
        if (currentPotentialTime != t)
        {
            updatePotential(t);
            currentPotentialTime = t;
        }
        potential = externalPotential;
        return true;
    }
    else
    {
        return false;
    }
}


const arma::cx_dmat&          FHExternalPotential::getExternalPotentialReference() const 
{
    return externalPotential;
}

arma::cx_double         FHExternalPotential::getExternalPotential(uint iSite, uint jSite) const
{
    if (potentialType == FHExternalPotentialType::EPnull)
    {
        return 0;
    }
    else 
    {
        return externalPotential(iSite, jSite);
    }

}

arma::cx_double         FHExternalPotential::getExternalPotential(double t, uint iSite, uint jSite)
{
    if (currentPotentialTime != t)
    {
        updatePotential(t);
    }

    return getExternalPotential(iSite, jSite);

}


arma::cx_dmat        FHExternalPotential::getExternalPotential(double t) const
{
    arma::cx_dmat extPot = externalPotential;

    return extPot;
}

void                    FHExternalPotential::updatePotential(double t)
{
    currentPotentialTime = t;
}


std::vector<double>         FHExternalPotential::getParameters() const
{
    std::vector<double> diagPot = std::vector<double>();

    if (potentialType == FHExternalPotentialType::EPconstant)
    {
        for (uint iS = 0; iS < externalPotential.n_cols; iS++)
        {
            diagPot.push_back(externalPotential.at(iS, iS).real());
        }
    }

    return diagPot;
}

FHExternalPotentialType FHExternalPotential::getExternalPotentialType() const
{
    return potentialType;
}

Json::Value         FHExternalPotential::writePotentialJson()
{
    Json::Value jpot;
    if (potentialType == FHExternalPotentialType::EPconstant)
    {
        jpot["type"] = "constant";
        Json::Value diagP(Json::arrayValue);
        for (uint iS = 0; iS < externalPotential.n_cols; iS++)
        {
            diagP.append(externalPotential.at(iS, iS).real());
        }
        jpot["diagonal"] = diagP;
    }
    else
    {
        jpot["type"] = "null";
    }
    return jpot;
}


void    FHExternalPotential::addSite(double pot)
{
//std::cout << "fhEst" << pot << std::endl;
    if ((potentialType == FHExternalPotentialType::EPnull ) && (pot == 0))
    {

    }
    else
    {
        //std::cout << "fhEst 1"  << std::endl;
        if (potentialType == FHExternalPotentialType::EPnull)
        {

            uint nSites = fhLattice->getnSites();
                      //  std::cout << "fhEst 2"  << nSites << std::endl;
            externalPotential = arma::cx_dmat(nSites+1, nSites+1, arma::fill::zeros);
            externalPotential.at(nSites, nSites) = pot;
            potentialType = FHExternalPotentialType::EPconstant;
        }   
        else
        {
            uint nSites = fhLattice->getnSites();
            //std::cout << "fhEst 3"  << nSites << std::endl;
            arma::cx_dmat newPot = arma::cx_dmat(nSites+1, nSites+1, arma::fill::zeros);
            newPot.submat(0, 0, arma::size(externalPotential)) = externalPotential;
            externalPotential = newPot;
            externalPotential.at(nSites, nSites) = pot;
        }

    }

}