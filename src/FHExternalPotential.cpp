#include "../include/FHExternalPotential.h"


FHExternalPotential::FHExternalPotential()
{
    timeDependantPotential = false;
    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPnull;
}

FHExternalPotential::FHExternalPotential(FHLattice* lat)
{

    //uint nSites = lat->getnSites();    
    //externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    
    externalPotential = arma::cx_dmat();

    fhLattice = lat;

    timeDependantPotential = false;
    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::EPnull;
}

bool                    FHExternalPotential::isThePotentialNull()
{
    return (potentialType == FHExternalPotentialType::EPnull);
}

bool                    FHExternalPotential::isThePotentialTimeDependant()
{
    return timeDependantPotential;
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


arma::cx_dmat*          FHExternalPotential::getExternalPotentialPointer()
{
    return &externalPotential;
}

arma::cx_double         FHExternalPotential::getExternalPotential(uint iSite, uint jSite)
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


arma::cx_dmat        FHExternalPotential::getExternalPotential(double t)
{
    arma::cx_dmat extPot = externalPotential;

    return extPot;
}

void                    FHExternalPotential::updatePotential(double t)
{
    currentPotentialTime = t;
}


std::vector<double>         FHExternalPotential::getParameters()
{
    return std::vector<double>();
}

FHExternalPotentialType FHExternalPotential::getExternalPotentialType()
{
    return potentialType;
}

Json::Value         FHExternalPotential::writePotentialJson()
{
    Json::Value jpot;
    jpot["type"] = "null";
    return jpot;
}


void    FHExternalPotential::addSite(double pot)
{

    if ((potentialType == FHExternalPotentialType::EPnull ) || (pot == 0))
    {

    }
    else
    {
        std::cout << "maybe use a constant potential" << std::endl;

    }

}