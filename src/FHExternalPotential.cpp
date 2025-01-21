#include "../include/FHExternalPotential.h"


FHExternalPotential::FHExternalPotential()
{
    nullPotential = true;
    timeDependantPotential = false;
    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::null;
}

FHExternalPotential::FHExternalPotential(FHLattice* lat)
{


    uint nSites = lat->getnSites();

    
    //externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    externalPotential = arma::cx_dmat();

    fhLattice = lat;

    nullPotential = true;
    timeDependantPotential = false;
    currentPotentialTime = 0;
    potentialType = FHExternalPotentialType::null;
}

bool                    FHExternalPotential::isThePotentialNull()
{
    return nullPotential;
}

bool                    FHExternalPotential::isThePotentialTimeDependant()
{
    return timeDependantPotential;
}

bool                    FHExternalPotential::getExternalPotential(double t, arma::cx_dmat& potential)
{
    if (!nullPotential)
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
    if (potentialType == FHExternalPotentialType::null)
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
    jpot["type"] = "None";
    return jpot;
}