#include "../include/FHExternalPotential.h"


FHExternalPotential::FHExternalPotential()
{
    nullPotential = true;
    timeDependantPotential = false;
    currentPotentialTime = 0;
}

FHExternalPotential::FHExternalPotential(const arma::dmat* sitePos)
{
    uint nSites = sitePos->n_rows;

    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    sitePositions = sitePos;

    nullPotential = true;
    timeDependantPotential = false;
    currentPotentialTime = 0;
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
    return externalPotential(iSite, jSite);
}

arma::cx_double         FHExternalPotential::getExternalPotential(double t, uint iSite, uint jSite)
{
    if (currentPotentialTime != t)
    {
        updatePotential(t);
    }
    return externalPotential(iSite, jSite);
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