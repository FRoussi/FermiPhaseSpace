#include "FHExternalPotential.h"


FHExternalPotential::FHExternalPotential()
{
    isThereAPotential = false;
    currentPotentialTime = 0;
}

FHExternalPotential::FHExternalPotential(const arma::dmat* sitePos)
{
    uint nSites = sitePos->n_cols;

    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    sitePositions = sitePos;

    isThereAPotential = false;
    currentPotentialTime = 0;
}

bool FHExternalPotential::getExternalPotential(double t, arma::cx_dmat& potential)
{
    if (isThereAPotential)
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

void FHExternalPotential::updatePotential(double t)
{
    currentPotentialTime = t;
}