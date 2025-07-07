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