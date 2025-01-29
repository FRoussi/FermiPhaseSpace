#include "../include/FHCoulombEnergy.h"


FHCoulombEnergy::FHCoulombEnergy()
{
    CoulombEnergyType = FHCoulombEnergyType::CEnull;
}

FHCoulombEnergy::FHCoulombEnergy(FHLattice* lat)
{
    CoulombEnergyType = FHCoulombEnergyType::CEnull;
    fhLattice = lat;
}

FHCoulombEnergy::FHCoulombEnergy(FHLattice* lat, double onSiteEnergy)
{
    CoulombEnergyType = FHCoulombEnergyType::CEonSite;
    fhLattice = lat;

    parameterSet1 = std::vector<double>();
    parameterSet2 = std::vector<double>();

    parameterSet1.push_back(onSiteEnergy);



    uint nSites = fhLattice->getnSites();
    U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    addCoulomb(0, onSiteEnergy);


}


FHCoulombEnergy::FHCoulombEnergy(FHLattice* lat, std::vector<double> CoulombDs, std::vector<double> CoulombEs)
{
    CoulombEnergyType = FHCoulombEnergyType::CEdistanceEnergy;
    fhLattice = lat;

    parameterSet1 = std::vector<double>(CoulombDs);
    parameterSet2 = std::vector<double>(CoulombEs);

    uint nSites = fhLattice->getnSites();
    U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    uint nParam = std::max(parameterSet1.size(), parameterSet2.size());

    for (uint iDist = 0; iDist < nParam; iDist++)
    {
        addCoulomb(parameterSet1[iDist], parameterSet2[iDist]);
    }

}



arma::cx_dmat*          FHCoulombEnergy::getCoulombUEnergyPointer()
{
    return &U;
}

arma::cx_dmat*          FHCoulombEnergy::getCoulombVEnergyPointer()
{
    return &V;
}


arma::cx_double         FHCoulombEnergy::getCoulombUEnergy(uint iSite, uint jSite)
{
    if (CoulombEnergyType == FHCoulombEnergyType::CEnull)
    {
        return 0;
    }
    else 
    {
        return U(iSite, jSite);
    }

}

arma::cx_double         FHCoulombEnergy::getCoulombVEnergy(uint iSite, uint jSite)
{
    if (CoulombEnergyType == FHCoulombEnergyType::CEnull)
    {
        return 0;
    }
    else 
    {
        return V(iSite, jSite);
    }

}




std::vector<double>         FHCoulombEnergy::getParameters(uint set)
{
    if (set == 1)
    {
        return parameterSet1;
    }
    else if (set == 2)
    {
        return parameterSet2;
    }
    else
    {
        return std::vector<double> ();
    }

}

FHCoulombEnergyType FHCoulombEnergy::getCoulombEnergyType()
{
    return CoulombEnergyType;
}

Json::Value         FHCoulombEnergy::writeTunnelinEnergyJson()
{
    Json::Value jpot;
    switch(CoulombEnergyType)
    {
        case(FHCoulombEnergyType::CEnull):
        {
            jpot["type"] = "null";
        }
        case(FHCoulombEnergyType::CEonSite):
        {
            jpot["type"] = "onSite";
            jpot["onSiteInteraction"] = parameterSet1[0];
        }
        case(FHCoulombEnergyType::CEdistanceEnergy):
        {
            jpot["type"] = "distanceEnergy";

            Json::Value dVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet1.size(); i++)
            {
                dVecJ.append(parameterSet1[i]);
            }
            jpot["distances"] = dVecJ;

            Json::Value eVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet2.size(); i++)
            {
                eVecJ.append(parameterSet2[i]);
            }
            jpot["energies"] = eVecJ;
        }
        default:
        {
            jpot["type"] = "null";
        }
    }

    return jpot;
}




void                        FHCoulombEnergy::addCoulomb(double distance, double energy)
{
    double ths = fhLattice->getSamePointDistance();

    uint nSites = fhLattice->getnSites();

    for (uint iSite=0; iSite<nSites; iSite++)
    {
        for (uint jSite = iSite+1; jSite<nSites; jSite++)
        {
            if (std::norm(fhLattice->getDistance(iSite, jSite) - distance) < ths)
            {
                U(iSite, jSite) += energy;
                U(jSite, iSite) += energy;
            }
        }
    }


}