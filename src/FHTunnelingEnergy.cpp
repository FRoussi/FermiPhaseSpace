#include "../include/FHTunnelingEnergy.h"


FHTunnelingEnergy::FHTunnelingEnergy()
{
    tunnelingEnergyType = FHTunnelingEnergyType::TEnull;
}

FHTunnelingEnergy::FHTunnelingEnergy(FHLattice* lat)
{
    tunnelingEnergyType = FHTunnelingEnergyType::TEnull;
    fhLattice = lat;
}

FHTunnelingEnergy::FHTunnelingEnergy(FHLattice* lat, double neighbourD, double tunnelingE)
{
    tunnelingEnergyType = FHTunnelingEnergyType::TEnearestNeigbors;
    fhLattice = lat;

    parameterSet1 = std::vector<double>();
    parameterSet2 = std::vector<double>();

    parameterSet1.push_back(neighbourD);
    parameterSet2.push_back(tunnelingE);


    uint nSites = fhLattice->getnSites();
    tunnelingEnergy = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    addTunneling(neighbourD, tunnelingE);


}


FHTunnelingEnergy::FHTunnelingEnergy(FHLattice* lat, std::vector<double> tunnelingDs, std::vector<double> tunnelingEs)
{
    tunnelingEnergyType = FHTunnelingEnergyType::TEdistanceEnergy;
    fhLattice = lat;

    parameterSet1 = std::vector<double>(tunnelingDs);
    parameterSet2 = std::vector<double>(tunnelingEs);

    uint nSites = fhLattice->getnSites();
    tunnelingEnergy = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    uint nParam = std::max(parameterSet1.size(), parameterSet2.size());

    for (uint iDist = 0; iDist < nParam; iDist++)
    {
        addTunneling(parameterSet1[iDist], parameterSet2[iDist]);
    }

}



arma::cx_dmat*          FHTunnelingEnergy::getTunnelingEnergyPointer()
{
    return &tunnelingEnergy;
}

arma::cx_double         FHTunnelingEnergy::getTunnelingEnergy(uint iSite, uint jSite)
{
    if (tunnelingEnergyType == FHTunnelingEnergyType::TEnull)
    {
        return 0;
    }
    else 
    {
        return tunnelingEnergy(iSite, jSite);
    }

}




std::vector<double>         FHTunnelingEnergy::getParameters(uint set)
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

FHTunnelingEnergyType FHTunnelingEnergy::getTunnelingEnergyType()
{
    return tunnelingEnergyType;
}

Json::Value         FHTunnelingEnergy::writeTunnelinEnergyJson()
{
    Json::Value jpot;
    switch(tunnelingEnergyType)
    {
        case(FHTunnelingEnergyType::TEnull):
        {
            jpot["type"] = "null";
        }
        case(FHTunnelingEnergyType::TEnearestNeigbors):
        {
            jpot["type"] = "nearestNeighbors";
            jpot["neighborDistance"] = parameterSet1[0];
            jpot["neighborEnergy"] = parameterSet2[0];
        }
        case(FHTunnelingEnergyType::TEdistanceEnergy):
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




void                        FHTunnelingEnergy::addTunneling(double distance, double energy)
{
    double ths = fhLattice->getSamePointDistance();

    uint nSites = fhLattice->getnSites();

    for (uint iSite=0; iSite<nSites; iSite++)
    {
        for (uint jSite = iSite+1; jSite<nSites; jSite++)
        {
            if (std::norm(fhLattice->getDistance(iSite, jSite) - distance) < ths)
            {
                tunnelingEnergy(iSite, jSite) += energy;
                tunnelingEnergy(jSite, iSite) += energy;
            }
        }
    }


}