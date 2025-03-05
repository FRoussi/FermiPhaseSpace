#include "../include/FHTunnelingInteraction.h"


FHTunnelingInteraction::FHTunnelingInteraction()
{
    tunnelingInteractionType = FHTunnelingInteractionType::TEnull;
}

FHTunnelingInteraction::FHTunnelingInteraction(FHLattice* lat)
{
    tunnelingInteractionType = FHTunnelingInteractionType::TEnull;
    fhLattice = lat;
}


FHTunnelingInteraction::FHTunnelingInteraction(FHLattice* lat, Json::Value tunnJson)
{
    tunnelingInteractionType = FHTunnelingInteractionType::TEnull;
    fhLattice = lat;
    std::string tunnType = tunnJson.get("type", "null").asString();
    if (tunnType == "nearestNeighbors")
    {

        tunnelingInteractionType = FHTunnelingInteractionType::TEnearestNeigbors;
        double tunnDist = tunnJson.get("neighborDistance", 1.).asDouble();
        double tunnE = tunnJson.get("neighborInteraction", -1.).asDouble();

    parameterSet1 = std::vector<double>();
    parameterSet2 = std::vector<double>();

    parameterSet1.push_back(tunnE);
    parameterSet2.push_back(tunnDist);

        uint nSites = fhLattice->getnSites();
        tunnelingInteraction = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

            addTunneling(tunnDist, tunnE);

    }

}

FHTunnelingInteraction::FHTunnelingInteraction(FHLattice* lat, double neighbourD, double tunnelingE)
{

    tunnelingInteractionType = FHTunnelingInteractionType::TEnearestNeigbors;
    fhLattice = lat;

    parameterSet1 = std::vector<double>();
    parameterSet2 = std::vector<double>();

    parameterSet1.push_back(tunnelingE);
    parameterSet2.push_back(neighbourD);


    uint nSites = fhLattice->getnSites();
    tunnelingInteraction = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    addTunneling(neighbourD, tunnelingE);


}


FHTunnelingInteraction::FHTunnelingInteraction(FHLattice* lat, std::vector<double> tunnelingDs, std::vector<double> tunnelingEs)
{
    tunnelingInteractionType = FHTunnelingInteractionType::TEdistanceInteraction0;
    fhLattice = lat;

    parameterSet1 = std::vector<double>(tunnelingEs);
    parameterSet2 = std::vector<double>(tunnelingDs);

    uint nSites = fhLattice->getnSites();
    tunnelingInteraction = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    uint nParam = std::max(parameterSet1.size(), parameterSet2.size());

    for (uint iDist = 0; iDist < nParam; iDist++)
    {
        addTunneling(parameterSet2[iDist], parameterSet1[iDist]);
    }

}



arma::cx_dmat*          FHTunnelingInteraction::getTunnelingInteractionPointer()
{
    return &tunnelingInteraction;
}

arma::cx_double         FHTunnelingInteraction::getTunnelingInteraction(uint iSite, uint jSite)
{
    if (tunnelingInteractionType == FHTunnelingInteractionType::TEnull)
    {
        return 0;
    }
    else 
    {
        return tunnelingInteraction(iSite, jSite);
    }

}




std::vector<double>         FHTunnelingInteraction::getParameters(uint set)
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

FHTunnelingInteractionType FHTunnelingInteraction::getTunnelingInteractionType()
{
    return tunnelingInteractionType;
}

Json::Value         FHTunnelingInteraction::writeTunnelingInteractionJson()
{
    Json::Value jpot;
    switch(tunnelingInteractionType)
    {
        case(FHTunnelingInteractionType::TEnull):
        {
            jpot["type"] = "null";
            break;
        }
        case(FHTunnelingInteractionType::TEnearestNeigbors):
        {
            jpot["type"] = "nearestNeighbors";
            jpot["neighborDistance"] = parameterSet2[0];
            jpot["neighborInteraction"] = parameterSet1[0];
            break;
        }
        case(FHTunnelingInteractionType::TEdistanceInteraction0):
        {
            jpot["type"] = "distanceInteraction";

            Json::Value dVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet2.size(); i++)
            {
                dVecJ.append(parameterSet2[i]);
            }
            jpot["distances"] = dVecJ;

            Json::Value eVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet1.size(); i++)
            {
                eVecJ.append(parameterSet1[i]);
            }
            jpot["energies"] = eVecJ;
            break;
        }
        default:
        {
            jpot["type"] = "null";
            break;
        }
    }

    return jpot;
}




void                        FHTunnelingInteraction::addTunneling(double distance, double Interaction)
{
    double ths = fhLattice->getSamePointDistance();

    uint nSites = fhLattice->getnSites();

    for (uint iSite=0; iSite<nSites; iSite++)
    {
        for (uint jSite = iSite+1; jSite<nSites; jSite++)
        {
            if (std::norm(fhLattice->getDistance(iSite, jSite) - distance) < ths)
            {
                tunnelingInteraction(iSite, jSite) += Interaction;
                tunnelingInteraction(jSite, iSite) += Interaction;
            }
        }
    }


}


FHTunnelingInteraction::~FHTunnelingInteraction()
{
    
}