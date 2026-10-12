#include "../include/FHCoulombInteraction.h"


FHCoulombInteraction::FHCoulombInteraction()
{
    CoulombInteractionType = FHCoulombInteractionType::CEnull;
}

FHCoulombInteraction::FHCoulombInteraction(FHLattice* lat)
{
    CoulombInteractionType = FHCoulombInteractionType::CEnull;
    fhLattice = lat;
}

FHCoulombInteraction::FHCoulombInteraction(FHLattice* lat, Json::Value coulJson)
{
    CoulombInteractionType = FHCoulombInteractionType::CEnull;
    fhLattice = lat;
    std::string coulType = coulJson.get("type", "null").asString();
    if (coulType == "onSite")
    {
        CoulombInteractionType = FHCoulombInteractionType::CEonSite;

        double onSiteE = coulJson.get("onSiteInteraction", 1.).asDouble();

        parameterSet1 = std::vector<double>();
        parameterSet2 = std::vector<double>();

        parameterSet1.push_back(onSiteE);
        parameterSet2.push_back(0);

        uint nSites = fhLattice->getnSites();
        U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

        addCoulomb(0, onSiteE);
    }
}

FHCoulombInteraction::FHCoulombInteraction(FHLattice* lat, double onSiteInteraction)
{

    CoulombInteractionType = FHCoulombInteractionType::CEonSite;
    fhLattice = lat;

    parameterSet1 = std::vector<double>();
    parameterSet2 = std::vector<double>();

    parameterSet1.push_back(onSiteInteraction);
    parameterSet2.push_back(0);



    uint nSites = fhLattice->getnSites();
    U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    addCoulombinU(0, onSiteInteraction);


}


FHCoulombInteraction::FHCoulombInteraction(FHLattice* lat, std::vector<double> CoulombDs, std::vector<double> CoulombEs)
{
    CoulombInteractionType = FHCoulombInteractionType::CEdistanceInteraction;
    fhLattice = lat;

    parameterSet1 = std::vector<double>(CoulombEs);
    parameterSet2 = std::vector<double>(CoulombDs);

    uint nSites = fhLattice->getnSites();
    U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    uint nParam = std::max(parameterSet1.size(), parameterSet2.size());

    for (uint iDist = 0; iDist < nParam; iDist++)
    {
        addCoulomb(parameterSet2[iDist], parameterSet1[iDist]);
    }

}



arma::cx_dmat*          FHCoulombInteraction::getCoulombUInteractionPointer()
{
    return &U;
}

arma::cx_dmat*          FHCoulombInteraction::getCoulombVInteractionPointer()
{
    return &V;
}


arma::cx_double         FHCoulombInteraction::getCoulombUInteraction(uint iSite, uint jSite) const
{
    if (CoulombInteractionType == FHCoulombInteractionType::CEnull)
    {
        return 0;
    }
    else 
    {
        return U(iSite, jSite);
    }

}

arma::cx_double         FHCoulombInteraction::getCoulombVInteraction(uint iSite, uint jSite) const
{
    if (CoulombInteractionType == FHCoulombInteractionType::CEnull)
    {
        return 0;
    }
    else 
    {
        return V(iSite, jSite);
    }

}




std::vector<double>         FHCoulombInteraction::getParameters(uint set) const
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

FHCoulombInteractionType FHCoulombInteraction::getCoulombInteractionType() const
{
    return CoulombInteractionType;
}



Json::Value         FHCoulombInteraction::writeCoulombInteractionJson()
{
    Json::Value jpot;
    switch(CoulombInteractionType)
    {
        case(FHCoulombInteractionType::CEnull):
        {
            jpot["type"] = "null";
            break;
        }
        case(FHCoulombInteractionType::CEonSite):
        {
            jpot["type"] = "onSite";
            jpot["onSiteInteraction"] = parameterSet1[0];
            break;
        }
        case(FHCoulombInteractionType::CEdistanceInteraction):
        {
            jpot["type"] = "distanceInteraction";

            Json::Value dVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet1.size(); i++)
            {
                dVecJ.append(parameterSet2[i]);
            }
            jpot["distances"] = dVecJ;

            Json::Value eVecJ = Json::Value(Json::arrayValue);
            for (uint i = 0; i<parameterSet2.size(); i++)
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

void                         FHCoulombInteraction::addSite(arma::cx_dvec siteInteractions)
{
    uint nSites = fhLattice->getnSites() + 1;
    if (siteInteractions.size() != nSites)
    {
        std::cout << "!!error : FHCoulombInteraction - addSite - Wrong vector size" << std::endl;
    }
    else
    {
       
        if (CoulombInteractionType == FHCoulombInteractionType::CEnull)
        {

            U = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
            V = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        }
        else if (CoulombInteractionType == FHCoulombInteractionType::CEonSite)
        {

            V = arma::cx_dmat(nSites, nSites);
            arma::cx_dmat newU = arma::cx_dmat(nSites, nSites);
            newU.submat(0, 0, arma::size(U)) = U;
            U = newU;
        }
        else
        {

            arma::cx_dmat newU = arma::cx_dmat(nSites+1, nSites+1);
            newU.submat(0, 0, arma::size(U)) = U;
            U = newU;

            arma::cx_dmat newV = arma::cx_dmat(nSites+1, nSites+1);
            newV.submat(0, 0, arma::size(V)) = V;
            V = newV;
        }


        U.submat(0, nSites-1, nSites-1, nSites-1) += siteInteractions/arma::cx_double(2,0);
        U.submat(nSites-1, 0, nSites-1, nSites-1) += siteInteractions.st()/arma::cx_double(2,0);


        V.submat(0, nSites-1, nSites-1, nSites-1) += siteInteractions/arma::cx_double(2,0);
        V.submat(nSites-1, 0, nSites-1, nSites-1) += siteInteractions.st()/arma::cx_double(2,0);
        V.at(nSites-1,nSites-1) = 0;


        CoulombInteractionType = FHCoulombInteractionType::CEgiven;
    
    }
}

void                        FHCoulombInteraction::addCoulomb(double distance, double Interaction)
{
    addCoulombinU(distance, Interaction);
    addCoulombinV(distance, Interaction);
}

void                        FHCoulombInteraction::addCoulombinU(double distance, double Interaction)
{
    double ths = fhLattice->getSamePointDistance();

    uint nSites = fhLattice->getnSites();

    for (uint iSite=0; iSite<nSites; iSite++)
    {
        for (uint jSite = 0; jSite<nSites; jSite++)
        {
            if (std::norm(fhLattice->getDistance(iSite, jSite) - distance) < ths)
            {
                U.at(iSite, jSite) += Interaction;
            }
        }
    }


}




void                        FHCoulombInteraction::addCoulombinV(double distance, double Interaction, bool fillDiagonal)
{
    double ths = fhLattice->getSamePointDistance();

    uint nSites = fhLattice->getnSites();

    for (uint iSite=0; iSite<nSites; iSite++)
    {
        for (uint jSite = 0; jSite<nSites; jSite++)
        {
            if ((iSite != jSite) || fillDiagonal)
            {
                if (std::norm(fhLattice->getDistance(iSite, jSite) - distance) < ths)
                {
                    V.at(iSite, jSite) += Interaction;
                }
            }
        }
    }
}


FHCoulombInteraction::~FHCoulombInteraction()
{

}