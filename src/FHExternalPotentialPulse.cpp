#include "../include/FHExternalPotentialPulse.h"


FHExternalPotentialPulse::FHExternalPotentialPulse(FHLattice* lat, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    phi = 0;
    polarization = polar;
    withTunnelingCoefficients = false;


    if (A0 !=0)
    {
        potentialType = FHExternalPotentialType::EPtimeDependant;
    }

}


FHExternalPotentialPulse::FHExternalPotentialPulse(FHLattice* lat, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar, const arma::cx_dmat& tunMat): FHExternalPotential(lat)
{
    uint nSites = lat->getnSites();
    externalPotential = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    polarization = polar;
    phi = 0;
    withTunnelingCoefficients = false;

    if (A0 !=0)
    {
        potentialType = FHExternalPotentialType::EPtimeDependant;
    }

    setTunnelingMatrix(tunMat);

}


void FHExternalPotentialPulse::setTunnelingMatrix(const arma::cx_dmat& tunMat)
{
    withTunnelingCoefficients = true;
    tunnelingMatrix = &tunMat;
}


arma::cx_dmat FHExternalPotentialPulse::getExternalPotential(double t)
{

    arma::cx_dmat extPotential = externalPotential;
    arma::cx_dvec A_t = arma::cx_double(A0 * std::exp( - pow((t-t0)/td, 2)/2 ) * std::cos(w0 * (t-t0) + phi),0 )* polarization;

    std::complex<double> zeroComp = std::complex<double>(0,0);


    for (uint iSite = 0; iSite < fhLattice->getnSites() ; iSite++)
    {
        for (uint jSite = 0; jSite < fhLattice->getnSites() ; jSite++)
        {
            if (withTunnelingCoefficients)
            {
                if (tunnelingMatrix->at(iSite, jSite) != zeroComp)
                {
                    extPotential(iSite,jSite) = tunnelingMatrix->at(iSite, jSite) * std::exp(jj*arma::dot(A_t,(fhLattice->getSitePosition(iSite) - fhLattice->getSitePosition(jSite))));  
                }
                else
                {
                    extPotential(iSite,jSite) = zeroComp;
                }

            }
        }
    }

    return extPotential;
}


void FHExternalPotentialPulse::updatePotential(double t)
{
    arma::cx_dvec A_t = arma::cx_double(A0 * std::exp( - pow((t-t0)/td, 2)/2 ) * std::cos(w0 * (t-t0) + phi),0 )* polarization;

    std::complex<double> zeroComp = std::complex<double>(0,0);

    for (uint iSite = 0; iSite < fhLattice->getnSites() ; iSite++)
    {
        for (uint jSite = 0; jSite < fhLattice->getnSites() ; jSite++)
        {
            if (withTunnelingCoefficients)
            {
                if (tunnelingMatrix->at(iSite, jSite) != zeroComp)
                {
                    externalPotential(iSite,jSite) = tunnelingMatrix->at(iSite, jSite) * std::exp(jj*arma::dot(A_t,(fhLattice->getSitePosition(iSite) - fhLattice->getSitePosition(jSite))));  
                }
                else
                {
                    externalPotential(iSite,jSite) = zeroComp;
                }

            }
        }
    }
}

Json::Value                 FHExternalPotentialPulse::writePotentialJson()
{
    Json::Value jpot;
    jpot["type"] = "pulse";
    jpot["A0"] = A0;
    jpot["t0"] = t0;
    jpot["td"] = td;
    jpot["w0"] = w0;
    jpot["phi"] = phi;
    return jpot;
}

FHExternalPotentialPulse::~FHExternalPotentialPulse()
{
    
}