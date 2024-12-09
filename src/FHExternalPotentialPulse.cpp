#include "../include/FHExternalPotentialPulse.h"


FHExternalPotentialPulse::FHExternalPotentialPulse(FHLattice* lat, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar): FHExternalPotential(lat)
{
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    phi = 0;
    polarization = polar;
    withTunnelingCoefficients = false;

    if (A0 !=0)
    {
        nullPotential = false;
        timeDependantPotential = true;
    }

}


FHExternalPotentialPulse::FHExternalPotentialPulse(FHLattice* lat, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar, const arma::cx_dmat* tunMat): FHExternalPotential(lat)
{
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    polarization = polar;
    phi = 0;
    withTunnelingCoefficients = false;

    if (A0 !=0)
    {
        nullPotential = false;
        timeDependantPotential = true;
    }

    setTunnelingMatrix(tunMat);

}


void FHExternalPotentialPulse::setTunnelingMatrix(const arma::cx_dmat* tunMat)
{
    withTunnelingCoefficients = true;
    tunnelingMatrix = tunMat;
}


arma::cx_dmat FHExternalPotentialPulse::getExternalPotential(double t)
{

    arma::cx_dmat extPotential = externalPotential;
    arma::cx_dvec A_t = arma::cx_double(A0 * std::exp( - pow((t-t0)/td, 2)/2 ) * std::cos(w0 * (t-t0) + phi),0 )* polarization;

    std::complex<double> zeroComp = std::complex<double>(0,0);

    for (uint iSite = 0; iSite < sitePositions->n_rows ; iSite++)
    {
        for (uint jSite = 0; jSite < sitePositions->n_rows ; jSite++)
        {
            if (withTunnelingCoefficients)
            {
                if (tunnelingMatrix->at(iSite, jSite) != zeroComp)
                {
                    extPotential(iSite,jSite) = tunnelingMatrix->at(iSite, jSite) * std::exp(jj*arma::dot(A_t,(sitePositions->row(iSite) - sitePositions->row(jSite))));  
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

    for (uint iSite = 0; iSite < sitePositions->n_rows ; iSite++)
    {
        for (uint jSite = 0; jSite < sitePositions->n_rows ; jSite++)
        {
            if (withTunnelingCoefficients)
            {
                if (tunnelingMatrix->at(iSite, jSite) != zeroComp)
                {
                    externalPotential(iSite,jSite) = tunnelingMatrix->at(iSite, jSite) * std::exp(jj*arma::dot(A_t,(sitePositions->row(iSite) - sitePositions->row(jSite))));  
                }
                else
                {
                    externalPotential(iSite,jSite) = zeroComp;
                }

            }
        }
    }
}