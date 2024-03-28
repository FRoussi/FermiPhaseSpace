#include "FHExternalPotentialPulse.h"


FHExternalPotentialPulse::FHExternalPotentialPulse(const arma::dmat* sitePos, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar): FHExternalPotential(sitePos)
{
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    phi = 0;
    polarization = polar;
    withTunnelingCoefficients = false;
}


FHExternalPotentialPulse::FHExternalPotentialPulse(const arma::dmat* sitePos, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar, const arma::cx_dmat* tunMat): FHExternalPotential(sitePos)
{
    A0 = intensityA0;
    t0 = centert0;
    td = widthtd;
    w0 = frequencyw0;
    polarization = polar;
    phi = 0;
    withTunnelingCoefficients = false;

    setTunnelingMatrix(tunMat);

}


void FHExternalPotentialPulse::setTunnelingMatrix(const arma::cx_dmat* tunMat)
{
    withTunnelingCoefficients = true;
    tunnelingMatrix = tunMat;
}

void FHExternalPotentialPulse::updatePotential(double t)
{
    arma::cx_dvec A_t = arma::cx_double(A0 * std::exp( - pow((t-t0)/td, 2)/2 ) * std::cos(w0 * (t-t0) + phi),0 )* polarization;

    std::complex<double> zeroComp = std::complex<double>(0,0);

    for (uint iSite; iSite < sitePositions->n_rows ; iSite++)
    {
        for (uint jSite; jSite < sitePositions->n_rows ; jSite++)
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