#pragma once

#include "FHExternalPotential.h"



class FHExternalPotentialPulse : public FHExternalPotential
{
public:

    FHExternalPotentialPulse(const arma::dmat* sitePos, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar);
    FHExternalPotentialPulse(const arma::dmat* sitePos, double intensityA0, double centert0, double widthtd, double frequencyw0, arma::dvec polar, const arma::cx_dmat* tunMat);

    void setTunnelingMatrix(const arma::cx_dmat* tunMat);

    void updatePotential(double t);

private :

    bool withTunnelingCoefficients;



    double A0;
    double t0;
    double td;
    double w0;
    double phi;
    arma::dvec polarization;

    std::complex<double> jj = std::complex<double>(0,1);

    const arma::cx_dmat* tunnelingMatrix;

};