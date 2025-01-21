#include "../include/FHExternalPotentialConstant.h"


FHExternalPotentialConstant::FHExternalPotentialConstant(FHLattice* lat, arma::dvec Pot): FHExternalPotential(lat)
{
    externalPotential.diag() = arma::cx_dmat(Pot, arma::dvec(Pot.size(), arma::fill::zeros));
}


FHExternalPotentialConstant::FHExternalPotentialConstant(FHLattice* lat, uint iSite, double pot): FHExternalPotential(lat)
{
    externalPotential(iSite, iSite) = std::complex<double>(pot);
}


FHLattice::~FHLattice()
{

}
