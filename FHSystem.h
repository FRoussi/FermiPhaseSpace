#pragma once

#include <vector>
#include <armadillo>

std::complex<double> U_os_nij(double U0, double alpha, double d);

class FHSystem
{
public:
    FHSystem();
    FHSystem(int x);
    FHSystem(int x, int y, int z);
    FHSystem(std::string fileName);

    void setInteraction_OnSite(double U);
    void setTunneling_BetweenNeighbours(double J);

    unsigned int getnSites();
    double getCoordinate(int dim, int pos);

   double getDistance(int i, int j);
    bool areNeighbour(int i, int j);
    arma::cx_double getInteractionEnergy(int i, int j);
    arma::cx_double getTunnelingEnergy(int i, int j);
    arma::cx_dmat* getUMatrixPointer();
    arma::cx_dmat* getJMatrixPointer();

    ~FHSystem();

private:
    unsigned int nDimentions;
    unsigned int nSites;
    
    arma::dmat sitePositions;

    arma::dmat distanceMatrix;
    std::vector<std::vector<bool>> neighbourMatrix;

    arma::cx_dmat UMatrix;
    arma::cx_dmat JMatrix;

   void computeMatrices();


};