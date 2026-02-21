#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <omp.h>




#include "FHSystem.h"



class FHfTWAStateSearch
{
public:

    enum fTWAStateSearchMethods
    {
        MC1,
    };

    FHfTWAStateSearch(FHSystem* fhs);




    void                        findDensity();
    void                        initiateStateSearch();




    ~FHfTWAStateSearch();

protected:


// Fermi Hubbard System  :
FHSystem*               fhSystem;                     //<- Pointer to the Fermi-Hubbard model
Json::Value             jsonFile;   

// fTWA parameters :    
    uint                nPoints;
    uint                nTries;

    uint                nMaxPoints;
    uint                nMaxTries;


// Initial conditions :
    std::vector<bool>           posPartU;
    std::vector<bool>           posPartD;

    unsigned int nPartU;
    unsigned int nPartD;

    arma::dvec                  diagNu0;
    arma::dvec                  diagNd0;
    arma::dvec                  diagMu0;
    arma::dvec                  diagMd0;

    arma::cx_dmat               rhou0;           // Density matrices used for the next step:
    arma::cx_dmat               rhod0;

    arma::cx_dmat               rhou_mean;           // Initial rho matrices before adding the noise:
    arma::cx_dmat               rhod_mean;

    arma::cx_dmat               rhou;           // Next Step:
    arma::cx_dmat               rhod;


    double                      maxE;
    double                      meanE;
    double                      minE;

    double                      kT;


// Initial noise Variance matrix :
    arma::dmat                  sqrtNoiseVarianceMatrixU;
    arma::dmat                  sqrtNoiseVarianceMatrixD;

    const arma::cx_dmat*        U;                     // Link to Energy matrices
    const arma::cx_dmat*        tM;
    const arma::cx_dmat*        J;
    const arma::cx_dmat*        V;
    const arma::cx_dmat*        extPot;

   // arma::dmat                  V_1p;               //<- rotation matrix for the momentum space
    //arma::dmat                  V_1p_1;             //<- Its transpose

    arma::cx_dmat               halfDiag;           //<- Diagonal matrix of 0.5
    arma::cx_dvec               halfVec;            //<- Vector of 0.5

    void derivative(double currentTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);
    void energy(double currentTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_double& EtH, arma::cx_double& EU, arma::cx_double& Epe);

    void makeTry( const arma::cx_dmat& rhou0, const arma::cx_dmat& rhod0,  arma::cx_dmat& rhou, arma::cx_dmat& rhod);


};