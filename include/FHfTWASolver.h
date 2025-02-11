#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <omp.h>




#include "FHSolver.h"



class FHfTWASolver : public FHSolver
{
public:

    enum fTWAIntegrator
    {
        RungeKutta2,
        RungeKutta4_38,
    };

    FHfTWASolver(FHSystem* fhs);


    void                        setIntegratorOrder(int order);
    void                        setfTWAParameters(unsigned int nTraj);
    void                        setMethodParameters(Json::Value param);



    void                        initiateSolver();

    void                        step_RK2(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double d_t);
    void                        step_RK4_38(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double d_t);
    void                        runSimulation();


    ~FHfTWASolver();

protected:




// fTWA parameters :    
    unsigned int                nTraj;

    // Integrator   
    int                         integratorOrder;
    fTWAIntegrator              integrator;



// Initial conditions :
    std::vector<bool>           posPartU;
    std::vector<bool>           posPartD;

    unsigned int nPartU;
    unsigned int nPartD;

    arma::dvec                  diagNu0;
    arma::dvec                  diagNd0;
    arma::dvec                  diagMu0;
    arma::dvec                  diagMd0;

    arma::cx_dmat               rhou0_mean;           // Initial rho matrices before adding the noise:
    arma::cx_dmat               rhod0_mean;

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

    void derivative(double currentTime, const arma::cx_dmat& nu, const arma::cx_dmat& nd, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);


};