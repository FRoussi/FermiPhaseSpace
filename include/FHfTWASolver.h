#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>
#include <omp.h>

#include "FHSystem.h"
#include "ResultWriter.h"



class FHfTWASolver
{
public:

    enum fTWAIntegrator
    {
        RungeKutta2,
        RungeKutta4_38,
    };

    FHfTWASolver(FHSystem* fhs);

    void                        setResultFile(std::string workName, std::string folderNameName);
    void                        setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void                        setfTWAParameters(unsigned int nTraj);
    void                        setTalkative(bool yn);
    void                        initiateSolver();

    void                        step_RK2(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt);
    void                        step_RK4_38(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt);
    void                        runSimulation();

    const arma::cx_double       jj = arma::cx_double(0,1);

    void setComputeMomenta(bool yn);

    ~FHfTWASolver();

private:

// Fermi Hubbard System  :
    FHSystem*                   fhSystem;

// results in :         
    bool                        writeResults;
    std::string                 workName;              //<- A name for the work, for the results file
    std::string                 folderName;            //<- A name for the work folder, for the results file

    bool                        computeMomenta;

// Work parameters
    bool talkative;

// fTWA parameters :    
    unsigned int                nTraj;

    // Integrator   
    fTWAIntegrator              integrator;

// Time parameter
    double                      tMax;
    double                      dt;
    unsigned int                nStepsPerTimePoint;
    unsigned int                nTimePoints;

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

    arma::dmat                  V_1p;               //<- rotation matrix for the momentum space
    arma::dmat                  V_1p_1;             //<- Its transpose

    arma::cx_dmat               halfDiag;           //<- Diagonal matrix of 0.5
    arma::cx_dvec               halfVec;            //<- Vector of 0.5

    void derivative(double currentTime, const arma::cx_dmat& nu, const arma::cx_dmat& nd, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);


};