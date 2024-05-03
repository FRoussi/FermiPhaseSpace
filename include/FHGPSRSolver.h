#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>
#include <omp.h>

#include "FHSystem.h"
#include "ResultWriter.h"



bool rsvd(arma::cx_dmat& U, arma::dvec& s, arma::cx_dmat& V, arma::cx_dmat& D, unsigned int k);


class FHGPSRSolver
{
public:

    enum GPSRDriftIntegrator
    {
        RungeKutta2,
        RungeKutta4_38
    };
    enum GPSRDiffusionIntegrator
    {
        analytical,
        numerical_svd,
        numerical_rsvd
    };

    FHGPSRSolver(FHSystem* fhs);

    void                setResultFile(std::string workName, std::string folderNameName);
    void                setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void                setGPSRParameters(unsigned int nTraj);
    void                initiateSolver();

    void                step_Drift_RK2(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double dt, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift);
    void                step_Drift_RK4_38(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double dt, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift);

    void                step_Diffusion(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double dt, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff);
    void                step_Diffusion_svd(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double dt, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff);

    arma::cx_dmat       diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0);
    arma::cx_dmat       quarterDiffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0);
    arma::cx_dmat       diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dvec UMat);
    arma::cx_dmat       diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat UMat);

    arma::cx_dmat       sqrtDecomposition(arma::cx_dmat& D);
    arma::cx_dmat       sqrtDecomposition_fromA(arma::cx_dmat& A);
    
    void                step(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt);
    void                runSimulation();

    const std::complex<double> jj       = std::complex<double>(0,1);
    const std::complex<double> sqrt_pjj = std::sqrt(jj);
    const std::complex<double> sqrt_mjj = std::sqrt(-jj);

    ~FHGPSRSolver();

private:

// Fermi Hubbard System  :
    FHSystem*                       fhSystem;

// results in :
    bool                            writeResults;
    std::string                     workName;                                           //<- A name for the work, for the results file
    std::string                     folderName;                                         //<- A name for the work folder, for the results file

// GPSR parameters :
    unsigned int                    nTraj;

    GPSRDriftIntegrator             driftIntegrator;
    GPSRDiffusionIntegrator         diffusionIntegrator;

// Time parameter
    double                          tMax;
    double                          dt;
    unsigned int                    nStepsPerTimePoint;
    unsigned int                    nTimePoints;

// Initial conditions :
    std::vector<bool>               posPartU;
    std::vector<bool>               posPartD;

    arma::dvec                      diagNu0;
    arma::dvec                      diagNd0;

// Initial rho matrices :
    arma::cx_dmat                   rhou0;
    arma::cx_dmat                   rhod0;

// Link to Energy matrices
    const arma::cx_dmat*            U;                     
    const arma::cx_dmat*            tM;
    const arma::cx_dmat*            J;
    const arma::cx_dmat*            V;
    const arma::cx_dmat*            extPot;

// Usefull things :
    arma::cx_dmat                   sqrtU;

    void derivative(double currentSimulationTime, const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);



};