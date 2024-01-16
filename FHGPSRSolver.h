#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>
#include <omp.h>

#include "FHSystem.h"
#include "ResultWriter.h"



class FHGPSRSolver
{
public:
    FHGPSRSolver(FHSystem* fhs);


    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void setGPSRParameters(unsigned int nTraj);
    void initiateSolver();

    void step_Drift_RK2(const arma::cx_dmat * rhou, const arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Drift, arma::cx_dmat* drhod_Drift);
    void step_Drift_RK4_38(const arma::cx_dmat * rhou, const arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Drift, arma::cx_dmat* drhod_Drift);
    void step_Diffusion(const arma::cx_dmat * rhou, const arma::cx_dmat* rhod, double U0, double dt, arma::cx_dmat * drhou_Diff, arma::cx_dmat* drhod_Diff);
    void step_Diffusion_svd(const arma::cx_dmat * rhou, const arma::cx_dmat* rhod, double U0, double dt, arma::cx_dmat * drhou_Diff, arma::cx_dmat* drhod_Diff);
    arma::cx_dmat diffusionMatrix(const arma::cx_dmat * rhou, const arma::cx_dmat* rhod, double U0);
    arma::cx_dmat sqrtDecomposition(arma::cx_dmat D);
    void step(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt);
    void runSimulation();

    const std::complex<double> jj = arma::cx_double(0,1);
    const std::complex<double> sqrt_mjj = std::sqrt(-jj);
    const std::complex<double> sqrt_pjj = std::sqrt(jj);

    

    ~FHGPSRSolver();

private:



// Fermi Hubbard System  :
    FHSystem* fhSystem;

// GPSR parameters :
    unsigned int nTraj;

// Time parameter
    double tMax;
    double dt;
    unsigned int nStepsPerTimePoint;
    unsigned int nTimePoints;

// Initial conditions :
    std::vector<bool> posPartU;
    std::vector<bool> posPartD;

    arma::dvec diagNu0;
    arma::dvec diagNd0;

// Initial rho matrices :
    arma::cx_dmat rhou0;
    arma::cx_dmat rhod0;



// Link to Energy matrices
  arma::cx_dmat* U;
  arma::cx_dmat* J;

// Usefull things :
  arma::cx_dmat sqrtU;

  void derivative(const arma::cx_dmat* nut, const arma::cx_dmat* ndt, arma::cx_dmat* dnut_dt, arma::cx_dmat* dndt_dt);



};