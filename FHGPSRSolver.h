#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <armadillo>
#include <omp.h>

#include "FHSystem.h"



class FHGPSRSolver
{
public:
    FHGPSRSolver(FHSystem* fhs);


    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void setGPSRParameters(unsigned int nTraj);
    void initiateSolver(std::vector<unsigned int> posPartU, std::vector<unsigned int> posPartD );

    void step_Drift(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Drift, arma::cx_dmat* drhod_Drift);
    void step_Diffusion(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Diff, arma::cx_dmat* drhod_Diff);
    void step(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt);
    void runSimulation();

    const arma::cx_double jj = arma::cx_double(0,1);
    const arma::cx_double sqrt_mjj = std::sqrt(-jj);
    const arma::cx_double sqrt_pjj = std::sqrt(jj);

    

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

  void derivative(arma::cx_dmat* nut, arma::cx_dmat* ndt, arma::cx_dmat* dnut_dt, arma::cx_dmat* dndt_dt);



};