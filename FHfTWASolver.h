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
    FHfTWASolver(FHSystem* fhs);

    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void setfTWAParameters(unsigned int nTraj);
    void initiateSolver();

    void step(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt);
    void runSimulation();

    const arma::cx_double jj = arma::cx_double(0,1);

    ~FHfTWASolver();

private:



// Fermi Hubbard System  :
    FHSystem* fhSystem;

// fTWA parameters :
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

// Initial rho matrices before adding the noise:
    arma::cx_dmat rhou0_mean;
    arma::cx_dmat rhod0_mean;

// Initial noise Variance matrix :
    arma::dmat sqrtNoiseVarianceMatrixU;
    arma::dmat sqrtNoiseVarianceMatrixD;


// Link to Energy matrices
    arma::cx_dmat* U;
    arma::cx_dmat* J;

// Usefull things :
    arma::cx_dmat halfDiag;

    void derivative(arma::cx_dmat* nut, arma::cx_dmat* ndt, arma::cx_dmat* dnut_dt, arma::cx_dmat* dndt_dt);


};