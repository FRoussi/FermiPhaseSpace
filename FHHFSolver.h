#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

#include "FHSystem.h"
#include "ResultWriter.h"




class FHHFSolver
{
public:
    FHHFSolver(FHSystem* fhs);

    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoints = 1);
    void initiateSolver();

    void step_RK1(double dt);
    void step_RK2(double dt);
    void step_RK4(double dt);
    void step_RK45(double dt);
    void runSimulation();

    const arma::cx_double jj = arma::cx_double(0,1);

    ~FHHFSolver();

private:



// Fermi Hubbard System  :
    FHSystem* fhSystem;


// Time parameter
    double tMax;
    double dt;
    unsigned int nStepsPerTimePoint;
    unsigned int nTimePoints;

// Initial conditions :
    std::vector<bool> posPartU;
    std::vector<bool> posPartD;

// Link to Energy matrices
    arma::cx_dmat* U;
    arma::cx_dmat* J;

// Initial density matrices :
    arma::cx_dmat Nu0;
    arma::cx_dmat Nd0;

// Current density matrices :
    arma::cx_dmat Nut;
    arma::cx_dmat Ndt;


void derivative(arma::cx_dmat* nut, arma::cx_dmat* ndt, arma::cx_dmat* dnut_dt, arma::cx_dmat* dndt_dt);

};