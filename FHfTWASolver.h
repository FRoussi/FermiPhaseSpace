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

    void setResultFile(std::string workName, std::string folderNameName);
    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoint = 1);
    void setfTWAParameters(unsigned int nTraj);
    void initiateSolver();

    void step_RK2(arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt);
    void step_RK4_38(arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt);
    void runSimulation();

    const arma::cx_double jj = arma::cx_double(0,1);

    ~FHfTWASolver();

private:



// Fermi Hubbard System  :
    FHSystem* fhSystem;

// results in :
    bool writeResults;
    std::string workName;                                           //<- A name for the work, for the results file
    std::string folderName;                                         //<- A name for the work folder, for the results file

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

    void derivative(const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);


};