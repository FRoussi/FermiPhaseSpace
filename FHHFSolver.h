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

    enum HFIntegrator
    {
        RungeKutta1,
        RungeKutta2,
        RungeKutta4,
        RungeKutta4_38,
        RungeKutta45
    };

    FHHFSolver(FHSystem* fhs);

    void setResultFile(std::string workName, std::string folderNameName);
    void setTimeParameters(double t_max, double dt, unsigned int nStepsPerTimePoints = 1);
    void initiateSolver();

    void step_RK1(double dt);
    void step_RK2(double dt);
    void step_RK4(double dt);
    void step_RK4_38(double dt);
    void step_RK45(double dt);
    void runSimulation();

    const arma::cx_double jj = arma::cx_double(0,1);

    ~FHHFSolver();

private:


    FHSystem*               fhSystem;           //<- Fermi Hubbard System  

// results in :
    bool                    writeResults;
    std::string             workName;           //<- A name for the work, for the results file
    std::string             folderName;         //<- A name for the work folder, for the results file

// Time parameter
    double                  tMax;
    double                  dt;
    unsigned int            nStepsPerTimePoint;
    unsigned int            nTimePoints;


    HFIntegrator            integrator;            //<-Integrator 


    std::vector<bool>       posPartU;              // Initial conditions vector
    std::vector<bool>       posPartD;
    unsigned int            nPartU;                 //<- Number of up-spin particles
    unsigned int            nPartD;                 //<- Number of down-spin particles

    arma::cx_dmat*          U;                     // Link to Energy matrices
    arma::cx_dmat*          tM;

    arma::cx_dmat           Nu0;                  //<- Initial density matrices 
    arma::cx_dmat           Nd0;

    arma::cx_dmat           Nut;                  //<- Current density matrices
    arma::cx_dmat           Ndt;

    arma::dmat              V_1p;                   //<- rotation matrix for the momentum space
    arma::dmat              V_1p_1;                 //<- Its transpose

    void derivative(const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);

};