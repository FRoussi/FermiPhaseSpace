#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <omp.h>

//#define ARMA_DONT_USE_OPENMP
#include <armadillo>



#include "FHSolver.h"



class FHHFSolver : public FHSolver
{
public:

    enum HFIntegrator
    {
        RungeKutta1,
        RungeKutta1_Mem,
        RungeKutta2,
        RungeKutta4,
        RungeKutta4_38,
        RungeKutta45
    };

    FHHFSolver(FHSystem* fhs);

    void                        initiateSolver();
    void                        setComputeMomenta(bool yn);


    void                        step_RK1(double dt);
    void                        step_RK1_Mem(double dt, arma::cx_dmat & dNu_dt, arma::cx_dmat & dNd_dt);
    void                        step_RK2(double dt);
    void                        step_RK4(double dt);
    void                        step_RK4_38(double dt);
    void                        step_RK45(double dt);
    void                        runSimulation();


    ~FHHFSolver();

private:
        

// Time parameter
    double                  currentSimulationTime;

    HFIntegrator            integrator;            //<-Integrator 


    std::vector<bool>       posPartU;              // Initial conditions vector
    std::vector<bool>       posPartD;
    unsigned int            nPartU;                 //<- Number of up-spin particles
    unsigned int            nPartD;                 //<- Number of down-spin particles

    const arma::cx_dmat*    U;                     // Link to Energy matrices
    const arma::cx_dmat*    tM;
    const arma::cx_dmat*    J;
    const arma::cx_dmat*    V;
    const arma::cx_dmat*    extPot;

    arma::cx_dmat           Nu0;                  //<- Initial density matrices 
    arma::cx_dmat           Nd0;

    arma::cx_dmat           Nut;                  //<- Current density matrices
    arma::cx_dmat           Ndt;

    arma::dmat              V_1p;                   //<- rotation matrix for the momentum space
    arma::dmat              V_1p_1;                 //<- Its transpose

    void derivative(double currentTime, const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt);

};