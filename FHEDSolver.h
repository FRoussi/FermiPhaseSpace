#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
//#include <cmath>

#include "FHSystem.h"
#include "HermitianOperator.h"
#include "ResultWriter.h"

#define M_PI 3.14159265358979323846

std::vector<std::vector<bool>> getCombinations(unsigned int n, unsigned int p, unsigned int max_n);
unsigned int nChoosep(unsigned int n, unsigned int p);
std::vector<bool> positionToBinary(std::vector<unsigned int>, unsigned int nSites);
std::vector<unsigned int> binaryToPosition(std::vector<bool>);

unsigned int countSimilarities(std::vector<bool> phiH1, std::vector<bool> phiH2, unsigned int * lastIndexWrong1, unsigned int * lastIndexWrong2);
bool canHopHalf(std::vector<bool> phiH1, std::vector<bool> phiH2, unsigned int n_pH, unsigned int * idFrom, unsigned int * idTo);
unsigned int getParticleNumberBetween(std::vector<bool> posP, unsigned int from, unsigned int to);
arma::cx_dvec expmv(std::complex<double> t, arma::sp_cx_dmat A, arma::cx_dvec v, double tol = 1.e-7, unsigned int m = 30);

class FHEDSolver  	
{
public:
    FHEDSolver(FHSystem* fhSystem);

    void setTimeParameters(double t_max, double dt, unsigned int nPerTS);

    void initiateSolver();

    void buildHamiltonian();
    //void numberOperator();

    void step();
    void runSimulation();

    ~FHEDSolver();

private:
// Fermi Hubbard System  :
    FHSystem* fhSystem;

// Time parameter
    double tMax;
    double dt;
    unsigned int nStepsPerTimePoint;
    unsigned int nTimePoints;
// Other parameters :
    unsigned int nPartU;
    unsigned int nPartD;
    unsigned int basisSizeSpinUp;
    unsigned int basisSizeSpinDown;
    unsigned int basisSize;
    std::vector<std::vector<unsigned int>> basisArray;
    std::vector<std::vector<bool>> c_pu;
    std::vector<std::vector<bool>> c_pd;
// Initial conditions :
    std::vector<bool> posPartU;
    std::vector<bool> posPartD;
    arma::cx_dvec psi_0;

// Current Wave function :
    arma::cx_dvec psi;
// operators :
    arma::sp_cx_dmat H_U;
    arma::sp_cx_dmat H_J;
    arma::sp_cx_dmat H;

    arma::cx_dmat Hdense;
    //std::vector<arma::sp_cx_dmat> NiiuSet;
    //std::vector<arma::sp_cx_dmat> NiidSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijuSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijdSet;

    unsigned int getBasisIndex(std::vector<bool> posPU, std::vector<bool> posPD );


    std::complex<double> U_os_psil(unsigned int phi_l);


    void buildHamiltonian_U();
    void buildHamiltonian_J();
    void buildPairOperators();
    arma::sp_cx_dmat numberOperatorSpinUp(unsigned int iSite);
    arma::sp_cx_dmat numberOperatorSpinDown(unsigned int iSite);
    arma::sp_cx_dmat tunnelingOperatorSpinUp(unsigned int iSite, unsigned int jSite);
    arma::sp_cx_dmat tunnelingOperatorSpinDown(unsigned int iSite, unsigned int jSite);

};

