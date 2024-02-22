#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
//#include <cmath>

#include "FHSystem.h"


#include "ResultWriter.h"


#define M_PI 3.14159265358979323846

std::vector<std::vector<bool>> getCombinations(unsigned int n, unsigned int p, unsigned int max_n);
unsigned int nChoosep(unsigned int n, unsigned int p);
std::vector<bool> positionToBinary(std::vector<unsigned int>&, unsigned int nSites);
std::vector<unsigned int> binaryToPosition(std::vector<bool>&);

unsigned int countSimilarities(std::vector<bool>& phiH1, std::vector<bool>& phiH2, unsigned int& lastIndexWrong1, unsigned int& lastIndexWrong2);
bool canHopHalf(std::vector<bool>& phiH1, std::vector<bool>& phiH2, unsigned int n_pH, unsigned int& idFrom, unsigned int& idTo);
unsigned int getParticleNumberBetween(std::vector<bool>& posP, unsigned int from, unsigned int to);
arma::cx_dvec expmv(std::complex<double> t, arma::sp_cx_dmat& A, arma::cx_dvec& v, double tol = 1.e-7, unsigned int m = 30);

class FHEDSolver  	
{
public:
    FHEDSolver(FHSystem* fhSystem);

    void setResultFile(std::string workName, std::string folderNameName);
    void setTimeParameters(double t_max, double dt, unsigned int nPerTS);

    void initiateSolver();

    void buildHamiltonian();
    //void numberOperator();

    void step(double d_t);
    void runSimulation();

    ~FHEDSolver();

private:
// Fermi Hubbard System  :
    FHSystem* fhSystem;                                             //<- Pointer to the Fermi-Hubbard model

// results in :
    bool writeResults;
    std::string workName;                                           //<- A name for the work, for the results file
    std::string folderName;                                         //<- A name for the work folder, for the results file

// Time parameter
    double tMax;                                                    //<- Maximum time extend of the simulation 
    double dt;                                                      //<- Time step of calculated observables        
    unsigned int nStepsPerTimePoint;                                //<- Number of simulated time steps between the calculation of observable
    unsigned int nTimePoints;                                       //<- Number of time steps (=tMax/dt)
// Other parameters :
    unsigned int nPartU;                                            //<- Number of up-spin fermions
    unsigned int nPartD;                                            //<- Number of down-spin fermions
    unsigned int basisSizeSpinUp;                                   //<- Size of the basis of up-spin fermions (nPartU among nSites)
    unsigned int basisSizeSpinDown;                                 //<- Size of the basis of down-spin fermions (nPartD among nSites)
    unsigned int basisSize;                                         //<- Size of the Hilbert space basis  (basisSizeSpinUp X basisSizeSpinDown)
    std::vector<std::vector<unsigned int>> basisArray;              
    std::vector<std::vector<bool>> c_pu;
    std::vector<std::vector<bool>> c_pd;
// Initial conditions :
    std::vector<bool> posPartU;
    std::vector<bool> posPartD;
    arma::cx_dvec psi_0;                                            //<- Initial wave function

// Current Wave function :
    arma::cx_dvec psi;                                              //<- Current wave function
// operators :
    arma::sp_cx_dmat H_U;
    arma::sp_cx_dmat H_J;
    arma::sp_cx_dmat H;

    //arma::cx_dmat Hdense;
    //std::vector<arma::sp_cx_dmat> NiiuSet;
    //std::vector<arma::sp_cx_dmat> NiidSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijuSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijdSet;

    unsigned int getBasisIndex(std::vector<bool>& posPU, std::vector<bool>& posPD );          //<- Get the basis index of the basis vector corresponding to particle positions 'posPU' and 'posPD'


    std::complex<double> U_os_psil(unsigned int phi_l);                                     //<- Compute the interation energy of the basis vector 'phi_l'


    void buildHamiltonian_U();                                                              //<- Build and store the operator "Interaction part of the Hamiltonian"
    void buildHamiltonian_J();                                                              //<- Build and store the operator "Tunneling part of the Hamiltonian"
    void buildPairOperators();                                                              //<- Build and store all the usefull number and tunneling operators.
    arma::sp_cx_dmat numberOperatorSpinUp(unsigned int iSite);                              //<- Build the operator "Number of up-spin particles on site 'iSite'".
    arma::sp_cx_dmat numberOperatorSpinDown(unsigned int iSite);                            //<- Build the operator "Number of down-spin particles on site 'iSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinUp(unsigned int iSite, unsigned int jSite);       //<- Build the operator "Tunneling of up-spin particles between sites 'iSite' and 'jSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinDown(unsigned int iSite, unsigned int jSite);     //<- Build the operator "Tunneling of down-spin particles between sites 'iSite' and 'jSite'".

};

