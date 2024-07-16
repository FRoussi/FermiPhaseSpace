#pragma once

#include <bitset>
//#include <climits>
#include <iostream>
#include <string>
#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
//#include <cmath>


#include "FHEDSolver.h"
#include "FHSystem.h"
#include "ResultWriterTxt.h"




std::vector<std::vector<bool>> getCombinations(unsigned int n, unsigned int p, unsigned int max_n);
unsigned int nChoosep(unsigned int n, unsigned int p);



bool canHopHalf(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, unsigned int n_pH, unsigned int& idFrom, unsigned int& idTo);


class tJEDSolver  	
{
public:
    tJEDSolver(FHSystem* fhSystem);

    void setResultFile(std::string workName, std::string folderNameName);
    void setTimeParameters(double t_max, double dt, unsigned int nPerTS);

    bool initiateSolver();

    void buildHamiltonian();
    //void numberOperator();

    void step(double d_t);
    void runSimulation();

    ~tJEDSolver();

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
    std::vector<std::array<std::bitset<MAX_BITS>,2>> basisArray;              
    std::vector<std::bitset<MAX_BITS>> c_pu;
    std::vector<std::bitset<MAX_BITS>> c_pd_reduced;
// Initial conditions :
    std::bitset<MAX_BITS> posPartU;
    std::bitset<MAX_BITS> posPartD;
    arma::cx_dvec psi_0;                                            //<- Initial wave function

// Current Wave function :
    arma::cx_dvec psi;                                              //<- Current wave function
// operators :
    arma::sp_cx_dmat H_U;
    arma::sp_cx_dmat H_t;
    arma::sp_cx_dmat H;

    //arma::cx_dmat Hdense;
    //std::vector<arma::sp_cx_dmat> NiiuSet;
    //std::vector<arma::sp_cx_dmat> NiidSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijuSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijdSet;

    unsigned int getBasisIndex(std::bitset<MAX_BITS>& posPU, std::bitset<MAX_BITS>& posPD );          //<- Get the basis index of the basis vector corresponding to particle positions 'posPU' and 'posPD'


    std::complex<double> U_os_psil(unsigned int phi_l);                                     //<- Compute the interation energy of the basis vector 'phi_l'


    void buildHamiltonian_U();                                                              //<- Build and store the operator "Interaction part of the Hamiltonian"
    void buildHamiltonian_tM();                                                              //<- Build and store the operator "Tunneling part of the Hamiltonian"
    void buildPairOperators();                                                              //<- Build and store all the usefull number and tunneling operators.
    arma::sp_cx_dmat numberOperatorSpinUp(unsigned int iSite);                              //<- Build the operator "Number of up-spin particles on site 'iSite'".
    arma::sp_cx_dmat numberOperatorSpinDown(unsigned int iSite);                            //<- Build the operator "Number of down-spin particles on site 'iSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinUp(unsigned int iSite, unsigned int jSite);       //<- Build the operator "Tunneling of up-spin particles between sites 'iSite' and 'jSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinDown(unsigned int iSite, unsigned int jSite);     //<- Build the operator "Tunneling of down-spin particles between sites 'iSite' and 'jSite'".
    arma::sp_cx_dmat exchangeOperatorSpinUp(unsigned int iSite, unsigned int jSite);       //<- Build the operator "Exchange of up-spin and spin-down particles between sites 'iSite' and 'jSite'".

};

