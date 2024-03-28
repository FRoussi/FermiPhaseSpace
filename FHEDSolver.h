#pragma once

#include <bitset>
//#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
//#include <cmath>


#include "FHSystem.h"
#include "ResultWriter.h"


#define MAX_BITS 32
#define M_PI 3.14159265358979323846


std::vector<std::vector<bool>>      getCombinations(unsigned int n, unsigned int p, unsigned int max_n);
unsigned int                        nChoosep(unsigned int n, unsigned int p);

std::bitset<MAX_BITS>               bitsetFromVector(const std::vector<bool>& vec);

void                                bitsetswap(std::bitset<MAX_BITS>& bs1, std::bitset<MAX_BITS>& bs2);

void                                fusionSortVector(std::vector<std::bitset<MAX_BITS>>& vec);

std::bitset<MAX_BITS>               positionToBinary(std::vector<uint>&);

std::vector<uint>                   binaryToPosition(std::bitset<MAX_BITS>&);

uint                                getParticleNumberBetween(std::bitset<MAX_BITS>& posP, uint from, uint to);

uint                                countDoublons(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2);

uint                                countSimilarities(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, uint nSites, uint& lastIndexWrong1, uint& lastIndexWrong2);

std::bitset<MAX_BITS>               expandBitSetUnderConstraint(std::bitset<MAX_BITS> reducedBS, std::bitset<MAX_BITS> constraintBS, uint sizeReducedBS = MAX_BITS);


bool                                canHopHalf(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, unsigned int n_pH, unsigned int& idFrom, unsigned int& idTo);


arma::cx_dvec                       expmv(std::complex<double> t, arma::sp_cx_dmat& A, arma::cx_dvec& v, double tol = 1.e-7, unsigned int m = 30);



class FHEDSolver  	
{
public:
    FHEDSolver(FHSystem* fhSystem);

    void setResultFile(std::string workName, std::string folderNameName);
    void setTimeParameters(double t_max, double dt, unsigned int nPerTS);
    void setComputeMomenta(bool yn);

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
    std::vector<std::array<std::bitset<MAX_BITS>,2>> basisArray;              
    std::vector<std::bitset<MAX_BITS>> c_pu;
    std::vector<std::bitset<MAX_BITS>> c_pd;
// Initial conditions :
    std::bitset<MAX_BITS> posPartU;
    std::bitset<MAX_BITS> posPartD;
    arma::cx_dvec psi_0;                                            //<- Initial wave function

// Current Wave function :
    arma::cx_dvec psi;                                              //<- Current wave function
// operators :
    arma::sp_cx_dmat H_U;
    arma::sp_cx_dmat H_t;
    arma::sp_cx_dmat H_J;
    arma::sp_cx_dmat H;

    //arma::cx_dmat Hdense;
    //std::vector<arma::sp_cx_dmat> NiiuSet;
    //std::vector<arma::sp_cx_dmat> NiidSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijuSet;
    std::vector<std::vector<arma::sp_cx_dmat>> NijdSet;

    //Momentum space :
    bool computeMomenta;
    uint nOfMomentumStates;
    arma::dmat momentaVec;
    arma::dvec momentaNRJ;

    unsigned int getBasisIndex(std::bitset<MAX_BITS>& posPU, std::bitset<MAX_BITS>& posPD );          //<- Get the basis index of the basis vector corresponding to particle positions 'posPU' and 'posPD'


    std::complex<double> U_os_psil(unsigned int phi_l);                                     //<- Compute the interation energy of the basis vector 'phi_l'


    void buildHamiltonian_U();                                                              //<- Build and store the operator "Interaction part of the Hamiltonian"
    void buildHamiltonian_tM();                                                              //<- Build and store the operator "Tunneling part of the Hamiltonian"
    void buildHamiltonian_J();                                                              //<- Build and store the operator "Exchange coupling part of the Hamiltonian"
    void buildPairOperators();                                                              //<- Build and store all the usefull number and tunneling operators.
    arma::sp_cx_dmat numberOperatorSpinUp(unsigned int iSite);                              //<- Build the operator "Number of up-spin particles on site 'iSite'".
    arma::sp_cx_dmat numberOperatorSpinDown(unsigned int iSite);                            //<- Build the operator "Number of down-spin particles on site 'iSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinUp(unsigned int iSite, unsigned int jSite);       //<- Build the operator "Tunneling of up-spin particles between sites 'iSite' and 'jSite'".
    arma::sp_cx_dmat tunnelingOperatorSpinDown(unsigned int iSite, unsigned int jSite);     //<- Build the operator "Tunneling of down-spin particles between sites 'iSite' and 'jSite'".

};

