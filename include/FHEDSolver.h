#pragma once

#include <bitset>
//#include <climits>
#include <iostream>
#include <string>
#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
//#include <cmath>

#include "FHSolver.h"




#define MAX_BITS 32
#define M_PI 3.14159265358979323846


std::vector<std::vector<bool>>      getCombinations(unsigned int n, unsigned int p, unsigned int max_n);
unsigned int                        nChoosep(unsigned int n, unsigned int p);

bool                        operator<(const std::bitset<MAX_BITS>& x, const std::bitset<MAX_BITS>& y);

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


arma::cx_dvec                       expmv(std::complex<double> t, arma::sp_dmat& A, arma::cx_dvec& v, double tol = 1.e-7, unsigned int m = 30);



class FHEDSolver : public FHSolver
{
public:


    enum storeOperators{
        SOnone,
        SOneeded,
        SOall
    };

    FHEDSolver(FHSystem* fhSystem);


    void initiateSolver();

    double currentSimulationTime;

    void buildHamiltonian();

    void updateHamiltonian(double time);
    //void numberOperator();

    void step(double d_t);
    void runSimulation();


    arma::cx_dvec getTBGroundState();
    arma::cx_dvec getMFGroundState();

    ~FHEDSolver();

private:


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
    arma::sp_dmat H_U;
    arma::sp_dmat H_t;
   // arma::sp_dmat H_J;
    arma::sp_dmat H_externalPotential;
    arma::sp_dmat H;

    // Current one-particle density matrix :
    arma::cx_dmat Nut;
    arma::cx_dmat Ndt;
    bool oPDensityMatUpToDate;


    //arma::cx_dmat Hdense;
    //std::vector<arma::sp_cx_dmat> NiiuSet;
    //std::vector<arma::sp_cx_dmat> NiidSet;
    //std::vector<arma::dvec> NiiuSet;
    //std::vector<arma::dvec> NiidSet;
    std::vector<std::vector<arma::sp_dmat>> NijuSet;
    std::vector<std::vector<arma::sp_dmat>> NijdSet;

    //Momentum space :
    uint nOfMomentumStates;
    arma::dmat momentaVec;
    arma::dvec momentaNRJ;

    storeOperators storeNiju;
    storeOperators storeNijd;

    uint getBasisIndex(std::bitset<MAX_BITS>& posPU, std::bitset<MAX_BITS>& posPD );          //<- Get the basis index of the basis vector corresponding to particle positions 'posPU' and 'posPD'


    std::complex<double> U_os_psil(uint phi_l);                                     //<- Compute the interation energy of the basis vector 'phi_l'



    void buildHamiltonian_U();                                                              //<- Build and store the operator "Interaction part of the Hamiltonian"
    void buildHamiltonian_tM();                                                              //<- Build and store the operator "Tunneling part of the Hamiltonian"
    //void buildHamiltonian_J();                                                              //<- Build and store the operator "Exchange coupling part of the Hamiltonian"
    void buildHamiltonian_ExternalPotential(double t);                                      //<- Build and store the operator "External Potential part of the Hamiltonian"
    void buildPairOperators();                                                              //<- Build and store all the usefull number and tunneling operators.

    arma::sp_dmat buildOperatorFromDiagonal(uint bSize, arma::dvec diag);

    arma::dvec numberVectorOperatorSpinUp(uint iSite);                              //<- Build the operator "Number of up-spin particles on site 'iSite'".
    arma::dvec numberVectorOperatorSpinDown(uint iSite);                            //<- Build the operator "Number of down-spin particles on site 'iSite'".
    arma::sp_dmat numberOperatorSpinUp(uint iSite);                              //<- Build the operator "Number of up-spin particles on site 'iSite'".
    arma::sp_dmat numberOperatorSpinDown(uint iSite);                            //<- Build the operator "Number of down-spin particles on site 'iSite'".
    arma::sp_dmat tunnelingOperatorSpinUp(uint iSite, uint jSite);       //<- Build the operator "Tunneling of up-spin particles between sites 'iSite' and 'jSite'".
    arma::sp_dmat tunnelingOperatorSpinDown(uint iSite, uint jSite);     //<- Build the operator "Tunneling of down-spin particles between sites 'iSite' and 'jSite'".


    void        computeOneParticleDensityMatrix();

};

