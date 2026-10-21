#pragma once

#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

#include "FHLattice.h"
#include "FHCoreHole.h"
#include "FHTunnelingInteraction.h"
#include "FHCoulombInteraction.h"
#include "FHExternalPotential.h"
#include "FHExternalPotentialPulse.h"
#include "json/json.h"




#define FPS_VERSION_MAJOR 0
#define FPS_VERSION_MINOR 5
#define FPS_VERSION_PATCH 0


std::complex<double> U_os_nij(double U0, double alpha, double d);

enum FHSystemType {
    SysVanilla,
    SysCoreHole
};

class FHSystem
{
public:

    enum FHInitialConditionType {
        spinWave,
        chargeWave,
        givenDiagonal,
        givenDensityMatrix,
        givenWaveFunctionBasisVector,
        givenWaveFunction,
        tbGroundState,
        mfGroundState,
        groundState
    };


private:
    
    FHSystemType                    systemType{FHSystemType::SysVanilla};
    FHLattice                       lattice;

// coreHole
    FHCoreHole                      coreHole;

// sites
    uint                            nSites{0};

// PARTICLES
    unsigned int                    upSpinParticleNumber{0};
    unsigned int                    downSpinParticleNumber{0};

// ENERGY
    FHCoulombInteraction            coulombInteraction;
    FHTunnelingInteraction          tunnelingInteraction;
    FHExternalPotential             externalPotential;

// INITIAL CONDITIONS
    FHInitialConditionType          initialCondition;
    std::vector<bool>               initConditions_upSpin;
    std::vector<bool>               initConditions_downSpin;

// Rotated representations
    // Momentum space :
    bool TBGroundStateComputed{false};
    arma::cx_dmat          tbgs_Nu;
    arma::cx_dmat          tbgs_Nd;
    arma::cx_dmat          tbgs_V;
    arma::dvec             tbgs_E;

    // MF ground state :
    bool MFGroundStateComputed{false};
    arma::cx_dmat           mfgs_Nu;
    arma::cx_dmat           mfgs_Nd;
    arma::cx_dmat           mfgs_Vu;
    arma::cx_dmat           mfgs_Vd;
    arma::dvec              mfgs_Eu;
    arma::dvec              mfgs_Ed;
    uint                    iterationsForMFState;
    double                  residueForMFState;

// SYSTEM
    std::string                     systemName;
    std::string                     systemFolder;
    Json::Value                     originalJsonSystem;

    bool talkative{false};



public:
// CONSTRUCTORS

    FHSystem(bool talk = true);
    FHSystem(int x, bool talk = true);
    FHSystem(int x, int y, int z, bool talk = true);
  //  FHSystem(std::string fileName, std::string folderName = "");
    FHSystem(Json::Value rootSystem, std::string folderName, bool talk = true);
    FHSystem(FHLattice lattice);


// GETTERS
    FHSystemType                        getSystemType() const;
    FHLattice*                          getLattice() ;
    FHCoreHole*                         getCoreHole() ;
    // site positions
    unsigned int                        getnSites() const;
    unsigned int                        getnDimension() const;
    arma::drowvec                       getSitePosition(uint iSite) const;
    double                              getCoordinate(uint dim, uint pos) const;
    arma::drowvec                       getLatticeSitePosition(uint iSite) const;
    double                              getLatticeCoordinate(uint dim, uint pos) const;
    // site names
    bool                                sitesHaveNames() const;
    std::string                         getSiteName(uint site) const;
    // system symmetry
    bool                                isPeriodic() const;
    arma::dmat                          getLatticeMatrix() const;
    arma::drowvec                       getLatticeZeroVector() const;

    // 
    bool                                isSystemTimeDependant() const;
    bool                                isExternalPotentialNull() const;
    bool                                isTBGroundStateComputed() const;
    bool                                isMFGroundStateComputed() const;

   // double                              computeDistance(uint iSite, uint jSite);
    double                              getDistance(uint iSite, uint jSite) const;
    bool                                areNeighbour(uint iSite, uint jSite) const;

    FHInitialConditionType              getInitialCondition() const;
    unsigned int                        getUpSpinParticleNumber() const;
    unsigned int                        getDownSpinParticleNumber() const;

    FHCoulombInteractionType            getCoulombInteractionType() const;
    FHTunnelingInteractionType          getTunnelingInteractionType() const;
    FHExternalPotentialType             getExternalPotentialType() const;
    //FHCoulombInteractionType            getSameSpinInteractionType();
   // tJExchangeCouplingType              getExchangeCouplingEnergyType();

    std::complex<double>                getCoulombInteraction(int iSite, int jSite) ;
    std::complex<double>                getSameSpinInteractionEnergy(int iSite, int jSite) ;
    std::complex<double>                getTunnelingInteraction(int iSite, int jSite) ;
    //std::complex<double>                getExchangeCouplingEnergy(int iSite, int jSite);
    std::complex<double>                getExternalPotential(int iSite, double time) ;
    std::complex<double>                getExternalPotential(int iSite, int jSite, double time) ;
    arma::cx_dmat                       getExternalPotential(double time);
    void                                updateExternalPotential(double time);

    const arma::cx_dmat&                getUMatrix() const;
    const arma::cx_dmat&                getVMatrix() const;
    //const arma::cx_dmat*                getJMatrixPointer();
    const arma::cx_dmat&                gettMatrix() const;
    const arma::cx_dmat&                getExtPot() const;
    const arma::cx_dmat*                get_tbgs_Nu_Pointer();
    const arma::cx_dmat*                get_tbgs_Nd_Pointer();
    const arma::cx_dmat*                get_tbgs_V_Pointer();
    const arma::dvec*                   get_tbgs_E_Pointer();
    const arma::cx_dmat*                get_mfgs_Nu_Pointer();
    const arma::cx_dmat*                get_mfgs_Nd_Pointer();
    const arma::cx_dmat*                get_mfgs_Vu_Pointer();
    const arma::cx_dmat*                get_mfgs_Vd_Pointer();
    const arma::dvec*                   get_mfgs_Eu_Pointer();
    const arma::dvec*                   get_mfgs_Ed_Pointer();


    std::vector<bool>                   getInitialCondition_upSpin();
    std::vector<bool>                   getInitialCondition_downSpin();

    std::string                         getName() const;
    std::string                         getFolder() const;
    std::string                         getFPSVersion() const;



// SETTERS      
    void                                setLattice(FHLattice newLat);

    void                                setPositions(const arma::dmat& positions, std::string sitesName = "C");
    void                                setPositions(const arma::dmat& positions, std::vector<std::string> siteNames);

    bool                                setSiteTypes(arma::uvec sTypes);
    void                                setSiteTypeNames(std::vector<std::string> typeNames);  

    void                                setSymmetry(arma::dmat latticeMat); 
    void                                setLatticeZero(arma::drowvec latticeZ);
    void                                setNeighbourDistance(double neigbD); 
    void                                setSamePointDistance(double samePointD); 


    void                                setInitialConditionsDiagonal(const std::vector<std::vector<unsigned int>>& initConditions);
    void                                setInitialConditions(FHInitialConditionType initialCondType);
    void                                setInitialConditions(FHInitialConditionType initialCondType, std::vector<unsigned int> partNumber );

    void                                setBasisInteractions();

   // void                                setOppositeSpinInteractionType(FHOppositeSpinInteractionType UType);
   // void                                setTunnelingType(FHTunnelingType tType);
   // void                                setSameSpinInteractionType(FHSameSpinInteractionType VType);
    //void                                setExchangeCouplingEnergyType(tJExchangeCouplingType JType);

   // void                                setUInteraction_OnSite(double U0, bool add = false);
   // void                                setUInteraction_GivenDiagonal(arma::dvec& Ua, bool add = false);
    //void                                setUInteraction_GivenMatrix(arma::dmat& Uab, bool add = false);
    //void                                setUInteraction_BySiteType(arma::dvec& Us, bool add = false);
   // void                                setTunneling_BetweenNeighbours(double tH1, bool add = false);
    //void                                setTunneling_atDistance(double tH1, double dist, bool add = false);
    //void                                setTunneling_SamePosition(double tH1, bool add = false);
    //void                                setTunneling_BySiteType(arma::dmat& tHs1s2, double d, bool add = false);
   // void                                setTunneling_GivenMatrix(arma::dmat tH, bool add = false);
   // void                                setExchangeCoupling(double J, bool add = false);

    void                                setName(std::string name);

    Json::Value                         writeFHSystemJson();
    Json::Value                         getOriginalJsonSystem();
    Json::Value                         writeTBGroundStateJson();
    Json::Value                         writeMFGroundStateJson();
    //void                                writeFHSystemBis(std::string folderName = "");



// External Potential
    void                                makePotential_pulse();
   // void                                readPotential(std::string fileName, std::string folderName = "");


// WannierTools

    void writeWannierToolsDatFile(std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsInBandFile(std::string wtIBFile, std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsDosFile(std::string wtDosFile, std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsInBandUnfoldFile(std::string wtIBUFile, std::string wtDatFile, std::string folderName = "");
    void WannierToolsSystemInofStream(std::ofstream& wannierFileStream);


// INITIATION
    void                                initiateNullInteraction();
    void                                computeTBGroundState(bool withExternalPotential = false);
    void                                computeMFGroundState(bool withExternalPotential = false);


    ~FHSystem();




//    bool isIn(std::string w, std::vector<std::string> v);

};