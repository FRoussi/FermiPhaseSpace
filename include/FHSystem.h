#pragma once

#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

#include "FHLattice.h"
#include "FHExternalPotential.h"
#include "FHExternalPotentialPulse.h"
#include "json/json.h"






std::complex<double> U_os_nij(double U0, double alpha, double d);

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
        quenched,
        groundState
    };

    enum FHOppositeSpinInteractionType {
        Unone,
        UonSiteConstant,
        UonSiteBySiteType,
        UonSiteGivenDiagonal,
        UdistanceFunction,
        UGivenMatrix
    };

    enum FHSameSpinInteractionType {
        Vnone,
        VonSiteConstant,
        VonSiteBySiteType,
        VonSiteGivenDiagonal,
        VdistanceFunction,
        VGivenMatrix
    };

    enum tJExchangeCouplingType {
        tJJnone,
        tJNearestNeighboursConstant,
        tJGivenMatrix,
        tJBySiteType
    };

    enum FHTunnelingType {
        tHnone,
        nearestNeighboursConstant,
        nearestNeighboursBySiteType,
        tunnelingDistanceFunction,
        tunnelingGivenMatrix
    };

    enum ExternalPotentialType {
        ExtPotNone,
        ExtPotConstant,
        ExtPotPulse
    };

// CONSTRUCTORS

    FHSystem(bool talk = true);
    FHSystem(int x, bool talk = true);
    FHSystem(int x, int y, int z, bool talk = true);
  //  FHSystem(std::string fileName, std::string folderName = "");
    FHSystem(Json::Value rootSystem, std::string folderName, bool talk = true);


// GETTERS
FHLattice                               getLattice();
    // site positions
    unsigned int                        getnSites();
    unsigned int                        getnDimension();
    arma::drowvec                       getSitePosition(uint iSite);
    double                              getCoordinate(uint dim, uint pos);
    arma::drowvec                       getLatticeSitePosition(uint iSite);
    double                              getLatticeCoordinate(uint dim, uint pos);
    // site names
    bool                                sitesHaveNames();
    std::string                         getSiteName(uint site);
    // system symmetry
    bool                                isPeriodic();
    arma::dmat                          getLatticeMatrix();
    arma::drowvec                       getLatticeZeroVector();

    // 
    bool                                isSystemTimeDependant();
    bool                                isExternalPotentialNull();

    double                              computeDistance(uint iSite, uint jSite);
    double                              getDistance(uint iSite, uint jSite);
    bool                                areNeighbour(uint iSite, uint jSite);

    FHInitialConditionType              getInitialCondition();
    unsigned int                        getUpSpinParticleNumber();
    unsigned int                        getDownSpinParticleNumber();

    FHOppositeSpinInteractionType       getOppositeSpinInteractionType();
    FHTunnelingType                     getTunnelingType();
    FHSameSpinInteractionType           getSameSpinInteractionType();
    tJExchangeCouplingType              getExchangeCouplingEnergyType();

    std::complex<double>                getInteractionEnergy(int iSite, int jSite);
    std::complex<double>                getSameSpinInteractionEnergy(int iSite, int jSite);
    std::complex<double>                getTunnelingEnergy(int iSite, int jSite);
    std::complex<double>                getExchangeCouplingEnergy(int iSite, int jSite);
    std::complex<double>                getExternalPotential(int iSite, double time);
    std::complex<double>                getExternalPotential(int iSite, int jSite, double time);
    arma::cx_dmat                       getExternalPotential(double time);
    void                                updateExternalPotential(double time);

    const arma::cx_dmat*                getUMatrixPointer();
    const arma::cx_dmat*                getVMatrixPointer();
    const arma::cx_dmat*                getJMatrixPointer();
    const arma::cx_dmat*                gettMatrixPointer();
    const arma::cx_dmat*                getExtPotPointer();


    std::string                         getName();
    std::string                         getFolder();

    std::vector<bool>                   getInitialCondition_upSpin();
    std::vector<bool>                   getInitialCondition_downSpin();


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

    void                                setOppositeSpinInteractionType(FHOppositeSpinInteractionType UType);
    void                                setTunnelingType(FHTunnelingType tType);
    void                                setSameSpinInteractionType(FHSameSpinInteractionType VType);
    void                                setExchangeCouplingEnergyType(tJExchangeCouplingType JType);

    void                                setUInteraction_OnSite(double U0, bool add = false);
    void                                setUInteraction_GivenDiagonal(arma::dvec& Ua, bool add = false);
    void                                setUInteraction_GivenMatrix(arma::dmat& Uab, bool add = false);
    void                                setUInteraction_BySiteType(arma::dvec& Us, bool add = false);
    void                                setTunneling_BetweenNeighbours(double tH1, bool add = false);
    void                                setTunneling_atDistance(double tH1, double dist, bool add = false);
    void                                setTunneling_SamePosition(double tH1, bool add = false);
    void                                setTunneling_BySiteType(arma::dmat& tHs1s2, double d, bool add = false);
    void                                setTunneling_GivenMatrix(arma::dmat tH, bool add = false);
    void                                setExchangeCoupling(double J, bool add = false);

    void                                setName(std::string name);

    Json::Value                         writeFHSystemJson();
    Json::Value                         getOriginalJsonSystem();
    //void                                writeFHSystemBis(std::string folderName = "");



// External Potential
    void                                makePotential_pulse();
    void                                readPotential(std::string fileName, std::string folderName = "");


// WannierTools

    void writeWannierToolsDatFile(std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsInBandFile(std::string wtIBFile, std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsDosFile(std::string wtDosFile, std::string wtDatFile, std::string folderName = "");
    void writeWannierToolsInBandUnfoldFile(std::string wtIBUFile, std::string wtDatFile, std::string folderName = "");
    void WannierToolsSystemInofStream(std::ofstream& wannierFileStream);


// INITIATION
    void                                initiateEmptyEnergyMatrices();


    ~FHSystem();

private:
    

    FHLattice                       lattice;


// PARTICLES
    unsigned int                    upSpinParticleNumber;
    unsigned int                    downSpinParticleNumber;


// ENERGY
    FHOppositeSpinInteractionType   UInteractionType;
    FHSameSpinInteractionType       VInteractionType;
    tJExchangeCouplingType          tJExCType;
    FHTunnelingType                 tunnelingType;
    ExternalPotentialType           externalPotentialType;

    arma::cx_dmat                   UMatrix;
    arma::cx_dmat                   VMatrix;
    arma::cx_dmat                   JMatrix;
    arma::cx_dmat                   tMatrix;
    FHExternalPotential             externalPotential;

// INITIAL CONDITIONS
    FHInitialConditionType          initialCondition;
    std::vector<bool>               initConditions_upSpin;
    std::vector<bool>               initConditions_downSpin;

// SYSTEM
    std::string                     systemName;
    std::string                     systemFolder;
    Json::Value                     originalJsonSystem;



    bool talkative;


//    bool isIn(std::string w, std::vector<std::string> v);

};