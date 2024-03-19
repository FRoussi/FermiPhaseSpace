#pragma once

#include <vector>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

#include "fnctarma.h"
#include "utilitaries.h"



std::complex<double> U_os_nij(double U0, double alpha, double d);

class FHSystem
{
public:

    enum FHInitialConditionType {
        givenDiagonal,
        givenDensityMatrix,
        givenWaveFunctionBasisVector,
        givenWaveFunction,
        quenched
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
        tJNearestNeighbourConstant,
        tJGivenMatrix,
        tJBySiteType
    };

    enum FHTunnelingType {
        tHnone,
        nearestNeighbourConstant,
        nearestNeighbourBySiteType,
        tunnelingDistanceFunction,
        tunnelingGivenMatrix
    };

// CONSTRUCTORS

    FHSystem();
    FHSystem(int x);
    FHSystem(int x, int y, int z);
    FHSystem(std::string fileName, std::string folderName = "");


// GETTERS

    unsigned int            getnSites();
    FHInitialConditionType  getInitialCondition();
    unsigned int            getUpSpinParticleNumber();
    unsigned int            getDownSpinParticleNumber();
    bool                    sitesHaveNames();
    arma::drowvec           getSitePosition(int site);
    std::string             getSiteName(int site);
    double                  getCoordinate(int dim, int pos);


    FHOppositeSpinInteractionType      getOppositeSpinInteractionType();
    FHTunnelingType                    getTunnelingType();
    FHSameSpinInteractionType          getSameSpinInteractionType();
    tJExchangeCouplingType             getExchangeCouplingEnergyType();
    double                  getDistance(int i, int j);
    bool                    areNeighbour(int i, int j);
    std::complex<double>    getInteractionEnergy(int i, int j);
    std::complex<double>    getSameSpinInteractionEnergy(int i, int j);
    std::complex<double>    getTunnelingEnergy(int i, int j);
    std::complex<double>    getExchangeCouplingEnergy(int i, int j);
    std::complex<double>    getExternalPotential(int i, double t);
    arma::cx_dmat*          getUMatrixPointer();
    arma::cx_dmat*          getVMatrixPointer();
    arma::cx_dmat*          getJMatrixPointer();
    arma::cx_dmat*          gettMatrixPointer();
    arma::cx_dvec*          getExtPotPointer();

    std::string             getName();
    std::string             getFolder();

    std::vector<bool>       getInitialCondition_upSpin();
    std::vector<bool>       getInitialCondition_downSpin();

// SETTERS

    void                    setInteraction_OnSite(double U);
    void                    setTunneling_BetweenNeighbours(double tH);
    void                    setExchangeCoupling(double J);
    void                    setOppositeSpinInteractionType(FHOppositeSpinInteractionType UType);
    void                    setTunnelingType(FHTunnelingType tType);
    void                    setSameSpinInteractionType(FHSameSpinInteractionType VType);
    void                    setExchangeCouplingEnergyType(tJExchangeCouplingType JType);
    void                    setName(std::string name);
    bool                    setSiteTypes(arma::uvec sTypes);
    void                    setSiteTypeNames(std::vector<std::string> typeNames);         
    void                    setPositions(const arma::dmat& positions);
    void                    setInitialConditions(const std::vector<std::vector<unsigned int>>& initConditions);
    void                    setInitialConditions(FHInitialConditionType initialCondType, std::vector<unsigned int> partNumber );


    void                    writeFHSystem(std::string folderName = "");

    void initiateEmptyEnergyMatrices();


    ~FHSystem();

private:
    
    arma::dmat                      sitePositions;
    arma::uvec                      siteTypes;
    std::vector<std::string>        siteTypeNames;

    unsigned int                    upSpinParticleNumber;
    unsigned int                    downSpinParticleNumber;

    arma::dmat                      distanceMatrix;
    std::vector<std::vector<bool>>  neighbourMatrix;
    std::vector<std::vector<bool>>  sameSiteMatrix;

// ENERGY
    FHOppositeSpinInteractionType   UInteractionType;
    FHSameSpinInteractionType       VInteractionType;
    tJExchangeCouplingType          tJExCType;
    FHTunnelingType                 tunnelingType;

    arma::cx_dmat                   UMatrix;
    arma::cx_dmat                   VMatrix;
    arma::cx_dmat                   JMatrix;
    arma::cx_dmat                   tMatrix;
    arma::cx_dvec                   externalPotential;

// INITIAL CONDITIONS
    FHInitialConditionType          initialCondition;
    std::vector<bool>               initConditions_upSpin;
    std::vector<bool>               initConditions_downSpin;

    std::string                     systemName;
    std::string                     systemFolder;


    void                            computeMatrices();


//    bool isIn(std::string w, std::vector<std::string> v);

};