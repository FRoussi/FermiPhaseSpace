#pragma once

#include <vector>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

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

    enum FHInteractionType {
        onSiteConstant,
        onSiteBySiteType,
        onSiteGivenDiagonal,
        distanceFunction,
        matrixGiven
    };

    enum tJExchangeCouplingType {
        tJNearestNeighbourConstant
    };

    enum FHTunnelingType {
        nearestNeighbourConstant,
        nearestNeighbourBySiteType,
        tunnelingDistanceFunction,
        tunnelingMatrixGiven
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

    double                  getDistance(int i, int j);
    bool                    areNeighbour(int i, int j);
    std::complex<double>    getInteractionEnergy(int i, int j);
    std::complex<double>    getTunnelingEnergy(int i, int j);
    std::complex<double>    getExchangeCouplingEnergy(int i, int j);
    arma::cx_dmat*          getUMatrixPointer();
    arma::cx_dmat*          getJMatrixPointer();
    arma::cx_dmat*          gettMatrixPointer();

    std::string             getName();
    std::string             getFolder();

    std::vector<bool>       getInitialCondition_upSpin();
    std::vector<bool>       getInitialCondition_downSpin();

// SETTERS

    void                    setInteraction_OnSite(double U);
    void                    setTunneling_BetweenNeighbours(double J);
    void                    setName(std::string name);
    bool                    setSiteTypes(arma::uvec sTypes);
    void                    setSiteTypeNames(std::vector<std::string> typeNames);                



    void                    writeFHSystem(std::string folderName = "");

    void        initiateSystem(const arma::dmat& positions, const std::vector<std::vector<unsigned int>>& initConditions );
    void        initiateSystem(const arma::dmat& positions, FHInitialConditionType initialCondType, std::vector<unsigned int> partNumber );

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

    arma::cx_dmat                   UMatrix;
    arma::cx_dmat                   JMatrix;
    arma::cx_dmat                   tMatrix;

    FHInitialConditionType          initialCondition;
    std::vector<bool>               initConditions_upSpin;
    std::vector<bool>               initConditions_downSpin;

    std::string                     systemName;
    std::string                     systemFolder;


    void        computeMatrices();



//    bool isIn(std::string w, std::vector<std::string> v);


};