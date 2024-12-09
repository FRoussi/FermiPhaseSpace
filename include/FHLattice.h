#pragma once

#pragma once

#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

#include "json/json.h"



arma::dmat armaMatfromVectorsVector(std::vector<std::vector<double>>& vecVec);

template <typename T>       bool isIn(T w, std::vector<T>& v);

template <typename T>       bool isIn(T w, std::vector<T>& v, uint& firstIndex);

class FHLattice  	
{

    enum FHLatticeType {
        none,
        cubic,
        bulkCenteredCubic,
        faceCenteredCubic,
        triangular,
        graphene,
        moire,
        kagome
    };

public:
    FHLattice(bool talk = true);
    FHLattice(int x, bool talk = true);
    FHLattice(int x, int y, int z, bool talk = true);
  //  FHSystem(std::string fileName, std::string folderName = "");
    FHLattice(Json::Value rootSystem, std::string folderName, bool talk = true);

    FHLattice(FHLattice* fhLattice);

    // GETTERS
    // site positions
    const unsigned int                        getnSites();
    const unsigned int                        getnDimension();
    const arma::drowvec                       getSitePosition(uint iSite);
    const arma::dmat*                         getSitePositions();
    const double                              getCoordinate(uint dim, uint pos);
    const arma::drowvec                       getLatticeSitePosition(uint iSite);
    const double                              getLatticeCoordinate(uint dim, uint pos);
    const FHLatticeType                       getLatticeType();
    // site names
    bool                                sitesHaveNames();
    uint                                getSiteNameNumber();
    std::string                         getSiteName(uint iSite);
    arma::uvec                          getSiteTypes();
    uint                                getSiteType(uint iSite);
    // system symmetry
    bool                                isPeriodic();
    arma::dmat                          getLatticeMatrix();
    arma::drowvec                       getLatticeZeroVector();


    double                              computeDistance(uint iSite, uint jSite);
    double                              getDistance(uint iSite, uint jSite);
    bool                                areNeighbour(uint iSite, uint jSite);
    bool                                areSameSite(uint iSite, uint jSite);

    std::string                         getLatticeName();

    // SETTERS      

    void                                setPositions(const arma::dmat& positions, std::string sitesName = "C");
    void                                setPositions(const arma::dmat& positions, std::vector<std::string> siteNames);

    bool                                setSiteTypes(arma::uvec sTypes);
    void                                setSiteTypeNames(std::vector<std::string> typeNames);  

    void                                setSymmetry(arma::dmat latticeMat); 
    void                                setLatticeZero(arma::drowvec latticeZ);
    void                                setNeighbourDistance(double neigbD); 
    void                                setSamePointDistance(double samePointD); 

    void                                setLatticeName(std::string newName);

    void                computeMatrices();

    // XYZ file
    void                readPositions(std::string xyzFile);
    void                writePositions(std::string xyzFile, std::string folderName = "");

    // Json
    Json::Value             writeFHLatticeJson();

    ~FHLattice();


private:

    // SITES
    arma::dmat                      sitePositions;
    arma::uvec                      siteTypes;
    std::vector<std::string>        siteTypeNames;
    std::vector<std::string>        siteNames;

// LATTICE
    bool                            isPeriodicBool;
    arma::dmat                      latticeMatrix;
    arma::dmat                      inverseLatticeMatrix;
    arma::drowvec                   latticeZeroVector;

    // TYPE
    FHLatticeType                   latticeType;
    std::vector<int>                intParameters;
    std::vector<double>             floatParameters; 

// Connections
    double                          neighbourDistance = 1;
    double                          samePointDistance = 0.001;
    arma::dmat                      distanceMatrix;
    std::vector<std::vector<bool>>  neighbourMatrix;
    std::vector<std::vector<bool>>  sameSiteMatrix;


    std::string latticeName;
    bool talkative;


};