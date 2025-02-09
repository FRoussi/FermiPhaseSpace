#pragma once

#pragma once

#include <vector>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>



#include "json/json.h"

using std::vector;

vector<vector<double>> makeBCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a);
vector<vector<double>> makeBCCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a);
vector<vector<double>> makeFCCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a);
vector<vector<double>> makeHoneyCombs(unsigned int nCombs, double a);

void getNextIndex(unsigned int& i, unsigned int& j);

arma::dmat armaMatfromVectorsVector(vector<vector<double>>& vecVec);

template <typename T>       bool isIn(T w, vector<T>& v);

template <typename T>       bool isIn(T w, vector<T>& v, uint& firstIndex);

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

class FHLattice  	
{



public:
    FHLattice(bool talk = true);
    FHLattice(int x, bool talk = true);
    FHLattice(int x, int y, int z, bool talk = true);
  //  FHSystem(std::string fileName, std::string folderName = "");
    FHLattice(Json::Value rootSystem, std::string folderName, bool talk = true);
    FHLattice(std::string latticeTypeString, vector<double> parameters);
    FHLattice(FHLatticeType latType, vector<double> parameters);

    FHLattice(FHLattice* fhLattice);

    // GETTERS
    // site positions
    const unsigned int                          getnSites();
    const unsigned int                          getnDimension();
    const arma::drowvec                         getSitePosition(uint iSite);
    const arma::dmat*                           getSitePositions();
    const double                                getCoordinate(uint dim, uint pos);
    const arma::drowvec                         getLatticeSitePosition(uint iSite);
    const double                                getLatticeCoordinate(uint dim, uint pos);
    const FHLatticeType                         getLatticeType();
    std::string                                 getLatticeName();
    // site names
    bool                                        sitesHaveNames();
    uint                                        getSiteNameNumber();
    std::vector<std::string>                    getSiteNames();
    std::string                                 getSiteName(uint iSite);
    arma::uvec                                  getSiteTypes();
    uint                                        getSiteType(uint iSite);
    std::vector<std::string>                    getSiteTypeNames();
    // system symmetry      
    bool                                        isPeriodic();
    arma::dmat                                  getLatticeMatrix();
    arma::drowvec                               getLatticeZeroVector();


    double                                      getNeighbourDistance();
    double                                      getSamePointDistance();
    double                                      computeDistance(uint iSite, uint jSite);
    double                                      getDistance(uint iSite, uint jSite);
    bool                                        areNeighbour(uint iSite, uint jSite);
    bool                                        areSameSite(uint iSite, uint jSite);

    const arma::dmat*                           getDistanceMatrix();
    const std::vector<std::vector<bool>>*       getNeighbourMatrix();



    // SETTERS      

    void                                        setPositions(const arma::dmat& positions, std::string sitesName = "C");
    void                                        setPositions(const arma::dmat& positions, std::vector<std::string> siteNames);

    bool                                        setSiteTypes(arma::uvec sTypes);
    void                                        setSiteTypeNames(vector<std::string> typeNames);  

    void                                        setSymmetry(arma::dmat latticeMat); 
    void                                        setLatticeZero(arma::drowvec latticeZ);
    void                                        setNeighbourDistance(double neigbD); 
    void                                        setSamePointDistance(double samePointD); 

    void                                        setLatticeName(std::string newName);

    void                                        computeMatrices();

    // XYZ file     
    void                                        readPositions(std::string xyzFile);
    void                                        writePositions(std::string xyzFile, std::string folderName = "");

    // Json     
    Json::Value                                 writeFHLatticeJson();

    ~FHLattice();


private:

    // SITES
    arma::dmat                      sitePositions;
    arma::uvec                      siteTypes;
    vector<std::string>             siteTypeNames;
    vector<std::string>             siteNames;

// LATTICE
    bool                            isPeriodicBool;
    arma::dmat                      latticeMatrix;
    arma::dmat                      inverseLatticeMatrix;
    arma::drowvec                   latticeZeroVector;

    // TYPE
    FHLatticeType                   latticeType;
    vector<int>                     intParameters;
    vector<double>                  floatParameters; 

// Connections
    double                          neighbourDistance = 1;
    double                          samePointDistance = 0.001;
    arma::dmat                      distanceMatrix;
    vector<vector<bool>>            neighbourMatrix;
    vector<vector<bool>>            sameSiteMatrix;

// Lattice type names
    std::map<FHLatticeType, std::string> latticeTypeNames;
    std::map<std::string, FHLatticeType> latticeTypeNamesReverse;
    void init_latticeNameMap();



    std::string                     latticeName;
    bool                            talkative;


};