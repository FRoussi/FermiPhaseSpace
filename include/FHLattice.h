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
using std::array;

vector<array<double,3>> makeBCGrid(uint nx, uint ny, uint nz, double a);
arma::dmat findBCGridBox(uint nx, uint ny, uint nz, double a);
vector<array<double,3>> makeBCCGrid(uint nx, uint ny, uint nz, double a);
vector<array<double,3>> makeFCCGrid(uint nx, uint ny, uint nz, double a);
vector<array<double,3>> makeHoneyCombs(uint nCombs, double a);
arma::dmat findHoneyCombGridBox(uint nCombs, double a);
vector<array<double,3>> makeGraphene(int nx, int ny, double a);
arma::dmat findGrapheneGridBox(int nx, int ny, double a);
vector<array<double,3>> makeMoireHoneyCombs(uint nCombs, double a, double b, double theta);
void getNextIndexMoireHoneyCombs(int &n, int &m, int &l);
vector<array<double,3>> makeMoireGraphene(int nx, int ny, double a, double b, double theta);
vector<array<double,3>> makeTetragones(int nx, int ny, int nz, double a);
arma::dmat findTetragonesGridBox(int nx, int ny, int nz, double a);

void getNextIndex(uint& i, uint& j);

arma::dmat armaMatfromVectorsPosition(vector<array<double,3>>& vecPos);

template <typename T>       bool isIn(T w, vector<T>& v);

template <typename T>       bool isIn(T w, vector<T>& v, uint& firstIndex);

    enum FHLatticeType {
        none,
        cubic,
        bulkCenteredCubic,
        faceCenteredCubic,
        triangular,
        honeyCombs,
        graphene,
        moireHC,
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
    const uint                                  getnSites();
    const uint                                  getnDimension();
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
    void                                        computeInverseLatticeMatrix();
    void                                        setNeighbourDistance(double neigbD); 
    void                                        setSamePointDistance(double samePointD); 

    void                                        setLatticeName(std::string newName);

    void                                        computeMatrices();
    double                                      evaluateNearestNeighbor();

    // XYZ file     
    void                                        readPositions(std::string xyzFile);
    void                                        writePositions(std::string xyzFile, std::string folderName = "");

    // Json     
    Json::Value                                 writeFHLatticeJson();

    ~FHLattice();


private:

    // SITES
    arma::dmat                              sitePositions;
    arma::uvec                              siteTypes;
    vector<std::string>                     siteTypeNames;
    vector<std::string>                     siteNames;

// LATTICE      
    bool                                    isPeriodicBool;
    arma::dmat                              latticeMatrix;
    arma::dmat                              inverseLatticeMatrix;
    arma::drowvec                           latticeZeroVector;

// TYPE     
    FHLatticeType                           latticeType;
    vector<int>                             intParameters;
    vector<double>                          floatParameters; 

// Connections      
    double                                  neighbourDistance = 1;
    double                                  samePointDistance = 0.0001;
    arma::dmat                              distanceMatrix;
    vector<vector<bool>>                    neighbourMatrix;
    vector<vector<bool>>                    sameSiteMatrix;

// Lattice type names
//MapDoubleEntry<int>           latticeNamesMap;
    //MapDoubleEntry<FHLatticeType>           latticeNamesMap;
    std::map<FHLatticeType, std::string>    latticeTypeNames;
    std::map<std::string, FHLatticeType>    latticeTypeNamesReverse;

// Names
    std::string                             latticeName;
    bool                                    talkative;



    void                                    init_latticeNameMap();

};