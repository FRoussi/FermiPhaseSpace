#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <armadillo>

#include "json/json.h"


using std::vector;
using std::array;

/// @brief The structural type of a lattice.
enum class FHLatticeType {
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
arma::dmat armaMatfromVectorsPosition(const vector<array<double,3>>& vecPos);


/// @brief A lattice of sites: positions, names, symmetry and connectivity.
///
/// Sites are stored as rows of an (nSites x nDimensions) matrix.
/// "Lattice" coordinates of a site are its position expressed in the basis of
/// the lattice symmetry vectors (see getLatticeSitePosition).
class FHLattice
{
public:
    FHLattice(bool talk = true);
    FHLattice(int x, bool talk = true);
    FHLattice(int x, int y, int z, bool talk = true);
    FHLattice(Json::Value rootSystem, std::string folderName, bool talk = true);
    FHLattice(std::string latticeTypeString, std::vector<double> parameters);
    FHLattice(FHLatticeType latType, std::vector<double> parameters);
    FHLattice(FHLattice* fhLattice);

    ~FHLattice() = default;

    // GETTERS
    // site positions
    uint                                getnSites() const;
    uint                                getnDimension() const;
    arma::drowvec                       getSitePosition(uint iSite) const;
    const arma::dmat*                   getSitePositions() const;
    double                              getCoordinate(uint dim, uint iSite) const;
    arma::drowvec                       getLatticeSitePosition(uint iSite) const;
    double                              getLatticeCoordinate(uint dim, uint iSite) const;
    FHLatticeType                       getLatticeType() const;
    std::string                         getLatticeName() const;
    // site names
    bool                                sitesHaveNames() const;
    uint                                getSiteNameNumber() const;
    std::vector<std::string>            getSiteNames() const;
    std::string                         getSiteName(uint iSite) const;
    arma::uvec                          getSiteTypes() const;
    uint                                getSiteType(uint iSite) const;
    std::vector<std::string>            getSiteTypeNames() const;
    // system symmetry
    bool                                isPeriodic() const;
    arma::dmat                          getLatticeMatrix() const;
    arma::drowvec                       getLatticeZeroVector() const;

    // distances and connectivity
    double                              getNeighbourDistance() const;
    double                              getSamePointDistance() const;
    double                              computeDistance(uint iSite, uint jSite) const;
    double                              getDistance(uint iSite, uint jSite) const;
    bool                                areNeighbour(uint iSite, uint jSite) const;
    bool                                areSameSite(uint iSite, uint jSite) const;
    double                              evaluateNearestNeighbor() const;
    const arma::dmat*                   getDistanceMatrix() const;
    const std::vector<std::vector<bool>>* getNeighbourMatrix() const;

    // SETTERS
    void                                setPositions(const arma::dmat& positions, std::string siteName = "C");
    void                                setPositions(const arma::dmat& positions, std::vector<std::string> siteNames);
    bool                                setSiteTypes(arma::uvec siteTypes);
    void                                setSiteTypeNames(std::vector<std::string> typeNames);
    void                                setSymmetry(arma::dmat latticeMat);
    void                                setLatticeZero(arma::drowvec latticeZero);
    void                                computeInverseLatticeMatrix();
    void                                setNeighbourDistance(double neighbourD);
    void                                setSamePointDistance(double samePointD);
    void                                setLatticeName(std::string newName);
    void                                computeMatrices();

    // XYZ file
    void                                readPositions(std::string xyzFile);
    void                                writePositions(std::string xyzFile, std::string folderName = "");

    // Json
    Json::Value                         writeFHLatticeJson() const;

private:
    // Sites
    arma::dmat                          sitePositions;
    arma::uvec                          siteTypes;
    std::vector<std::string>            siteTypeNames;
    std::vector<std::string>            siteNames;

    // Lattice symmetry
    bool                                isPeriodicBool = false;
    arma::dmat                          latticeMatrix  = arma::dmat(3, 3, arma::fill::zeros);
    arma::dmat                          inverseLatticeMatrix = arma::dmat(3, 3, arma::fill::zeros);
    arma::drowvec                       latticeZeroVector = arma::drowvec(3, arma::fill::zeros);

    // Type and parameters
    FHLatticeType                       latticeType = FHLatticeType::none;
    std::vector<int>                    intParameters;
    std::vector<double>                 floatParameters;

    // Connectivity
    double                              neighbourDistance = 1;
    double                              samePointDistance = 0.0001;
    arma::dmat                          distanceMatrix;
    std::vector<std::vector<bool>>      neighbourMatrix;
    std::vector<std::vector<bool>>      sameSiteMatrix;

    // Lattice type names
    std::map<FHLatticeType, std::string> latticeTypeNames;
    std::map<std::string, FHLatticeType> latticeTypeNamesReverse;

    std::string                         latticeName;
    bool                                talkative = true;

    void                                init_latticeNameMap();
};
