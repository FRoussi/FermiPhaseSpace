#include "../include/FHLattice.h"


FHLattice::FHLattice(bool talk)
{

    talkative = talk;

    latticeType = FHLatticeType::cubic;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);

    arma::dmat pos = arma::dmat(2,  3, arma::fill::zeros);
    pos(1,0) = 1;


    setPositions(pos);
   

}


FHLattice::FHLattice(int x, bool talk)
{
    talkative = talk;

   latticeType = FHLatticeType::cubic;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);

    arma::dmat pos = arma::dmat(1,  x, arma::fill::zeros);
    for (int d=0; d<x; d++)
    {
        pos[d] = (double)d;
    }

    int nPart = std::floor(x/2);
    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>(nPart,0);
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>(nPart,0);
    for(int itPart = 0; itPart <nPart; itPart++)
    {
        initCond_up[itPart] = itPart * 2;
        initCond_down[itPart] = itPart * 2 + 1;
    }

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;


    setPositions(pos);
  

}


FHLattice::FHLattice(int x, int y, int z, bool talk)
{

    talkative = talk;

   latticeType = FHLatticeType::cubic;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);

    arma::dmat pos = arma::dmat(3,  x * y * z, arma::fill::zeros);
    int tabIndex = 0;
    for (int i_z = 0; i_z < z; i_z++)
    {
        for (int i_y = 0; i_y < y; i_y++)
        {
            for (int i_x = 0; i_x < x; i_x++)
            {   
                pos(0,tabIndex) = double(i_x);
                pos(1,tabIndex) = double(i_y);
                pos(2,tabIndex) = double(i_z); 
                tabIndex +=1;
            }
        }       
    }



    setPositions(pos);
     

}



FHLattice::FHLattice(Json::Value rootSystem, std::string folderName, bool talk)
{

    talkative = talk;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);

// Positions
    std::string positionFile = rootSystem.get("positions", "").asString();

    if (folderName !="")
    {
        positionFile = folderName + "/" + positionFile;
    }

    // std::cout << positionFile << std::endl;
    readPositions(positionFile);


// Periodicity
    isPeriodicBool = rootSystem.get("periodic", false).asBool();


    Json::Value latMat = rootSystem["latticeMatrix"];
    if (latMat != Json::Value::null)
    {
        for (uint i = 0; i<3; i++)
        {
            Json::Value latRow = latMat[i]; 
            for (uint j = 0; j<3; j++)
            {
                latticeMatrix(i,j) = latRow[j].asDouble();
            }
        }
        bool isInv = arma::inv(inverseLatticeMatrix, latticeMatrix);
        if (!isInv)
        {
            std::cout << "Caution : the lattice matrix is singular." << std::endl;
        }
    }
    

    Json::Value latRow = rootSystem["latticeZero"];
    if (latRow != Json::Value::null)
    {
        for (uint j = 0; j<3; j++)
        {
            latticeZeroVector(j) = latRow[j].asDouble();
        }
    }

    



}



/// @brief getter - number of sites
/// @return number of sites
const unsigned int                                FHLattice::getnSites()
{
    return sitePositions.n_rows;
}

/// @brief getter - number of dimension (we always work in dim 3 though)
/// @return number of dimensions (if not 3, weird)
const unsigned int                                FHLattice::getnDimension()
{
    return sitePositions.n_cols;
}

/// @brief getter - position row vector of site 'iSite'
/// @param iSite a site number 
/// @return row vector of coordinates
const arma::drowvec                               FHLattice::getSitePosition(uint iSite)
{
    return sitePositions.row(iSite);
}

const arma::dmat*                         FHLattice::getSitePositions()
{
    return &sitePositions;
}

/// @brief getter - coordinate of a site
/// @param dim coordinate dimension
/// @param iSite site number
/// @return coordinate of site 'iSite'in dimension 'dim'
const double                                      FHLattice::getCoordinate(uint dim, uint iSite)
{ 
    if ((dim < sitePositions.n_cols) && (iSite < sitePositions.n_rows))
    {
        return sitePositions(iSite, dim);
    } 
    else
    {
        return 0;
    }
}

/// @brief getter - relative position row vector of site 'iSite'
/// @param iSite a site number 
/// @return row vector of relative coordinates
const arma::drowvec                       FHLattice::getLatticeSitePosition(uint iSite)
{
    //std::cout << sitePositions.row(iSite) << std::endl;
    //std::cout << inverseLatticeMatrix << std::endl;
   // std::cout << latticeZeroVector << std::endl;
    arma::drowvec pos = sitePositions.row(iSite) - latticeZeroVector;
   // std::cout << inverseLatticeMatrix * pos.t() << std::endl;
    return (inverseLatticeMatrix * pos.t()).t();
}

const double                              FHLattice::getLatticeCoordinate(uint dim, uint iSite)
{
    if ((dim < sitePositions.n_cols) && (iSite < sitePositions.n_rows))
    {
        arma::drowvec v = getLatticeSitePosition(iSite);
        return v[dim];
    } 
    else
    {
        return 0;
    }
}

const FHLattice::FHLatticeType                       FHLattice::getLatticeType()
{
    return latticeType;
}

bool                                        FHLattice::sitesHaveNames()
{
        return (siteTypes.n_rows == sitePositions.n_rows);
}

uint                                FHLattice::getSiteNameNumber()
{
    return (uint)siteNames.size();
}

std::string                                 FHLattice::getSiteName(uint iSite)
{
    //return siteTypeNames[siteTypes[site]];
    return siteNames[iSite];
}

    arma::uvec                          FHLattice::getSiteTypes()
    {
        return siteTypes;
    }
    uint                                FHLattice::getSiteType(uint iSite)
    {
        return siteTypes[iSite];
    }

bool                                        FHLattice::isPeriodic()
{
    return isPeriodicBool;
}

arma::dmat                                  FHLattice::getLatticeMatrix()
{
    return latticeMatrix;
}

arma::drowvec                        FHLattice::getLatticeZeroVector()
{
    return latticeZeroVector;
}

double                                      FHLattice::computeDistance(uint iSite, uint jSite)
{
    arma::drowvec dVec = sitePositions.row(iSite) - sitePositions.row(jSite); 
    double d = arma::norm(dVec);
    //std::cout << isPeriodic() << std::endl;
    if (isPeriodic())
    {
        for (int zSym = -1; zSym < 2; zSym++)                             //<- symmetries in z
        {
            for (int ySym = -1; ySym < 2; ySym++)                         //<- symmetries in y
            {
                for (int xSym = -1; xSym < 2; xSym++)                     //<- symmetries in x
                {
                    arma::drowvec dVecSym = dVec + xSym * latticeMatrix.col(0).t() + ySym * latticeMatrix.col(1).t() + zSym * latticeMatrix.col(2).t(); //<- new vector with one of the symmetric lattice
                    d = std::min(d, arma::norm(dVecSym));                   //<- we take the minimum of them all
                    //std::cout << xSym << "  " << ySym << "  " << zSym << "  " << arma::norm(dVecSym) << "   " << d << std::endl;
                }
            }
        }
            
    }

    return d;
}

double                                      FHLattice::getDistance(uint i, uint j)
{
    return distanceMatrix(i,j);
}

bool                                        FHLattice::areNeighbour(uint i, uint j)
{
    return neighbourMatrix[i][j];
}

    bool                                FHLattice::areSameSite(uint i, uint j)
    {
        return sameSiteMatrix[i][j];
    }




// SETTERS

void                    FHLattice::setPositions(const arma::dmat& positions, std::string sitesName)
{
    sitePositions = positions;
    siteNames = std::vector<std::string>(sitePositions.n_rows, sitesName);
    computeMatrices();
}

void                    FHLattice::setPositions(const arma::dmat& positions, std::vector<std::string> sN)
{
    sitePositions = positions;
    siteNames = sN;
    computeMatrices();
}


bool                    FHLattice::setSiteTypes(arma::uvec sTypes)
{
    siteTypes = sTypes;
    return sitesHaveNames();
}

void                    FHLattice::setSiteTypeNames(std::vector<std::string> typeNames)
{
    siteTypeNames = typeNames;
}


void                    FHLattice::setSymmetry(arma::dmat symVectors)
{
    if ((symVectors.n_cols != 3)&& (symVectors.n_rows != 3))
    {
        std::cout << "We want a 3x3 lattice matrix" << std::endl;
    }
    else
    {
        isPeriodicBool = true;
        latticeMatrix = symVectors;
        inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
        bool isInv = arma::inv(inverseLatticeMatrix, latticeMatrix);
        if (!isInv)
        {
            std::cout << "Caution : the lattice matrix is singular." << std::endl;
        }
    }

    if (latticeZeroVector.n_rows !=3)
    {
        latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    }

}

void                                FHLattice::setLatticeZero(arma::drowvec latticeZ)
{
    latticeZeroVector = latticeZ;
}


void                                FHLattice::setNeighbourDistance(double neigbD)
{
    neighbourDistance = neigbD;
}


void                                FHLattice::setSamePointDistance(double samePointD)
{
    samePointDistance = samePointD;
}


void                    FHLattice::readPositions(std::string xyzFile)
{

    uint nPart = 0;
    std::vector<std::vector<double>> positions3;

    siteNames = std::vector<std::string>();
    std::string line;
    std::string coord;
     

    std::ifstream posFile(xyzFile);
    if (posFile.is_open())
    {
        // Particle number :
        std::getline(posFile, line);

        nPart = std::stoi(line);
        // Comment :
        std::getline(posFile, line);
        // Particles :
        while(std::getline(posFile, line))
        {
            std::stringstream ss(line);
            // Get the site name
            ss >> coord;
            siteNames.push_back(coord);
            // Get its coordinates :
            std::vector<double> coordinates;
            while (ss >> coord)//(std::getline (ss, coord, '\t')) 
            {
                double d = std::stod(coord);
                coordinates.push_back(d);
            }                      
            coordinates.resize(3);
            positions3.push_back(coordinates);
        }

        posFile.close();
    }
    else
    {
        std::cout << "Problem with the particle file " << xyzFile << std::endl;
    }


    if (nPart == positions3.size())
    {
        if (talkative)
        {
            std::cout << "   -- Position file read, " << nPart << " particles." << std::endl;
        }
        setPositions(armaMatfromVectorsVector(positions3));
    }
    else
    {
        std::cout << "Problem in the xyz file : " << nPart << " particles announced, " << positions3.size() << " found." << std::endl;
    }

}


/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void                 FHLattice::writePositions(std::string xyzFile, std::string folderName)
{
        // Writing of everything on the file :
    std::string fileAdress =  xyzFile + ".xyz";
    //std::cout << "Writing positions in: "<< fileAdress << std::endl;
    if (folderName != "")
    {
        fileAdress = folderName + "/" + fileAdress;
    }

   std::ofstream streamPosition(fileAdress);
    if (streamPosition.is_open())
    {   
        streamPosition << sitePositions.n_rows << std::endl;
        for (uint iSite = 0; iSite < sitePositions.n_rows ; iSite ++)
        {
            streamPosition << std::endl << siteNames[iSite] << "\t" << sitePositions(iSite, 0) << "\t" << sitePositions(iSite, 1) << "\t" << sitePositions(iSite, 2) ;
        }
        streamPosition.close();
    }
}



/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void                    FHLattice::computeMatrices()
{

    unsigned int nSites = sitePositions.n_rows;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    sameSiteMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            distanceMatrix(i, j) = computeDistance(i,j);
            if (abs(distanceMatrix(i, j))< samePointDistance)
            {
                sameSiteMatrix[i][j] = true;
            }
            if (abs(distanceMatrix(i, j) - neighbourDistance) < samePointDistance)
            {
                neighbourMatrix[i][j] = true;
            }
        }
    }


    //std::cout << distanceMatrix << std::endl;
}





FHLattice::~FHLattice()
{

}









// OTHER FUNCTIONS

/// @brief Convert a vector of vector of double into the corresponding arma::dmat. The vector'vector and the matrix are seen in row-major order (reading/writing order).
/// @param vecVec A vector of vector of double 
/// @return The arma::dmat version of the vector
arma::dmat armaMatfromVectorsVector(std::vector<std::vector<double>>& vecVec)
{
    arma::dmat mattrix = arma::dmat();
    if (vecVec.size() > 0)
    {
        if (vecVec[0].size() > 0)
        {
            mattrix = arma::dmat(vecVec.size(), vecVec[0].size(), arma::fill::zeros);
            for (uint iRow = 0; iRow < vecVec.size(); iRow ++)
            {
                for (uint iCol = 0; iCol < std::min(vecVec[iRow].size(), vecVec[0].size()); iCol ++)
                {
                    mattrix(iRow, iCol) = vecVec[iRow][iCol];
                }
            }
        }
    }

    return mattrix;
}

/// @brief Check if a vector contains an element.
/// @param w the element to look for
/// @param v the vector of possibilities
/// @return true if w is in v
template <typename T>       bool isIn(T w, std::vector<T>& v)
{

        bool isInB = false;
        unsigned int index = 0;
        while ( (!isInB) && (index < v.size()) )
        {
            if (v[index] == w)  isInB = true;
            index ++;
        }

        return isInB;
}

/// @brief Check if a vector contains an element.
/// @param w the element to look for
/// @param v the vector of possibilities
/// @param firstIndex index of the first found element
/// @return true if w is in v
template <typename T>       bool isIn(T w, std::vector<T>& v, uint& firstIndex)
{
    bool isInB = false;
    unsigned int index = 0;
    while ( (!isInB) && (index < v.size()) )
    {
        if (v[index] == w)  isInB = true;
        index ++;
    }
    firstIndex = index-1;
    return isInB;
}



