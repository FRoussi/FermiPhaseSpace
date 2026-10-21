#include "../include/FHLattice.h"

#include <algorithm>
#include <cstring>

using std::vector;
using std::array;

// Grid construction helpers, defined at the bottom of this file.
// They are implementation details of FHLattice, hence not exposed in the header.

FHLattice::FHLattice(bool talk) :
    m_Talkative(talk),
    m_LatticeType(FHLatticeType::cubic),
    m_IntParameters({1, 1, 1}),
    m_FloatParameters({1.0})
{

    init_latticeNameMap();

    arma::dmat pos = arma::dmat(2,  3, arma::fill::zeros);
    pos(1,0) = 1;

    setPositions(pos);

}


FHLattice::FHLattice(int x, bool talk) :
    m_Talkative(talk),
    m_LatticeType(FHLatticeType::cubic),
    m_IntParameters({x, 1, 1}),
    m_FloatParameters({1.0})
{

    init_latticeNameMap();

    arma::dmat pos = arma::dmat(1,  x, arma::fill::zeros);
    for (int d=0; d<x; d++)
    {
        pos[d] = (double)d;
    }

    int nPart = std::floor(x/2);
    std::vector<uint> initCond_up = std::vector<uint>(nPart, 0);
    std::vector<uint> initCond_down = std::vector<uint>(nPart, 0);
    for(int itPart = 0; itPart <nPart; itPart++)
    {
        initCond_up[itPart] = itPart * 2;
        initCond_down[itPart] = itPart * 2 + 1;
    }

    std::vector<std::vector<uint> > initCond = std::vector<std::vector<uint> >(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;


    setPositions(pos);

    m_LatticeMatrix = findBCGridBox(x, 1, 1, 1);
    computeInverseLatticeMatrix();

}


FHLattice::FHLattice(int x, int y, int z, bool talk) :
    m_Talkative(talk),
    m_LatticeType(FHLatticeType::cubic),
    m_IntParameters({x, y, z}),
    m_FloatParameters({1.0})
{

    init_latticeNameMap();

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
     
    m_LatticeMatrix = findBCGridBox(x, y, z, 1);
    computeInverseLatticeMatrix();

}



FHLattice::FHLattice(Json::Value rootSystem, std::string folderName, bool talk) :
    m_Talkative(talk),
    m_LatticeType(FHLatticeType::none),
    m_FloatParameters({1.0})
{

    init_latticeNameMap();

// Positions
    std::string positionFile = rootSystem.get("positionsFile", "").asString();

    if (folderName !="")
    {
        positionFile = folderName + "/" + positionFile;
    }

    // std::cout << positionFile << std::endl;
    readPositions(positionFile);


    std::string latticeTypeString = rootSystem.get("latticeType", "none").asString();
    auto lti = m_LatticeTypeNamesReverse.find(latticeTypeString);

    if (lti != m_LatticeTypeNamesReverse.end())
    {
        m_LatticeType = lti->second;
    }

// Periodicity
    m_IsPeriodic = rootSystem.get("periodic", false).asBool();


    Json::Value latMat = rootSystem["latticeMatrix"];
    if (latMat != Json::Value::null)
    {
        for (uint i = 0; i<3; i++)
        {
            Json::Value latRow = latMat[i]; 
            for (uint j = 0; j<3; j++)
            {
                m_LatticeMatrix(i,j) = latRow[j].asDouble();
            }
        }
        computeInverseLatticeMatrix();
    }
    

    Json::Value latRow = rootSystem["latticeZero"];
    if (latRow != Json::Value::null)
    {
        for (uint j = 0; j<3; j++)
        {
            m_LatticeZeroVector(j) = latRow[j].asDouble();
        }
    }


    if (m_IsPeriodic)
    {
        computeMatrices();
    }



//distanceMatrix.print();

}

FHLattice::FHLattice(std::string latticeTypeString, vector<double> parameters) 
{
    
    vector<array<double,3>> pos = makeBCGrid(2,1,1,1);

    init_latticeNameMap();
    m_LatticeType = FHLatticeType::none;

    auto itLatType = m_LatticeTypeNamesReverse.find(latticeTypeString);
    if (itLatType != m_LatticeTypeNamesReverse.end())
    {
        m_LatticeType = itLatType->second;
    }




switch (m_LatticeType)
{
    case FHLatticeType::none:
        break;
    case FHLatticeType::cubic:
    {

        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const int nz = static_cast<int>(parameters[2]);
        const double a = parameters[3];
        pos = makeBCGrid(nx, ny, nz, a);
        m_LatticeMatrix = findBCGridBox(nx, ny, nz, a);
        m_IntParameters = {nx, ny, nz};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "CubicGrid";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    case FHLatticeType::bulkCenteredCubic:
    {
        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const int nz = static_cast<int>(parameters[2]);
        const double a = parameters[3];
        pos = makeBCCGrid(nx, ny, nz, a);
        m_LatticeMatrix = findBCGridBox(nx, ny, nz, a);
        m_IntParameters = {nx, ny, nz};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "BCCGrid";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a/std::sqrt(3);
        break;
    }
    case FHLatticeType::faceCenteredCubic:
    {
        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const int nz = static_cast<int>(parameters[2]);
        const double a = parameters[3];
        pos = makeFCCGrid(nx, ny, nz, a);
        m_LatticeMatrix = findBCGridBox(nx, ny, nz, a);
        m_IntParameters = {nx, ny, nz};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "FCCGrid";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a/std::sqrt(2);
        break;
    }
    case FHLatticeType::honeyCombs:
    {
        int nx = static_cast<int>(parameters[0]);
        double a = parameters[1];
        pos = makeHoneyCombs(nx,a);
        m_LatticeMatrix = findHoneyCombGridBox(nx, a);
        m_IntParameters = {nx};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "HoneyCombs";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    case FHLatticeType::triangular:
    {
        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const int nz = static_cast<int>(parameters[2]);
        const double a = parameters[3];
        pos = makeTetragones(nx, ny, nz, a);
        m_LatticeMatrix = findTetragonesGridBox(nx, ny, nz, a);
        m_IntParameters = {nx, ny, nz};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "Tetragones";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    case FHLatticeType::graphene:
    {
        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const double a = parameters[2];
        pos = makeGraphene(nx, ny, a);
        m_LatticeMatrix = findGrapheneGridBox(nx, ny, a);
        m_IntParameters = {nx, ny};
        m_FloatParameters = {a};
        computeInverseLatticeMatrix();
        m_LatticeName = "Graphene";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    case FHLatticeType::moireHC:
    {
        int nCombs = static_cast<int>(parameters[0]);
        double a = parameters[1];
        double b = parameters[2];
        double theta = parameters[3];
        pos = makeMoireHoneyCombs(nCombs,  a, b, theta);
        m_LatticeMatrix = findHoneyCombGridBox(nCombs, a);
        m_IntParameters = {nCombs};
        m_FloatParameters = {a, b, theta};
        computeInverseLatticeMatrix();
        m_LatticeName = "MoireHC";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    case FHLatticeType::moire:
    {
        const int nx = static_cast<int>(parameters[0]);
        const int ny = static_cast<int>(parameters[1]);
        const double a = parameters[2];
        const double b = parameters[3];
        const double theta = parameters[4];
        pos = makeMoireGraphene(nx, ny, a, b, theta);
        m_LatticeMatrix = findGrapheneGridBox(nx, ny, a);
        m_IntParameters = {nx, ny};
        m_FloatParameters = {a, b, theta};
        computeInverseLatticeMatrix();
        m_LatticeName = "Moire";
        m_SamePointDistance = 1.0E-6 * a;
        m_NeighbourDistance = a;
        break;
    }
    default:
        break;
    }


    setPositions(armaMatfromVectorsPosition(pos), "C");

    computeMatrices();
}


FHLattice::FHLattice(FHLatticeType latType, vector<double> parameters) :
    m_LatticeType(latType)
{

    vector<array<double,3>> pos = makeBCGrid(2,1,1,1);
    init_latticeNameMap();

switch (latType)
{
case FHLatticeType::none:
    break;
    case FHLatticeType::cubic:
        pos = makeBCGrid(static_cast<int>(parameters[0]), static_cast<int>(parameters[1]), static_cast<int>(parameters[2]), parameters[3]);
        m_LatticeName = "CubicGrid";
    break;
    case FHLatticeType::bulkCenteredCubic:
        pos = makeBCCGrid(static_cast<int>(parameters[0]),parameters[1], parameters[2], 1.);
        m_LatticeName = "BCCGrid";
    break;
    case FHLatticeType::faceCenteredCubic:
        pos = makeFCCGrid(static_cast<int>(parameters[0]),parameters[1], parameters[2], 1.);
        m_LatticeName = "FCCGrid";
    break;
    case FHLatticeType::graphene:
        pos = makeHoneyCombs(static_cast<int>(parameters[0]),parameters[1]);
        m_LatticeName = "HoneyCombs";
    break;

    default:
    break;
}

    setPositions(armaMatfromVectorsPosition(pos), "C");

    computeMatrices();

}


FHLattice::FHLattice(const FHLattice& fhLattice) :
    m_Talkative(fhLattice.m_Talkative),
    m_LatticeType(fhLattice.m_LatticeType),
    m_IntParameters(fhLattice.m_IntParameters),
    m_FloatParameters(fhLattice.m_FloatParameters),
    m_NeighbourDistance(fhLattice.m_NeighbourDistance),
    m_SamePointDistance(fhLattice.m_SamePointDistance),
    m_IsPeriodic(fhLattice.m_IsPeriodic),
    m_LatticeMatrix(fhLattice.m_LatticeMatrix),
    m_LatticeZeroVector(fhLattice.m_LatticeZeroVector),
    m_LatticeName(fhLattice.m_LatticeName)
{
    init_latticeNameMap();
    m_SitePositions = arma::dmat(fhLattice.getSitePositions());

    m_SiteTypes = fhLattice.getSiteTypes();
    m_SiteTypeNames = fhLattice.getSiteTypeNames();
    m_SiteNames = fhLattice.getSiteNames();

    if (arma::det(m_LatticeMatrix) != 0)
    {
        computeInverseLatticeMatrix();
    }

    m_LatticeType = fhLattice.getLatticeType();

    computeMatrices();

}


/// @brief getter - number of sites
/// @return number of sites
uint                                FHLattice::getnSites() const
{
    return m_SitePositions.n_rows;
}

/// @brief getter - number of dimension (we always work in dim 3 though)
/// @return number of dimensions (if not 3, weird)
uint                                FHLattice::getnDimension() const
{
    return m_SitePositions.n_cols;
}

/// @brief getter - position row vector of site 'iSite'
/// @param iSite a site number 
/// @return row vector of coordinates
arma::drowvec                               FHLattice::getSitePosition(uint iSite) const 
{
    return m_SitePositions.row(iSite);
}

const arma::dmat&                         FHLattice::getSitePositions() const
{
    return m_SitePositions;
}

/// @brief getter - coordinate of a site
/// @param dim coordinate dimension
/// @param iSite site number
/// @return coordinate of site 'iSite'in dimension 'dim'
double                                      FHLattice::getCoordinate(uint dim, uint iSite) const 
{ 
    if ((dim < m_SitePositions.n_cols) && (iSite < m_SitePositions.n_rows))
    {
        return m_SitePositions(iSite, dim);
    } 
    else
    {
        return 0;
    }
}

/// @brief getter - relative position row vector of site 'iSite'
/// @param iSite a site number 
/// @return row vector of relative coordinates
arma::drowvec                       FHLattice::getLatticeSitePosition(uint iSite) const
{
    //std::cout << sitePositions.row(iSite) << std::endl;
    //std::cout << inverseLatticeMatrix << std::endl;
   // std::cout << latticeZeroVector << std::endl;
    arma::drowvec pos = m_SitePositions.row(iSite) - m_LatticeZeroVector;
   // std::cout << inverseLatticeMatrix * pos.t() << std::endl;
    return (m_InverseLatticeMatrix * pos.t()).t();
}

double                              FHLattice::getLatticeCoordinate(uint dim, uint iSite) const
{
    if ((dim < m_SitePositions.n_cols) && (iSite < m_SitePositions.n_rows))
    {
        arma::drowvec v = getLatticeSitePosition(iSite);
        return v[dim];
    } 
    else
    {
        return 0;
    }
}

FHLatticeType                       FHLattice::getLatticeType() const
{
    return m_LatticeType;
}

bool                                        FHLattice::sitesHaveNames() const
{
    return (m_SiteTypes.n_rows == m_SitePositions.n_rows);
}

uint                                FHLattice::getSiteNameNumber() const
{
    return (uint)m_SiteNames.size();
}

std::string                                 FHLattice::getSiteName(uint iSite) const
{
    //return siteTypeNames[siteTypes[site]];
    return m_SiteNames[iSite];
}

vector<std::string>            FHLattice::getSiteNames() const
{
    return m_SiteNames;
}

arma::uvec                          FHLattice::getSiteTypes() const
{
    return m_SiteTypes;
}

uint                                FHLattice::getSiteType(uint iSite) const
{
    return m_SiteTypes[iSite];
}

vector<std::string>                         FHLattice::getSiteTypeNames() const {    return m_SiteTypeNames;}
bool                                        FHLattice::isPeriodic() const     {    return m_IsPeriodic;}
const arma::dmat&                           FHLattice::getLatticeMatrix() const {    return m_LatticeMatrix;}
const arma::drowvec&                        FHLattice::getLatticeZeroVector() const   {    return m_LatticeZeroVector;}
double                                      FHLattice::getNeighbourDistance() const
{
    return m_NeighbourDistance;
}
double                                   FHLattice::getSamePointDistance() const
{
    return m_SamePointDistance;
}

double                                      FHLattice::computeDistance(uint iSite, uint jSite) const 
{
    arma::drowvec dVec = m_SitePositions.row(iSite) - m_SitePositions.row(jSite); 
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
                    arma::drowvec dVecSym = dVec + xSym * m_LatticeMatrix.col(0).t() + ySym * m_LatticeMatrix.col(1).t() + zSym * m_LatticeMatrix.col(2).t(); //<- new vector with one of the symmetric lattice
                    d = std::min(d, arma::norm(dVecSym));                   //<- we take the minimum of them all
                    //std::cout << xSym << "  " << ySym << "  " << zSym << "  " << arma::norm(dVecSym) << "   " << d << std::endl;
                }
            }
        }
            
    }

    return d;
}


double                                      FHLattice::evaluateNearestNeighbor() const
{
    double minNeib = 0;
    if (m_DistanceMatrix.n_cols > 1)
    {
        minNeib = 1000;
        for (uint iSite = 0; iSite < m_DistanceMatrix.n_rows; iSite++)
        {
            for (uint jSite = iSite+1; jSite < m_DistanceMatrix.n_cols; jSite++)
            {
                if (m_DistanceMatrix.at(iSite, jSite) >  m_SamePointDistance)
                {
                    minNeib = std::min(minNeib, m_DistanceMatrix.at(iSite, jSite));
                }
            }
        }
    }
    return minNeib;
}


double                                      FHLattice::getDistance(uint i, uint j)    const   {   return m_DistanceMatrix(i,j);}
bool                                        FHLattice::areNeighbour(uint i, uint j)   const   {   return m_NeighbourMatrix[i][j];}
bool                                        FHLattice::areSameSite(uint i, uint j)      const   {   return m_SameSiteMatrix[i][j];}
std::string                                 FHLattice::getLatticeName()                 const   {   return m_LatticeName;}

const arma::dmat&                           FHLattice::getDistanceMatrix()              const   {   return m_DistanceMatrix;}
const std::vector<std::vector<bool>>&       FHLattice::getNeighbourMatrix()             const   {   return m_NeighbourMatrix;}


// SETTERS

void                    FHLattice::setPositions(const arma::dmat& positions, std::string sitesName)
{
    m_SitePositions = positions;
    m_SiteNames = vector<std::string>(m_SitePositions.n_rows, sitesName);
    computeMatrices();
}

void                    FHLattice::setPositions(const arma::dmat& positions, vector<std::string> sN)
{
    m_SitePositions = positions;
    m_SiteNames = sN;
    computeMatrices();
}


bool                    FHLattice::setSiteTypes(arma::uvec sTypes)
{
    m_SiteTypes = sTypes;
    return sitesHaveNames();
}

void                    FHLattice::setSiteTypeNames(vector<std::string> typeNames)
{
    m_SiteTypeNames = typeNames;
}


void                    FHLattice::setSymmetry(arma::dmat symVectors)
{
    if ((symVectors.n_cols != 3)&& (symVectors.n_rows != 3))
    {
        std::cout << "We want a 3x3 lattice matrix" << std::endl;
    }
    else
    {
        m_IsPeriodic = true;
        m_LatticeMatrix = symVectors;
        m_InverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
        bool isInv = arma::inv(m_InverseLatticeMatrix, m_LatticeMatrix);
        if (!isInv)
        {
            std::cout << "Caution : the lattice matrix is singular." << std::endl;
        }
    }

    if (m_LatticeZeroVector.n_rows !=3)
    {
        m_LatticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    }

}

void                                FHLattice::setLatticeZero(arma::drowvec latticeZ){    m_LatticeZeroVector = latticeZ;}
void                                FHLattice::setNeighbourDistance(double neigbD){    m_NeighbourDistance = neigbD;}
void                                FHLattice::setSamePointDistance(double samePointD){    m_SamePointDistance = samePointD;}
void                                FHLattice::setLatticeName(std::string newName){ m_LatticeName = newName;}

//// @brief Computes the inverse of the lattice matrix, if it exists.
void                                FHLattice::computeInverseLatticeMatrix()
{
    m_InverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    double rcond = 0;
    bool success = arma::inv(m_InverseLatticeMatrix, rcond, m_LatticeMatrix);
    if (!success)
    {
        std::cout << "! Problem ! Lattice matrice not inversible." << std::endl;
    }

}

/// @brief Reads the positions of the lattice sites from a .xyz file.
/// @param xyzFile The path to the file containing the site positions.
void                                FHLattice::readPositions(std::string xyzFile)
{

    uint nSites = 0;
    vector<array<double,3>> positions3;

    m_SiteNames = vector<std::string>();
    std::string line;
    std::string coord;
     

    std::ifstream posFile(xyzFile);
    if (posFile.is_open())
    {
        // Site number :
        std::getline(posFile, line);

        nSites = std::stoi(line);
        // Comment :
        std::getline(posFile, line);
        // Sites :
        while(std::getline(posFile, line))
        {
            std::stringstream ss(line);
            // Get the site name
            ss >> coord;
            m_SiteNames.push_back(coord);
            // Get its coordinates :
            vector<double> coordinates;
            while (ss >> coord)//(std::getline (ss, coord, '\t')) 
            {
                double d = std::stod(coord);
                coordinates.push_back(d);
            }                      
            coordinates.resize(3);

            positions3.push_back({coordinates[0], coordinates[1], coordinates[2]});
        }

        posFile.close();
    }
    else
    {
        std::cout << "Problem with the site position file " << xyzFile << std::endl;
    }


    if (nSites == positions3.size())
    {
        if (m_Talkative)
        {
            std::cout << "   -- Position file read, sites :   " << nSites << std::endl;
        }
        setPositions(armaMatfromVectorsPosition(positions3));
    }
    else
    {
        std::cout << "Problem in the xyz file : " << nSites << " sites announced, " << positions3.size() << " found." << std::endl;
    }



    computeMatrices();

}


/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void                 FHLattice::writePositions(std::string xyzFile, std::string folderName)
{
        // Writing of everything on the file :
    m_LatticeName = xyzFile;
    std::string fileAdress =  xyzFile + ".xyz";
   // std::cout << "Writing positions in: "<< fileAdress << std::endl;
    if (folderName != "")
    {
        fileAdress = folderName + "/" + fileAdress;
    }

   std::ofstream streamPosition(fileAdress);
    if (streamPosition.is_open())
    {   
        
        streamPosition << m_SitePositions.n_rows << std::endl;
        for (uint iSite = 0; iSite < m_SitePositions.n_rows ; iSite ++)
        {
            streamPosition << std::endl << m_SiteNames[iSite] << "\t" << m_SitePositions(iSite, 0) << "\t" << m_SitePositions(iSite, 1) << "\t" << m_SitePositions(iSite, 2) ;
        }
        streamPosition.close();
    }
}

/// @brief Writes the lattice information to a JSON object.
/// @return Json::Value The JSON object containing the lattice information.
Json::Value         FHLattice::writeFHLatticeJson() const
{
        Json::Value rootLattice;

        rootLattice["latticeType"] = m_LatticeTypeNames.at(m_LatticeType);


        rootLattice["positionsFile"] = m_LatticeName + ".xyz";

        // Periodicity
        rootLattice["periodic"] = m_IsPeriodic;
        Json::Value vecJson(Json::arrayValue);
        Json::Value matJson(Json::arrayValue);
        for (uint i = 0; i<3; i++)  vecJson[i] = m_LatticeZeroVector(i);
        for (uint i = 0; i<3; i++)
        {
            Json::Value vJson(Json::arrayValue);
            for (uint j = 0; j<3; j++)
            {
                vJson.append(m_LatticeMatrix.at(i,j));
            }
            matJson.append(vJson);
        }  
        rootLattice["latticeMatrix"] = matJson;
        rootLattice["latticeZero"] = vecJson;
        

        return rootLattice;

}




/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void                    FHLattice::computeMatrices()
{

//std::cout << "compute matrices, smpD =" << samePointDistance << ", nD =" << neighbourDistance << std::endl;

    uint nSites = m_SitePositions.n_rows;

    m_DistanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    m_NeighbourMatrix = vector<vector<bool>>(nSites,  vector<bool>(nSites, false));
    m_SameSiteMatrix = vector<vector<bool>>(nSites,  vector<bool>(nSites, false));


    // distance Matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            m_DistanceMatrix(i, j) = computeDistance(i,j);
        }
    }
    // same site matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            if (std::norm(m_DistanceMatrix(i, j))< m_SamePointDistance)
            {
                m_SameSiteMatrix[i][j] = true;
            }
        }
    }
    // neighbour matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            //std::cout << "i" << i << ", j=" << j << " , d=" << m_DistanceMatrix(i, j) << "|d-n|=" << std::norm(m_DistanceMatrix(i, j) - m_NeighbourDistance) <<  std::endl;
            if (std::norm(m_DistanceMatrix(i, j) - m_NeighbourDistance) < m_SamePointDistance)
            {
                m_NeighbourMatrix[i][j] = true;
            }
        }
    }




   // distanceMatrix.print();
}


/// @brief Initializes the mapping between lattice types and their string representations.
void FHLattice::init_latticeNameMap()
{
    m_LatticeTypeNames = std::map<FHLatticeType, std::string>();
    m_LatticeTypeNames[FHLatticeType::none] = "none";
    m_LatticeTypeNames[FHLatticeType::cubic] = "grid";
    m_LatticeTypeNames[FHLatticeType::bulkCenteredCubic] = "bcc";
    m_LatticeTypeNames[FHLatticeType::faceCenteredCubic] = "fcc";
    m_LatticeTypeNames[FHLatticeType::triangular] = "triangular";
    m_LatticeTypeNames[FHLatticeType::honeyCombs] = "honeyCombs";
    m_LatticeTypeNames[FHLatticeType::graphene] = "graphene";
    m_LatticeTypeNames[FHLatticeType::moireHC] = "moireHC";
    m_LatticeTypeNames[FHLatticeType::moire] = "moire";
    m_LatticeTypeNames[FHLatticeType::kagome] = "kagome";

    m_LatticeTypeNamesReverse = std::map<std::string, FHLatticeType>();

    for (const auto& entry : m_LatticeTypeNames)
    {
        m_LatticeTypeNamesReverse.emplace(entry.second, entry.first);
    }
}
 











// OTHER FUNCTIONS

// I - Grids

/// @brief Makes a simple 3D-grid of sites 
/// @param nx number of sites in the X dimension
/// @param ny number of sites in the Y dimension
/// @param nz number of sites in the Z dimension
/// @param a distance between nearest neighbours
/// @return a vector of site positions
vector<array<double,3>> makeBCGrid(uint nx, uint ny, uint nz, double a)
{
    
    uint nPoints = nx*ny*nz;
    //std::cout <<"makegrid" << nx << " " << ny << " " << nz << " " << a << std::endl; 
    vector<array<double,3>> pos = vector<array<double,3>>(nPoints);

    uint counter = 0;
    for (uint z = 0; z < nz; z++)
    {
        for (uint y = 0; y < ny; y++)
        {
            for (uint x = 0; x < nx; x++)
            {
                //std::cout <<"makegrid in" << x << " " << y << " " << z << " " << std::endl; 
                array<double,3> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
            }
        }
    }

    return pos;
}

arma::dmat findBCGridBox(uint nx, uint ny, uint nz, double a)
{
    arma::dvec diag = {nx * a, ny * a, nz * a};
    return arma::diagmat(diag);
}

/// @brief Makes a Body-centered Cubic grid of sites 
/// @param nx number of cubes in the X dimension
/// @param ny number of cubes in the Y dimension
/// @param nz number of cubes in the Z dimension
/// @param a side of a cube, so the distance between nearest neighbours is a/sqrt(3)
/// @return a vector of site positions
vector<array<double,3>> makeBCCGrid(uint nx, uint ny, uint nz, double a)
{
    
    uint nPoints = 2*(nx+1)*(ny+1)*(nz+1) + 2 *((nx+1)*(ny+1) + (nx+1)*(nz+1) + (ny+1)*(nz+1)) +1;
    vector<array<double,3>> pos = vector<array<double,3>>(nPoints);


    uint counter = 0;
    for (uint z = 0; z <= nz; z++)
    {
        for (uint y = 0; y <= ny; y++)
        {
            for (uint x = 0; x <= nx; x++)
            {
                array<double,3> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
                if ((x != nx) && (y != ny) && (z != nz) )
                {
                    newPos = {a*(x + 0.5), a*(y + 0.5), a*(z+0.5)};
                    pos[counter] = newPos;
                    counter +=1;
                } 
            }
        }
    }

    pos.resize(counter);

    return pos;
}


/// @brief Makes a Face-centered Cubic grid of sites 
/// @param nx number of cubes in the X dimension
/// @param ny number of cubes in the Y dimension
/// @param nz number of cubes in the Z dimension
/// @param a side of a cube, so the distance between nearest neighbours is a/sqrt(2)
/// @return a vector of site positions
vector<array<double,3>> makeFCCGrid(uint nx, uint ny, uint nz, double a)
{
    
    uint nPoints = 4*(nx+1)*(ny+1)*(nz+1) + 2 *((nx+1)*(ny+1) + (nx+1)*(nz+1) + (ny+1)*(nz+1)) +1;
    vector<array<double,3>> pos = vector<array<double,3>>(nPoints);

    uint counter = 0;
    for (uint z = 0; z <= nz; z++)
    {
        for (uint y = 0; y <= ny; y++)
        {
            for (uint x = 0; x <= nx; x++)
            {
                array<double,3> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
                if ((y != ny) && (z != nz) )
                {
                    newPos = {a*x , a*(y + 0.5), a*(z + 0.5)};
                    pos[counter] = newPos;
                    counter +=1;
                } 
                if ((x != nx) && (z != nz) )
                {
                    newPos = {a*(x + 0.5) , a*y, a*(z + 0.5)};
                    pos[counter] = newPos;
                    counter +=1;
                } 
                if ((x != nx) && (y != ny) )
                {
                    newPos = {a*(x + 0.5) , a*(y + 0.5), a*z};
                    pos[counter] = newPos;
                    counter +=1;
                } 
            }
        }
    }

    pos.resize(counter);

    return pos;
}




/// @brief Makes a graphene-shaped site vector
/// @param nCombs number of combds in the graphene sheet
/// @param a distance between nearest neighbours
/// @return a vector of site positions
vector<array<double,3>> makeHoneyCombs(uint nCombs, double a)
{

    double sqrt3 = std::sqrt(3);
    uint nOfPoints = 0;
    vector<array<double,3>> pos = vector<array<double,3>>();

    vector<array<double,3>> pos0 = { {-sqrt3/2 *a, -0.5 *a, 0},
                                {-sqrt3/2 *a, 0.5 *a, 0},
                                {0   *a,      1   *a, 0},
                                {0   *a,      -1  *a, 0},
                                {sqrt3/2 *a, -0.5 *a, 0},
                                {sqrt3/2 *a, 0.5  *a, 0}};
    vector<array<double,3>> posNew = pos0;


    double iVec[2] = {sqrt3 *a, 0 *a};
    double jVec[2] = {sqrt3/2 *a, 1.5 *a};

    uint basisIndex[2] = {0,0};


    for (uint bv = 0; bv < nCombs; bv++)
    {
        for (uint p = 0; p<6; p++)
        {
            posNew[p][0] = pos0[p][0] + iVec[0] * basisIndex[0] + jVec[0] * basisIndex[1];
            posNew[p][1] = pos0[p][1] + iVec[1] * basisIndex[0] + jVec[1] * basisIndex[1];
        }
        pos.insert(pos.end(), posNew.begin(), posNew.end());

        getNextIndex(basisIndex[0], basisIndex[1]);

        nOfPoints +=6;
    }

    // We remove the doublons :
    double samePointThreshold = 0.0001 * a;


    for (uint i1 = 0; i1< pos.size(); i1++)         //<- we look for all points equal to positions(i1)
    {
        uint nDoublonsFound = 0;                    //<- There can't be more than three times the same point
         
        uint i2 = i1+1;
        while ((i2<pos.size()) && (nDoublonsFound < 2))    //<- as long as there are points and more doublons to find
        {
            double d = 0;
            for (uint iDim = 0; iDim<pos[i1].size(); iDim++)
            {
                d += std::abs(pos[i1][iDim] - pos[i2][iDim]);
            }
            if (abs(d) < samePointThreshold)                //<- if i1 and i2 are the same point :
            {
                pos.erase(pos.begin() + i2);                //<- We remove the point i2
                nDoublonsFound ++;                          //<- We say we found one more doublon
            }
            else
            {
                i2 ++;                                      //<- else, we look at the next point
            }
                
        }

    }
    
    pos.shrink_to_fit();
    return pos;
}

arma::dmat findHoneyCombGridBox(uint nCombs, double a)
{
    double sqrt3 = std::sqrt(3);
    arma::dmat part2d = {{sqrt3 *a, sqrt3/2 *a}, {0, 1.5 * a}};

    double sqrtnCombs = std::sqrt((double)nCombs);
    uint fsnc = (uint)std::floor(sqrtnCombs);
    uint csnc = (uint)std::ceil(sqrtnCombs);

    uint nx = csnc;
    uint ny = fsnc + (uint)((nCombs - fsnc*fsnc)> fsnc);


    part2d.col(0) = part2d.col(0) * nx;
    part2d.col(1) = part2d.col(1) * ny; 

    arma::dmat box = arma::dmat(3,3,arma::fill::zeros);
    box.submat(0,0,1,1) = part2d;
    box(2,2) = 10*a;


    return box;
}


vector<array<double,3>> makeGraphene(int nx, int ny, double a)
{
    double sqrt3 = std::sqrt(3);
    uint nOfPoints =0;
    vector<array<double,3>> pos = vector<array<double,3>>();

    vector<array<double,3>> pos0 = { {0, -0.5 *a , 0},
                                    {0, 0.5 *a  , 0}};


    double iVec[2] = {sqrt3 *a, 0 *a};
    double jVec[2] = {sqrt3/2 *a, 1.5 *a};


    for (int j = 0; j < ny; j++)
    {
        vector<array<double,3>> posNewY = pos0;
        posNewY[0][0] += jVec[0] * j;
        posNewY[0][1] += jVec[1] * j;
        posNewY[1][0] += jVec[0] * j;
        posNewY[1][1] += jVec[1] * j;
        for (int i = 0; i<nx; i++)
        {
            vector<array<double,3>> posNew = posNewY;
            posNew[0][0] += iVec[0] * i;
            posNew[0][1] += iVec[1] * i;
            posNew[1][0] += iVec[0] * i;
            posNew[1][1] += iVec[1] * i;
            //std::cout << "graph vec " << i << " " << j << " v1(1)=" << posNew[0][1] << std::endl;
            pos.push_back(posNew[0]);
            pos.push_back(posNew[1]);
            nOfPoints +=2;
        }
    }

    pos.shrink_to_fit();
    return pos;
}

arma::dmat findGrapheneGridBox(int nx, int ny, double a)
{
    double sqrt3 = std::sqrt(3);
    arma::dmat mat3 = { {sqrt3 * a * nx,    sqrt3 / 2 * a * ny  , 0     },
                        {0             ,    1.5 * a * ny        , 0     },
                        {0             ,    0                   , 10*a  }};
    
    return mat3;
}

vector<array<double,3>> makeTetragones(int nx, int ny, int nz, double a)
{
    uint nOfPoints = 0;
    double sqrt3 = std::sqrt(3);
    vector<array<double,3>> pos = vector<array<double,3>>();
    double iVec[3] = {a     , 0 , 0};
    double jVec[3] = {1./2. *a, sqrt3 / 2. * a     , 0};
    double kVec[3] = {1./2. *a, 1. / (2*sqrt3) * a , sqrt(2) / sqrt3 * a};

    for (int k = 0; k<nz; k++)
    {
        for (int j = 0;j<ny; j++)
        {
            for (int i = 0; i<nx; i++)
            {
                array<double,3> newVec = {   i * iVec[0] + j * jVec[0] + k * kVec[0],   
                                            i * iVec[1] + j * jVec[1] + k * kVec[1],           
                                            i * iVec[2] + j * jVec[2] + k * kVec[2]};
                pos.push_back(newVec);
                nOfPoints +=1;
            } 
        } 
    }
    return pos;
}

arma::dmat findTetragonesGridBox(int nx,int ny, int nz, double a)
{
    double sqrt3 = std::sqrt(3);
    arma::dmat mat3 = { {a * nx, a / 2 * ny         , a / 2 * nz },
                        {0     , sqrt3 * a / 2 * ny , a / (2*sqrt3) * nz },
                        {0     , 0                  , sqrt(2) * a / sqrt3 * nz }};
    
    return mat3;
}

vector<array<double,3>> makeMoireHoneyCombs(uint nCombs, double a, double b, double theta)
{
    double samePointThreshold = 0.001 * a;
    double sqrt3 = std::sqrt(3);
    uint nOfPoints = 0;
    uint nLayers = 2;
    vector<array<double,3>> pos;

    arma::dmat pos0 =  { {-sqrt3/2   , -sqrt3/2  , 0, 0  , sqrt3/2   , sqrt3/2   },
                         {-0.5       , 0.5       , 1, -1 , -0.5      , 0.5       },
                         {0          , 0         , 0, 0  , 0         , 0         }};
    pos0 = a * pos0;


    arma::dmat graphVec = { { sqrt3 * a , sqrt3 * a / 2.},
                            {0          , 1.5},
                            {0          , 0}};

    
                                        
    double   thetaDeg = theta * arma::datum::pi/180;
    

    for (uint iLayer = 0; iLayer < nLayers; iLayer++)       //<- For each layer
    {

        uint nOfPointsLayer = 0;
        arma::dmat positionsLayer = arma::dmat(3,6*nCombs); 

        double cosT = std::cos(iLayer * thetaDeg);
        double sinT = std::sin(iLayer * thetaDeg);
        arma::dmat rotationLayer = {    {cosT   , -sinT , 0},
                                        {sinT   , cosT  , 0},
                                        {0      ,   0   , 1}};

        arma::dcolvec shiftLayer = {0,0,iLayer *b};

        arma::dmat graphVecLayer = rotationLayer * graphVec;
        arma::dmat pos0Layer = rotationLayer * pos0;

 

        int n = 0;
        int m = 0;
        int l = 0;
        arma::dcolvec basisIndex = {(double)(m-l), (double)l};
        

        for (uint iComb = 0; iComb < nCombs; iComb ++)                  //<- we build the combs
        {
            arma::dcolvec vectorComb = graphVecLayer * basisIndex + shiftLayer;
            arma::dmat posNew = pos0Layer.each_col() + vectorComb;
            positionsLayer.cols(nOfPointsLayer, nOfPointsLayer + 5) = posNew;
            nOfPointsLayer += 6;

            getNextIndexMoireHoneyCombs(n, m, l);
            basisIndex = {(double)(m-l), (double)l};
        }
    

        uint i1 = 1;
        while (i1 < positionsLayer.n_cols)          //<- We remove the doublons
        {
            uint nDoublonFound = 0;
            uint i2 = i1+1;
            while ((i2<positionsLayer.n_cols) && (nDoublonFound < 2))
            {
                if (arma::norm(positionsLayer.col(i1) - positionsLayer.col(i2), 1) < samePointThreshold)    //<- if 2 poits are the same
                {
                    positionsLayer.shed_col(i2);                                                            //<- we remove the second one
                    nDoublonFound += 1;                                                                     //<- we cannot have more than a 3-times doublons
                }
                else
                {
                    i2 += 1;
                }

            }

            i1 += 1;
        }
        nOfPoints += positionsLayer.n_cols;
        
        for (uint iPos = 0; iPos < positionsLayer.n_cols ; iPos++)
        {
            array<double,3> iVec = {positionsLayer(0,iPos), positionsLayer(1,iPos), positionsLayer(2,iPos)};
            pos.push_back(iVec);
        }
    }

    return pos;
}


vector<array<double,3>> makeMoireGraphene(int nx, int ny, double a, double b, double theta)
{
    vector<array<double,3>> pos;
    uint nLayers = 2;
    double   thetaDeg = theta * arma::datum::pi/180;

    for (uint iLayer = 0; iLayer < nLayers; iLayer++)
    {

        double cosT = std::cos(iLayer * thetaDeg);
        double sinT = std::sin(iLayer * thetaDeg);
        arma::dmat rotationLayer = {    {cosT   , -sinT , 0},
                                        {sinT   , cosT  , 0},
                                        {0      ,   0   , 1}};

        vector<array<double,3>> posLayer = makeGraphene(nx, ny, a);

        arma::dmat graphBox = findGrapheneGridBox(nx, ny, a);
        arma::dcolvec midBox = (graphBox.col(0) + graphBox.col(1))/2.;
        
        for (auto& point : posLayer)
        {
            arma::dcolvec aPos(point.data(), 3);
            aPos -= midBox;
            aPos = rotationLayer * aPos;
            aPos(2) = iLayer * b;

            std::memcpy(point.data(), aPos.memptr(), 3 * sizeof(double));
        }

        pos.insert(pos.end(), posLayer.begin(), posLayer.end());
    }


    return pos;
}

/// @brief Convert a vector of 3d-positions into the corresponding arma::dmat (one row per position).
/// @param vecPos the vector of positions
/// @return The arma::dmat version of the vector
arma::dmat armaMatfromVectorsPosition(const vector<array<double,3>>& vecPos)
{
    arma::dmat positions(vecPos.size(), 3, arma::fill::zeros);
    for (uint iRow = 0; iRow < vecPos.size(); iRow++)
    {
        for (uint iCol = 0; iCol < 3; iCol++)
        {
            positions(iRow, iCol) = vecPos[iRow][iCol];
        }
    }

    return positions;
}

/// @brief The suite of integer doubles S defined by  S(t+1) = getNextIndex(S(t))
/// and starting by (0,0) paves the positive plan in an ordely fashion. The
/// aim is for [S(0) .. S(t)] to stay "as square as possible" for every index t.  
void getNextIndex(uint& i, uint& j)
{
    if (i==j)
    {
        i=i+1;
        j = 0;
    }
    else if (i==(j+1))
    {
        j = i;
        i = 0;
    }
    else if (i>j)
    {
        j = j+1;
    }
    else
    {
        i=i+1;
    }

}

void getNextIndexMoireHoneyCombs(int &n, int &m, int &l)
{

    if (l == std::min(n, m+n))
    {
        if (m==(-n))
        {
            n = n+1;
            m = n;
            l = 0;
        }
        else
        {
            m = m-1;
            l = std::max(-n, m-n); 
        }
    }
    else
    {
        if (std::abs(m) == n)
        {
            l = l+1;
        }
        else
        {
            l = std::min(n, m+n);
        }
    }

}
