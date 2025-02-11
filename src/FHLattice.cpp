#include "../include/FHLattice.h"


FHLattice::FHLattice(bool talk)
{

    talkative = talk;

    latticeType = FHLatticeType::cubic;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::eye);
    latticeZeroVector = arma::drowvec(3, arma::fill::eye);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    init_latticeNameMap();

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
    init_latticeNameMap();

    arma::dmat pos = arma::dmat(1,  x, arma::fill::zeros);
    for (int d=0; d<x; d++)
    {
        pos[d] = (double)d;
    }

    int nPart = std::floor(x/2);
    vector<unsigned int> initCond_up =  vector<unsigned int>(nPart,0);
    vector<unsigned int> initCond_down =  vector<unsigned int>(nPart,0);
    for(int itPart = 0; itPart <nPart; itPart++)
    {
        initCond_up[itPart] = itPart * 2;
        initCond_down[itPart] = itPart * 2 + 1;
    }

    vector<vector<unsigned int>> initCond = vector<vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;


    setPositions(pos);
    latticeMatrix = findBCGridBox(x, 1, 1, 1);
    computeInverseLatticeMatrix();

}


FHLattice::FHLattice(int x, int y, int z, bool talk)
{

    talkative = talk;

   latticeType = FHLatticeType::cubic;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
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
     
    latticeMatrix = findBCGridBox(x, y, z, 1);
    computeInverseLatticeMatrix();

}



FHLattice::FHLattice(Json::Value rootSystem, std::string folderName, bool talk)
{


    talkative = talk;

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    init_latticeNameMap();

// Positions
    std::string positionFile = rootSystem.get("positionsFile", "").asString();

    if (folderName !="")
    {
        positionFile = folderName + "/" + positionFile;
    }

    // std::cout << positionFile << std::endl;
    readPositions(positionFile);


    latticeType = FHLatticeType::none;

    std::string latticeTypeString = rootSystem.get("latticeType", "none").asString();
    std::map< std::string, FHLatticeType>::iterator lti = latticeTypeNamesReverse.find(latticeTypeString);

    if (lti != latticeTypeNamesReverse.end())
    {
        latticeType = lti->second;
    }

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
        computeInverseLatticeMatrix();
    }
    

    Json::Value latRow = rootSystem["latticeZero"];
    if (latRow != Json::Value::null)
    {
        for (uint j = 0; j<3; j++)
        {
            latticeZeroVector(j) = latRow[j].asDouble();
        }
    }

    if (isPeriodic())
    {
        computeMatrices();
    }



}

FHLattice::FHLattice(std::string latticeTypeString, vector<double> parameters)
{
    
    vector<vector<double>> pos = makeBCGrid(2,1,1,1);

    isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    init_latticeNameMap();
    latticeType = FHLatticeType::none;

    std::map<std::string, FHLatticeType>::iterator itLatType = latticeTypeNamesReverse.find(latticeTypeString);
    if (itLatType != latticeTypeNamesReverse.end())
    {
        latticeType = itLatType->second;
    }


uint nx, ny, nz;
double a;

switch (latticeType)
{
    case none:
        break;
    case cubic:
    {
        nx = static_cast<int>(parameters[0]);
        ny = static_cast<int>(parameters[1]);
        nz = static_cast<int>(parameters[2]);
        a = (parameters[3]);
        pos = makeBCGrid(nx, ny, nz, a);
        latticeMatrix = findBCGridBox(nx, ny, nz, a);
        computeInverseLatticeMatrix();
        latticeName = "CubicGrid";
        break;
    }
    case bulkCenteredCubic:
    {
        nx = static_cast<int>(parameters[0]);
        ny = static_cast<int>(parameters[1]);
        nz = static_cast<int>(parameters[2]);
        a = (parameters[3]);
        pos = makeBCCGrid(nx, ny, nz, a);
        latticeMatrix = findBCGridBox(nx, ny, nz, a);
        computeInverseLatticeMatrix();
        latticeName = "BCCGrid";
        break;
    }
    case faceCenteredCubic:
    {
        nx = static_cast<int>(parameters[0]);
        ny = static_cast<int>(parameters[1]);
        nz = static_cast<int>(parameters[2]);
        a = (parameters[3]);
        pos = makeFCCGrid(nx, ny, nz, a);
        latticeMatrix = findBCGridBox(nx, ny, nz, a);
        computeInverseLatticeMatrix();
        latticeName = "FCCGrid";
        break;
    }
    case graphene:
    {
        nx = static_cast<int>(parameters[0]);
        a = (parameters[1]);
        pos = makeHoneyCombs(static_cast<int>(parameters[0]),parameters[1]);
        latticeMatrix = findHoneyCombGridBox(nx, a);
        computeInverseLatticeMatrix();
        latticeName = "Graphene";
        break;
    }
    default:
        break;
    }


    setPositions(armaMatfromVectorsVector(pos), "C");

    computeMatrices();
}

FHLattice::FHLattice(FHLatticeType latType, vector<double> parameters)
{

 vector<vector<double>> pos = makeBCGrid(2,1,1,1);

     isPeriodicBool = false;
    latticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    latticeZeroVector = arma::drowvec(3, arma::fill::zeros);
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    init_latticeNameMap();

    latticeType = latType;

switch (latType)
{
case none:
    break;
    case cubic:
        pos = makeBCGrid(static_cast<int>(parameters[0]), static_cast<int>(parameters[1]), static_cast<int>(parameters[2]), parameters[3]);
        latticeName = "CubicGrid";
    break;
    case bulkCenteredCubic:
        pos = makeBCCGrid(static_cast<int>(parameters[0]),parameters[1], parameters[2], 1.);
        latticeName = "BCCGrid";
    break;
    case faceCenteredCubic:
        pos = makeFCCGrid(static_cast<int>(parameters[0]),parameters[1], parameters[2], 1.);
        latticeName = "FCCGrid";
    break;
    case graphene:
        pos = makeHoneyCombs(static_cast<int>(parameters[0]),parameters[1]);
        latticeName = "Graphene";
    break;

    default:
    break;
}


    setPositions(armaMatfromVectorsVector(pos), "C");

    computeMatrices();



}


FHLattice::FHLattice(FHLattice* lat)
{
    uint nSites = lat->getnSites();
    uint nDim = lat->getnDimension();
    init_latticeNameMap();
    sitePositions = arma::dmat( lat->getSitePositions()->memptr(), nSites, nDim);

    siteTypes = lat->getSiteTypes();
    siteTypeNames = lat->getSiteTypeNames();
    siteNames = lat->getSiteNames();

    isPeriodicBool = lat->isPeriodic();
    latticeMatrix = lat->getLatticeMatrix();
    latticeZeroVector = lat->getLatticeZeroVector();

    if (arma::det(latticeMatrix) != 0)
    {
        inverseLatticeMatrix = arma::inv(latticeMatrix);
    }

    latticeType = lat->getLatticeType();

    neighbourDistance = lat->getNeighbourDistance();
    samePointDistance = lat->getSamePointDistance();

    computeInverseLatticeMatrix();
    computeMatrices();

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

const FHLatticeType                       FHLattice::getLatticeType()
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

vector<std::string>            FHLattice::getSiteNames()
{
    return siteNames;
}

arma::uvec                          FHLattice::getSiteTypes()
{
    return siteTypes;
}

uint                                FHLattice::getSiteType(uint iSite)
{
    return siteTypes[iSite];
}

vector<std::string>                         FHLattice::getSiteTypeNames()       {    return siteTypeNames;}
bool                                        FHLattice::isPeriodic()             {    return isPeriodicBool;}
arma::dmat                                  FHLattice::getLatticeMatrix()       {    return latticeMatrix;}
arma::drowvec                               FHLattice::getLatticeZeroVector()   {    return latticeZeroVector;}
double                                      FHLattice::getNeighbourDistance()
{
    return neighbourDistance;
}
double                                   FHLattice::getSamePointDistance()
{
    return samePointDistance;
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


double                                      FHLattice::evaluateNearestNeighbor()
{
    double minNeib = 0;
    if (distanceMatrix.n_cols > 1)
    {
        minNeib = 1000;
        for (uint iSite = 0; iSite < distanceMatrix.n_rows; iSite++)
        {
            for (uint jSite = iSite+1; jSite < distanceMatrix.n_cols; jSite++)
            {
                if (distanceMatrix(iSite, jSite) >  samePointDistance)
                {
                    minNeib = std::min(minNeib, distanceMatrix(iSite, jSite));
                }
            }
        }
    }
    return minNeib;
}


double                                      FHLattice::getDistance(uint i, uint j)      {   return distanceMatrix(i,j);}
bool                                        FHLattice::areNeighbour(uint i, uint j)     {   return neighbourMatrix[i][j];}
bool                                        FHLattice::areSameSite(uint i, uint j)      {   return sameSiteMatrix[i][j];}
std::string                                 FHLattice::getLatticeName()                 {   return latticeName;}

const arma::dmat*                           FHLattice::getDistanceMatrix()              {   return &distanceMatrix;}
const vector<vector<bool>>*                 FHLattice::getNeighbourMatrix()             {   return &neighbourMatrix;}


// SETTERS

void                    FHLattice::setPositions(const arma::dmat& positions, std::string sitesName)
{
    sitePositions = positions;
    siteNames = vector<std::string>(sitePositions.n_rows, sitesName);
    computeMatrices();
}

void                    FHLattice::setPositions(const arma::dmat& positions, vector<std::string> sN)
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

void                    FHLattice::setSiteTypeNames(vector<std::string> typeNames)
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

void                                FHLattice::setLatticeZero(arma::drowvec latticeZ){    latticeZeroVector = latticeZ;}
void                                FHLattice::setNeighbourDistance(double neigbD){    neighbourDistance = neigbD;}
void                                FHLattice::setSamePointDistance(double samePointD){    samePointDistance = samePointD;}
void                                FHLattice::setLatticeName(std::string newName){ latticeName = newName;}

void                                FHLattice::computeInverseLatticeMatrix()
{
    inverseLatticeMatrix = arma::dmat(3,3,arma::fill::zeros);
    double rcond = 0;
    bool success = arma::inv(inverseLatticeMatrix, rcond, latticeMatrix);
    if (!success)
    {
        std::cout << "! Problem ! Lattice matrice not inversible." << std::endl;
    }

}

void                                FHLattice::readPositions(std::string xyzFile)
{

    uint nSites = 0;
    vector<vector<double>> positions3;

    siteNames = vector<std::string>();
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
            siteNames.push_back(coord);
            // Get its coordinates :
            vector<double> coordinates;
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
        std::cout << "Problem with the site position file " << xyzFile << std::endl;
    }


    if (nSites == positions3.size())
    {
        if (talkative)
        {
            std::cout << "   -- Position file read, sites :   " << nSites << std::endl;
        }
        setPositions(armaMatfromVectorsVector(positions3));
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
    latticeName = xyzFile;
    std::string fileAdress =  xyzFile + ".xyz";
   // std::cout << "Writing positions in: "<< fileAdress << std::endl;
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

    Json::Value         FHLattice::writeFHLatticeJson()
    {
        Json::Value rootLattice;

        rootLattice["latticeType"] = latticeTypeNames[latticeType];


        rootLattice["positionsFile"] = latticeName + ".xyz";

        // Periodicity
        rootLattice["periodic"] = isPeriodicBool;
        Json::Value vecJson(Json::arrayValue);
        Json::Value matJson(Json::arrayValue);
        for (uint i = 0; i<3; i++)  vecJson[i] = latticeZeroVector(i);
        for (uint i = 0; i<3; i++)
        {
            Json::Value vJson(Json::arrayValue);
            for (uint j = 0; j<3; j++)
            {
                vJson.append(latticeMatrix(i,j));
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

    unsigned int nSites = sitePositions.n_rows;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = vector<vector<bool>>(nSites,  vector<bool>(nSites, false));
    sameSiteMatrix = vector<vector<bool>>(nSites,  vector<bool>(nSites, false));


    // distance Matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            distanceMatrix(i, j) = computeDistance(i,j);
        }
    }
    // same site matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            if (std::norm(distanceMatrix(i, j))< samePointDistance)
            {
                sameSiteMatrix[i][j] = true;
            }
        }
    }
    // neighbour matrix :
    for (uint i=0; i<nSites; i++)
    {
        for (uint j=0; j<nSites; j++)
        {
            //std::cout << "i" << i << ", j=" << j << " , d=" << distanceMatrix(i, j) << "|d-n|=" << std::norm(distanceMatrix(i, j) - neighbourDistance) <<  std::endl;
            if (std::norm(distanceMatrix(i, j) - neighbourDistance) < samePointDistance)
            {
                neighbourMatrix[i][j] = true;
            }
        }
    }




    //std::cout << distanceMatrix << std::endl;
}



void FHLattice::init_latticeNameMap()
{
    latticeTypeNames = std::map<FHLatticeType, std::string>();
    latticeTypeNames[FHLatticeType::none] = "none";
    latticeTypeNames[FHLatticeType::cubic] = "grid";
    latticeTypeNames[FHLatticeType::bulkCenteredCubic] = "bcc";
    latticeTypeNames[FHLatticeType::faceCenteredCubic] = "fcc";
    latticeTypeNames[FHLatticeType::triangular] = "triangular";
    latticeTypeNames[FHLatticeType::graphene] = "graphene";
    latticeTypeNames[FHLatticeType::moire] = "moire";
    latticeTypeNames[FHLatticeType::kagome] = "kagome";

    latticeTypeNamesReverse = std::map<std::string, FHLatticeType>();

    for (std::map<FHLatticeType, std::string>::iterator it = latticeTypeNames.begin(); it != latticeTypeNames.end(); it ++)
    {
        latticeTypeNamesReverse.emplace(it->second, it->first);
    }

 }
 


FHLattice::~FHLattice()
{

}









// OTHER FUNCTIONS

// I - Grids

/// @brief Makes a simple 3D-grid of sites 
/// @param nx number of sites in the X dimension
/// @param ny number of sites in the Y dimension
/// @param nz number of sites in the Z dimension
/// @param a distance between nearest neighbours
/// @return a vector of site positions
vector<vector<double>> makeBCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a)
{
    
    unsigned int nPoints = nx*ny*nz;
    //std::cout <<"makegrid" << nx << " " << ny << " " << nz << " " << a << std::endl; 
    vector<vector<double>> pos = vector<vector<double>>(nPoints);

    unsigned int counter = 0;
    for (unsigned int z = 0; z < nz; z++)
    {
        for (unsigned int y = 0; y < ny; y++)
        {
            for (unsigned int x = 0; x < nx; x++)
            {
                //std::cout <<"makegrid in" << x << " " << y << " " << z << " " << std::endl; 
                vector<double> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
            }
        }
    }

    return pos;
}

arma::dmat findBCGridBox(unsigned int nx, unsigned int ny, unsigned int nz, double a)
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
vector<vector<double>> makeBCCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a)
{
    
    unsigned int nPoints = 2*(nx+1)*(ny+1)*(nz+1) + 2 *((nx+1)*(ny+1) + (nx+1)*(nz+1) + (ny+1)*(nz+1)) +1;
    vector<vector<double>> pos = vector<vector<double>>(nPoints);


    unsigned int counter = 0;
    for (unsigned int z = 0; z <= nz; z++)
    {
        for (unsigned int y = 0; y <= ny; y++)
        {
            for (unsigned int x = 0; x <= nx; x++)
            {
                vector<double> newPos = {a*x, a*y, a*z};
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
vector<vector<double>> makeFCCGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a)
{
    
    unsigned int nPoints = 4*(nx+1)*(ny+1)*(nz+1) + 2 *((nx+1)*(ny+1) + (nx+1)*(nz+1) + (ny+1)*(nz+1)) +1;
    vector<vector<double>> pos = vector<vector<double>>(nPoints);

    unsigned int counter = 0;
    for (unsigned int z = 0; z <= nz; z++)
    {
        for (unsigned int y = 0; y <= ny; y++)
        {
            for (unsigned int x = 0; x <= nx; x++)
            {
                vector<double> newPos = {a*x, a*y, a*z};
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
vector<vector<double>> makeHoneyCombs(unsigned int nCombs, double a)
{

    double sqrt3 = std::sqrt(3);
    unsigned int nOfPoints = 0;
    vector<vector<double>> pos = vector<vector<double>>();

    vector<vector<double>> pos0 = { {-sqrt3/2 *a, -0.5 *a, 0},
                                {-sqrt3/2 *a, 0.5 *a, 0},
                                {0   *a,      1   *a, 0},
                                {0   *a,      -1  *a, 0},
                                {sqrt3/2 *a, -0.5 *a, 0},
                                {sqrt3/2 *a, 0.5  *a, 0}};
    vector<vector<double>> posNew = pos0;


    double iVec[2] = {sqrt3 *a, 0 *a};
    double jVec[2] = {sqrt3/2 *a, 1.5 *a};

    unsigned int basisIndex[2] = {0,0};


    for (unsigned int bv = 0; bv < nCombs; bv++)
    {
        for (unsigned int p = 0; p<6; p++)
        {
            posNew[p][0] = pos0[p][0] + iVec[0] * basisIndex[0] + jVec[0] * basisIndex[1];
            posNew[p][1] = pos0[p][1] + iVec[1] * basisIndex[0] + jVec[1] * basisIndex[1];
        }
        pos.insert(pos.end(), posNew.begin(), posNew.end());

        getNextIndex(basisIndex[0], basisIndex[1]);

        nOfPoints +=6;
    }

    // We remove the doublons :
    double samePointThreshold = 0.001 * a;


    for (unsigned int i1 = 0; i1< pos.size(); i1++)         //<- we look for all points equal to positions(i1)
    {
        unsigned int nDoublonsFound = 0;                    //<- There can't be more than three times the same point
         
        unsigned int i2 = i1+1;
        while ((i2<pos.size()) && (nDoublonsFound < 2))    //<- as long as there are points and more doublons to find
        {
            double d = 0;
            for (unsigned int iDim = 0; iDim<pos[i1].size(); iDim++)
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

arma::dmat findHoneyCombGridBox(unsigned int nCombs, double a)
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


/// @brief Convert a vector of vector of double into the corresponding arma::dmat. The vector'vector and the matrix are seen in row-major order (reading/writing order).
/// @param vecVec A vector of vector of double 
/// @return The arma::dmat version of the vector
arma::dmat armaMatfromVectorsVector(vector<vector<double>>& vecVec)
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
template <typename T>       bool isIn(T w, vector<T>& v)
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
template <typename T>       bool isIn(T w, vector<T>& v, uint& firstIndex)
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

/// @brief The suite of integer doubles S defined by  S(t+1) = getNextIndex(S(t))
/// and starting by (0,0) paves the positive plan in an ordely fashion. The
/// aim is for [S(0) .. S(t)] to stay "as square as possible" for every index t.  
void getNextIndex(unsigned int& i, unsigned int& j)
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

