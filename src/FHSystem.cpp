#include "../include/FHSystem.h"


FHSystem::FHSystem()
{

    arma::dmat pos = arma::dmat(1,  1, arma::fill::zeros);

    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>();
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>();

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;

    setPositions(pos);
    setInitialConditionsDiagonal(initCond);
    computeMatrices();

    initiateEmptyEnergyMatrices();

    UInteractionType        = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType           = FHTunnelingType::nearestNeighboursConstant;
    VInteractionType        = FHSameSpinInteractionType::Vnone;
    tJExCType               = tJExchangeCouplingType::tJJnone;
    externalPotentialType   = ExternalPotentialType::ExtPotNone;
    
}


FHSystem::FHSystem(int x)
{
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
    setInitialConditionsDiagonal(initCond);
    computeMatrices();

    initiateEmptyEnergyMatrices();

    UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType = FHTunnelingType::nearestNeighboursConstant;
    VInteractionType        = FHSameSpinInteractionType::Vnone;
    tJExCType               = tJExchangeCouplingType::tJJnone;
    externalPotentialType   = ExternalPotentialType::ExtPotNone;

}


FHSystem::FHSystem(int x, int y, int z)
{

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

    int nPart = std::floor( x * y * z/2);
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
       externalPotential = FHExternalPotential(&pos);
    setInitialConditionsDiagonal(initCond);
    computeMatrices();


    initiateEmptyEnergyMatrices();
    UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType = FHTunnelingType::nearestNeighboursConstant;
        VInteractionType        = FHSameSpinInteractionType::Vnone;
    tJExCType               = tJExchangeCouplingType::tJJnone;
    externalPotentialType   = ExternalPotentialType::ExtPotNone;

}


FHSystem::FHSystem(std::string fileName, std::string folderName)
{

    enum troolean
    {
        yes,
        no,
        dontknow
    };

    enum WhatToRead
    {
        nothing,
        positions,
        symmetries,
        energyParameter_U_constant,
        energyParameter_U_bySite,
        energyParameter_U_bySiteType,
        energyParameter_U_GivenMatrix,
        energyParameter_tH_NNConstant,
        energyParameter_tH_GivenMatrix,
        energyParameter_J,
        diagonalIC_UpSpin,
        diagonalIC_DownSpin,
        numberOfParticles_UpSpin,
        numberOfParticles_DownSpin
    };

    systemName = fileName;
    systemFolder = folderName;
    std::string fileAdressName = fileName; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + fileName;
    }

    UInteractionType = FHOppositeSpinInteractionType::Unone;
    VInteractionType = FHSameSpinInteractionType::Vnone;
    tunnelingType = FHTunnelingType::tHnone;
    tJExCType = tJExchangeCouplingType::tJJnone;
    externalPotentialType   = ExternalPotentialType::ExtPotNone;

    WhatToRead whatToR =  WhatToRead::nothing;  
    int nDim = -1;
    troolean orbitalHaveNames = troolean::dontknow;
    bool positionsReady = false;

    std::ifstream systFile (fileAdressName);
    if (systFile.is_open())
    {
        printf("== Reading %s \n", fileName.c_str());

        std::vector<unsigned int>           orbitalTypes    = std::vector<unsigned int>();
        std::vector<std::string>            orbitalNames    = std::vector<std::string>();
        std::vector<std::vector<double>>    posVecVec       = std::vector<std::vector<double>>();
        std::vector<std::vector<double>>    symVecVec       = std::vector<std::vector<double>>();
        std::vector<unsigned int>           initCond_up     = std::vector<unsigned int>();
        std::vector<unsigned int>           initCond_down   = std::vector<unsigned int>();
        std::vector<unsigned int>           particlesNumber(2, 0);

        std::vector<double> Jvector                         = std::vector<double>();
        std::vector<std::vector<double>> tHMatrix           = std::vector<std::vector<double>>();
        std::vector<std::vector<double>> UMatrix            = std::vector<std::vector<double>>();


        std::string line;
        while (std::getline(systFile, line))
        {
           if (line == "#POSITIONS")                                                //<- ready to read the positions
           {
                whatToR = WhatToRead::positions;
           }
           else if (line == "#POSITIONS_END")                                       //<- the positions stored in "pos" 
           {
                arma::dmat sitePos = armaMatfromVectorsVector(posVecVec);
                setPositions(sitePos);

                if (orbitalHaveNames == troolean::yes)
                {
                    std::vector<std::string>            orbitalNameTypes    = std::vector<std::string>();
                    uint numberOfNames = 0;
                    uint index = 0;
                    for (uint nameInt = 0; nameInt < orbitalNames.size(); nameInt++)
                    {
                        if (isIn( orbitalNames[nameInt], orbitalNameTypes, index))
                        {
                            orbitalTypes[nameInt] = index;
                        }
                        else
                        {
                            orbitalNameTypes.push_back(orbitalNames[nameInt]);
                            orbitalTypes[nameInt] = orbitalNameTypes.size() - 1;
                        }
                    }
                }

                initiateEmptyEnergyMatrices();
                positionsReady = true;
                whatToR = WhatToRead::nothing;
           }
           else if (line == "#SYMMETRIES")
           {
                whatToR = WhatToRead::symmetries;
           }
           else if (line == "#SYMMETRIES_END")
           {
                arma::dmat symmetryArmaMatrix = armaMatfromVectorsVector(symVecVec);
                setSymmetry(symmetryArmaMatrix);
                whatToR = WhatToRead::nothing;
           }
           else if (line == "#INITIALCONDITIONS_GIVENDIAGONAL")
           {
                initialCondition = FHInitialConditionType::givenDiagonal;
                whatToR = WhatToRead::diagonalIC_UpSpin;
           }
           else if (line == "#INITIALCONDITIONS_QUENCHED")
           {
                initialCondition = FHInitialConditionType::quenched;
                whatToR = WhatToRead::numberOfParticles_UpSpin;
           }
           else if (line == "#ENERGYPARAMETER_U_CONSTANT")
           {
                whatToR = WhatToRead::energyParameter_U_constant;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
           }
           else if (line == "#ENERGYPARAMETER_U_BYSITE")
           {
                whatToR = WhatToRead::energyParameter_U_bySite;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteGivenDiagonal;
           }
           else if (line == "#ENERGYPARAMETER_U_BYSITETYPE")
           {
                whatToR = WhatToRead::energyParameter_U_bySiteType;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteBySiteType;
           }
           else if (line == "#ENERGYPARAMETER_tH_NEARESTNEIGHBORS")
           {
                whatToR = WhatToRead::energyParameter_tH_NNConstant;
                tunnelingType = FHTunnelingType::nearestNeighboursConstant; 
           }
           else if (line == "#ENERGYPARAMETER_tH_GIVENMATRIX")
           {
                whatToR = WhatToRead::energyParameter_tH_GivenMatrix;
                tunnelingType = FHTunnelingType::tunnelingGivenMatrix;
           }
            else if (line == "#ENERGYPARAMETER_J")
           {
                whatToR = WhatToRead::energyParameter_J;
           }
           else switch (whatToR)
           {
            case WhatToRead::positions :  // Read positions
            {

                if (orbitalHaveNames == troolean::dontknow)
                {
                    std::stringstream st(line);
                    std::string coordt;
                    orbitalHaveNames = troolean::no;
                    std::getline (st, coordt, ' ');
                    try {
                        std::stod(coordt);
                    }catch(const std::invalid_argument& e)  {
                        std::cout << "Orbitals have names !" << std::endl;
                        orbitalHaveNames = troolean::yes;
                    }
                }

                std::stringstream ss(line);
                std::string coord;
                std::vector<double> coordinates = std::vector<double>();

                if (orbitalHaveNames == troolean::yes)
                {
                    std::getline (ss, coord, ' ');
                    orbitalNames.push_back(coord);
                }

                while (std::getline (ss, coord, ' ')) 
                {
                    double d = std::stod(coord);
                    coordinates.push_back(d);
                }

                if (nDim == -1)
                {
                    nDim = coordinates.size();
                    posVecVec.push_back(coordinates);
                }
                else
                {
                    if (coordinates.size() == nDim)
                    {
                        posVecVec.push_back(coordinates);
                    }
                    else
                    {
                        printf("Dimension problem with line \"%s\"", line.c_str());
                        return;
                    }
                }


                break;
            }
            case WhatToRead::symmetries :
            {
                std::stringstream ss(line);
                std::string coord;
                std::vector<double> coordinates = std::vector<double>();
                while (std::getline (ss, coord, ' ')) 
                {
                    double d = std::stod(coord);
                    coordinates.push_back(d);
                }
                if (coordinates.size() == nDim)
                {
                    symVecVec.push_back(coordinates);
                }
                else
                {
                    printf("Dimension problem with line \"%s\"", line.c_str());
                    return;
                }
            }
            case WhatToRead::energyParameter_U_constant :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::getline (ss, U_string, ' ');
                double U0 = std::stod(U_string);
                
                setOppositeSpinInteractionType(FHOppositeSpinInteractionType::UonSiteConstant);
                setUInteraction_OnSite(U0);
                break;
            }
            case WhatToRead::energyParameter_U_bySite :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::vector<double> Uvector = std::vector<double>();
                while (std::getline (ss, U_string, ' ')) 
                {
                     Uvector.push_back(std::stod(U_string));
                }     
                arma::dvec Udiag = arma::dvec(Uvector);
                setOppositeSpinInteractionType(FHOppositeSpinInteractionType::UonSiteGivenDiagonal);
                setUInteraction_GivenDiagonal(Udiag);
                break;
            }
            case WhatToRead::energyParameter_U_bySiteType :
            {

                std::stringstream ss(line);
                std::string site_name;
                std::string U_value;
                std::getline (ss, site_name, ' ');
                std::getline (ss, U_value, ' ');
                std::vector<std::string>::iterator itName = std::find(orbitalNames.begin(), orbitalNames.end(), site_name);
                std::vector<double> Uvector = std::vector<double>(orbitalNames.size());
                Uvector[itName - orbitalNames.begin()] = std::stod(U_value);          
                break;
            }
            case WhatToRead::energyParameter_U_GivenMatrix :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::vector<double> Uvector = std::vector<double>();
                while (std::getline (ss, U_string, ' ')) 
                {
                     Uvector.push_back(std::stod(U_string));
                }
                UMatrix.push_back(Uvector);
                if (UMatrix.size() == posVecVec.size())
                {
                    arma::dmat UM = armaMatfromVectorsVector(UMatrix);
                    setUInteraction_GivenMatrix(UM);
                }    
                break;
            }
            case WhatToRead::energyParameter_tH_NNConstant :
            {
                std::stringstream ss(line);
                std::string t_string;
                std::getline (ss, t_string, ' ');
                double tH1 = std::stod(t_string);
                setTunneling_BetweenNeighbours(tH1);
                break;
            }
            case WhatToRead::energyParameter_tH_GivenMatrix :
            {
                std::stringstream ss(line);
                std::string t_string;
                std::vector<double> tHvector = std::vector<double>();
                while (std::getline (ss, t_string, ' ')) 
                {
                     tHvector.push_back(std::stod(t_string));
                }
                tHMatrix.push_back(tHvector);
                if (tHMatrix.size() == posVecVec.size())
                {
                    arma::dmat tHM = armaMatfromVectorsVector(tHMatrix);
                    setTunneling_GivenMatrix(tHM);
                }    
                break;
            }
            case WhatToRead::energyParameter_J :
            {
                break;
            }
           case WhatToRead::diagonalIC_UpSpin :  // Read  up spin particles
           {
                std::stringstream ss(line);
                std::string siteN;
                
                while (std::getline(ss, siteN, ' '))            //<- We read the line word by word
                {
                    //std::cout << "Reading initial positions of up-spin particles, siteN = " << siteN << std::endl;
                    initCond_up.push_back(std::stoi(siteN));    //<- And fill the inital condition vector with the site numbers
                    particlesNumber[0]++;
                }
                whatToR = WhatToRead::diagonalIC_DownSpin;      //<- At the end, the up-spin particles positions have been read.
                break;
            }
            case WhatToRead::diagonalIC_DownSpin :              //<- Else, if the down-spin particles positions have not been read yet
            {
                std::stringstream ss(line);
                std::string siteN;
                while (std::getline(ss, siteN, ' ')) 
                {
                    // std::cout << "Reading initial positions of down-spin particles, siteN = " << siteN << std::endl;
                    initCond_down.push_back(std::stoi(siteN)); 
                    particlesNumber[1]++;
                }
                whatToR = WhatToRead::nothing; 
                break;
            }
            case WhatToRead::numberOfParticles_UpSpin :
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[0] = std::stoi(nPart);
                whatToR = WhatToRead::numberOfParticles_DownSpin;
            }
            break;
            case WhatToRead::numberOfParticles_DownSpin :
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[1] = std::stoi(nPart);
                whatToR = WhatToRead::nothing;
            }
            break;


           }
        }

        systFile.close();



    if (initialCondition == FHInitialConditionType::givenDiagonal)
    {
            std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
            initCond[0] = initCond_up;
            initCond[1] = initCond_down;

            setInitialConditionsDiagonal(initCond);
    }
    else if (initialCondition == FHInitialConditionType::quenched)
    {
        setInitialConditions(initialCondition, particlesNumber);
    }


    computeMatrices();

    }
    else
    {
        printf("Can't open %s", fileName.c_str());
    }


    
}




FHSystem::FHSystem(Json::Value rootSystem, std::string folderName)
{


// Positions
    std::string positionFile = rootSystem.get("positions", "").asString();

       

    if (folderName !="")
    {
        positionFile = folderName + "/" + positionFile;
    }

    std::cout << positionFile << std::endl;
    readPositions(positionFile);

    std::cout << sitePositions << std::endl;

    initiateEmptyEnergyMatrices();


// Initial condition
    std::string initCond = rootSystem.get("IC", "spinWave").asString();
    if (initCond == "spinWave")
    {
        initialCondition = FHInitialConditionType::spinWave;
        for (uint iPart = 0; iPart < sitePositions.n_rows ; iPart++ )
        {
            if (iPart%2 == 0)
            {
                initConditions_upSpin.push_back(true);
                initConditions_downSpin.push_back(false);
            }
            else
            {
                initConditions_upSpin.push_back(false);
                initConditions_downSpin.push_back(true);
            }
        }
    }
    else if (initCond == "chargeWave")
    {
        initialCondition = FHInitialConditionType::chargeWave;
        for (uint iPart = 0; iPart < sitePositions.n_rows ; iPart++ )
        {
            if (iPart%2 == 0)
            {
                initConditions_upSpin.push_back(true);
                initConditions_downSpin.push_back(true);
            }
            else
            {
                initConditions_upSpin.push_back(false);
                initConditions_downSpin.push_back(false);
            }
        }
    }
    else if (initCond == "quenched")
    {

    }
    else if (initCond == "givenInFile")
    {

    }

// Tunneling Energy
    std::string tunnEn = rootSystem["EnergyT"].asString();
    if (tunnEn == "nearestNeighbours")
    {
        std::cout << "Hello tunn "<< std::endl;
        tunnelingType = FHTunnelingType::nearestNeighboursConstant;
        double tH = rootSystem.get("parametersT", -7.).asDouble();
        std::cout << "Hello tunn tH =" << tH << std::endl;
        setTunneling_BetweenNeighbours(tH);
    }

// Opposite-spin interaction
    std::string UEn = rootSystem["EnergyU"].asString();
    if (UEn == "None")
    {
        UInteractionType = FHOppositeSpinInteractionType::Unone;
    }
    if (UEn == "onSite")
    {
       setOppositeSpinInteractionType(FHOppositeSpinInteractionType::UonSiteConstant);
        double U0 = rootSystem.get("parametersU", 1.).asDouble();                
        setUInteraction_OnSite(U0);
    }

    // Same-spin Interaction
    std::string VEn = rootSystem["EnergyV"].asString();
    if (VEn == "None")
    {
        VInteractionType = FHSameSpinInteractionType::Vnone;
    }

}



/// @brief getter - number of sites
/// @return number of sites
unsigned int                                FHSystem::getnSites()
{
    return sitePositions.n_rows;
}

unsigned int                                FHSystem::getnDimension()
{
    return sitePositions.n_cols;
}

arma::drowvec                               FHSystem::getSitePosition(uint site)
{
    return sitePositions.row(site);
}

/// @brief getter - coordinate of a site
/// @param dim coordinate dimension
/// @param pos site number
/// @return coordinate of site 'pos'in dimension 'dim'
double                                      FHSystem::getCoordinate(uint dim, uint iSite)
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


bool                                        FHSystem::sitesHaveNames()
{
        return (siteTypes.n_rows == sitePositions.n_rows);
}

std::string                                 FHSystem::getSiteName(uint site)
{
    return siteTypeNames[siteTypes[site]];
}


bool                                        FHSystem::isSymmetric()
{
    return (symmetryVectors.n_rows > 0);
}

arma::dmat                                  FHSystem::getSymmetryMatrix()
{
    return symmetryVectors;
}

bool                                        FHSystem::isSystemTimeDependant()
{
    return externalPotential.isThePotentialTimeDependant();
}
bool                                    FHSystem::isExternalPotentialNull()
{
    return externalPotential.isThePotentialNull();
}

double                                      FHSystem::computeDistance(uint iSite, uint jSite)
{
    arma::drowvec dVec = sitePositions.row(iSite) - sitePositions.row(jSite); 
    double d = arma::norm(dVec);
    if (isSymmetric())
    {
        uint symmetryNumber3 = pow(symmetryVectors.n_rows,3);                                                   //<- Number of vector to check : nSym^3

        for (uint lCombination = 0; lCombination < symmetryNumber3; lCombination++)                             //<- for every vector to check
        {
            arma::drowvec dVecSym = dVec;
            for (uint kDim = 0; kDim < sitePositions.n_cols ; kDim ++)                                          //<- vecSym = vec + sum_k [-1,0,1] * sym(k) 
            {
                dVecSym += ((((uint)std::floor(lCombination/pow(3,kDim))) % 3) - 1) * symmetryVectors.row(kDim);
            } 
            d = std::min(d, arma::norm(dVecSym));                                                               //<- We replace if its norm is smaller
        } 
    }

    return d;
}

double                                      FHSystem::getDistance(uint i, uint j)
{
    return distanceMatrix(i,j);
}

bool                                        FHSystem::areNeighbour(uint i, uint j)
{
    return neighbourMatrix[i][j];
}


/// @brief  getter - initial condition type
/// @return initial condition type
FHSystem::FHInitialConditionType            FHSystem::getInitialCondition()
{
    return initialCondition;
}

std::vector<bool>                           FHSystem::getInitialCondition_upSpin()
{
    return initConditions_upSpin;
}

std::vector<bool>                           FHSystem::getInitialCondition_downSpin()
{
    return initConditions_downSpin;
}

/// @brief getter - number of spin-up particles
/// @return number of spin-up particles
unsigned int                                FHSystem::getUpSpinParticleNumber()
{
    return upSpinParticleNumber;
}

/// @brief getter - number of spin-down particles
/// @return number of spin-down particle
unsigned int                                FHSystem::getDownSpinParticleNumber()
{
    return downSpinParticleNumber;
}



FHSystem::FHOppositeSpinInteractionType     FHSystem::getOppositeSpinInteractionType()
{
    return UInteractionType;
}

FHSystem::FHTunnelingType                   FHSystem::getTunnelingType() 
{
    return tunnelingType;
}

FHSystem::FHSameSpinInteractionType         FHSystem::getSameSpinInteractionType() 
{
    return VInteractionType;
}

FHSystem::tJExchangeCouplingType            FHSystem::getExchangeCouplingEnergyType() 
{
    return tJExCType;
}


std::complex<double>                        FHSystem::getInteractionEnergy(int i, int j)
{
    return UMatrix(i,j);
}

std::complex<double>                        FHSystem::getSameSpinInteractionEnergy(int i, int j)
{
    return VMatrix(i,j);
}

std::complex<double>                        FHSystem::getTunnelingEnergy(int i, int j)
{
    return tMatrix(i,j);
}

std::complex<double>                        FHSystem::getExchangeCouplingEnergy(int i, int j)
{
    return JMatrix(i,j);
}

std::complex<double>                        FHSystem::getExternalPotential(int i, double t)
{
    return externalPotential.getExternalPotential(t,i,i);
}

std::complex<double>                        FHSystem::getExternalPotential(int i, int j, double t)
{
    return externalPotential.getExternalPotential(t,i,j);
}

arma::cx_dmat                        FHSystem::getExternalPotential(double time)
{
    return externalPotential.getExternalPotential(time);
}

void                                FHSystem::updateExternalPotential(double time)
{
    externalPotential.updatePotential(time);
}


const arma::cx_dmat*                        FHSystem::getUMatrixPointer()
{
    return &UMatrix;
}

const arma::cx_dmat*                        FHSystem::getVMatrixPointer()
{
    return &VMatrix; 
}

const arma::cx_dmat*                        FHSystem::getJMatrixPointer()
{
    return &JMatrix;
}

const arma::cx_dmat*                        FHSystem::gettMatrixPointer()
{
    return &tMatrix;
}

const arma::cx_dmat*                        FHSystem::getExtPotPointer()
{
    return externalPotential.getExternalPotentialPointer();
}



std::string                                 FHSystem::getName()
{
    return systemName;
}

std::string                                 FHSystem::getFolder()
{
    return systemFolder;
}



// SETTERS

void                    FHSystem::setPositions(const arma::dmat& positions, std::string sitesName)
{
    sitePositions = positions;
    siteNames = std::vector<std::string>(sitePositions.n_rows, sitesName);
    computeMatrices();
}

void                    FHSystem::setPositions(const arma::dmat& positions, std::vector<std::string> sN)
{
    sitePositions = positions;
    siteNames = sN;
    computeMatrices();
}


bool                    FHSystem::setSiteTypes(arma::uvec sTypes)
{
    siteTypes = sTypes;
    return sitesHaveNames();
}

void                    FHSystem::setSiteTypeNames(std::vector<std::string> typeNames)
{
    siteTypeNames = typeNames;
}


void                    FHSystem::setSymmetry(arma::dmat symVectors)
{
    if (symVectors.n_cols != sitePositions.n_cols)
    {
        std::cout << "Symmetry vectors have a wrong dimension number." << std::endl;
    }
    else if (symVectors.n_rows > sitePositions.n_cols)
    {
        std::cout << "Symmetry vectors are too numerous for the system's dimension." << std::endl;
    }
    else
    {
        symmetryVectors = symVectors;
    }

}


/// @brief Set up the Fermi-Hubbard system from site positions and vectors of site number as initial condition
/// @param diagonalInitConditions vectors of sites on which there is a up-spin/down-spin particle
/// @return void
void FHSystem::setInitialConditionsDiagonal(const std::vector<std::vector<unsigned int>>& diagonalInitConditions )
{

    initialCondition =  FHSystem::FHInitialConditionType::givenDiagonal;

    upSpinParticleNumber = diagonalInitConditions[0].size();
    downSpinParticleNumber = diagonalInitConditions[1].size();

    initConditions_upSpin = std::vector<bool>(sitePositions.n_rows, false);
    initConditions_downSpin = std::vector<bool>(sitePositions.n_rows, false);



    for (unsigned int itSite = 0; itSite < diagonalInitConditions[0].size(); itSite++)
    {
        initConditions_upSpin[diagonalInitConditions[0][itSite]] = true;
    }
    for (unsigned int itSite = 0; itSite < diagonalInitConditions[1].size(); itSite++)
    {
        initConditions_downSpin[diagonalInitConditions[1][itSite]] = true;
    }

}


/// @brief Set up the Fermi-Hubbard system from site positions the number of particles as initial condition
/// @param initCondType a kind of initial condition (quenched, charge-wave, spin-wave ...)
/// @return void
void FHSystem::setInitialConditions(FHInitialConditionType initCondType)
{

    initialCondition = initCondType;

    initConditions_upSpin = std::vector<bool>(sitePositions.n_rows, false);
    initConditions_downSpin = std::vector<bool>(sitePositions.n_rows, false);
    upSpinParticleNumber = 0;
    downSpinParticleNumber = 0 ; 

    if (initialCondition == FHInitialConditionType::spinWave)
    {
        for (uint iSite = 0; iSite < sitePositions.n_rows; iSite ++)
        {
            if (iSite%2 == 0)
            {
                initConditions_upSpin[iSite] = true;
                upSpinParticleNumber++;
            }
            else
            {
                initConditions_downSpin[iSite] = true;
                downSpinParticleNumber++;
            }

        }
    }
    else if (initialCondition == FHInitialConditionType::spinWave)
    {
        for (uint iSite = 0; iSite < sitePositions.n_rows; iSite ++)
        {
            if (iSite%2 == 0)
            {
                initConditions_upSpin[iSite] = true;
                upSpinParticleNumber++;
                initConditions_downSpin[iSite] = true;
                downSpinParticleNumber++;
            }
        }
    }
    else
    {
        upSpinParticleNumber = sitePositions.n_rows/2;
        downSpinParticleNumber = sitePositions.n_rows/2;
    }

}

/// @brief Set up the Fermi-Hubbard system from site positions the number of particles as initial condition
/// @param initCondType a kind of initial condition (quenched, charge-wave, spin-wave ...)
/// @param partNumner vector of the numbers of up-spin, down-spin particle
/// @return void
void FHSystem::setInitialConditions(FHInitialConditionType initCondType, std::vector<unsigned int> partNumber )
{

    initialCondition = initCondType;
    initConditions_upSpin = std::vector<bool>(sitePositions.n_rows, false);
    initConditions_downSpin = std::vector<bool>(sitePositions.n_rows, false);

    upSpinParticleNumber = partNumber[0];
    downSpinParticleNumber = partNumber[1];
    if (initialCondition == FHInitialConditionType::spinWave)
    {
        for (uint iSite = 0; iSite < std::min((uint)sitePositions.n_rows, upSpinParticleNumber); iSite +=2)
        {
                initConditions_upSpin[iSite] = true;
        }
        for (uint iSite = 1; iSite < std::min((uint)sitePositions.n_rows, downSpinParticleNumber); iSite +=2)
        {
                initConditions_downSpin[iSite] = true;
        }

    }
    else if (initialCondition == FHInitialConditionType::spinWave)
    {
        for (uint iSite = 0; iSite < std::min((uint)sitePositions.n_rows, upSpinParticleNumber); iSite +=2)
        {
                initConditions_upSpin[iSite] = true;
        }
        for (uint iSite = 0; iSite < std::min((uint)sitePositions.n_rows, downSpinParticleNumber); iSite +=2)
        {
                initConditions_downSpin[iSite] = true;
        }
    }

}


void                    FHSystem::setOppositeSpinInteractionType(FHSystem::FHOppositeSpinInteractionType UType)
{
    UInteractionType = UType;
}

void                    FHSystem::setTunnelingType(FHSystem::FHTunnelingType tType)
{
    tunnelingType = tType;
}

void                    FHSystem::setSameSpinInteractionType(FHSystem::FHSameSpinInteractionType VType)
{
    VInteractionType = VType;
}

void                    FHSystem::setExchangeCouplingEnergyType(FHSystem::tJExchangeCouplingType JType)
{
    tJExCType = JType;
}



/// @brief setter - on-site interaction energy
/// @param U0 on-site interaction energy
void                    FHSystem::setUInteraction_OnSite(double U0, bool add)
{
    if (add)
    {
        UMatrix += U0 * arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::eye);
    }
    else
    {
        UMatrix = U0 * arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::eye);
    }

}

/// @brief setter - On-site interaction energy
/// @param U Vector of on-site interaction energy
void                    FHSystem::setUInteraction_GivenDiagonal(arma::dvec& Ua, bool add)
{
    if (Ua.n_rows == sitePositions.n_rows) 
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(arma::diagmat(Ua), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(arma::diagmat(Ua), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }      
    }
    else
    {
        std::cout << "Wrong size for the Ua interaction vector." << std::endl;
    }
}

/// @brief setter - Matrix interaction energy
/// @param Uab Matrix interaction energy
void                    FHSystem::setUInteraction_GivenMatrix(arma::dmat& Uab, bool add)
{
    if ((Uab.n_rows == sitePositions.n_rows) && (Uab.n_cols == sitePositions.n_rows))
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(Uab, arma::dmat(Uab.n_rows, Uab.n_cols, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(Uab, arma::dmat(Uab.n_rows, Uab.n_cols, arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the Uab interaction matrix." << std::endl;
    }

}


/// @brief setter - Matrix interaction energy
/// @param Us Matrix interaction energy by site type
void                    FHSystem::setUInteraction_BySiteType(arma::dvec& Us, bool add)
{
    if (Us.n_rows == siteTypeNames.size()) 
    {
        arma::dvec Uvec = arma::dvec(sitePositions.n_rows, arma::fill::zeros);
        for (uint iSite = 0; iSite < sitePositions.n_rows; iSite++)
        {
            Uvec[iSite] = Us[siteTypes[iSite]];
        }
        if (add)
        {
            UMatrix += arma::cx_dvec(arma::diagmat(Uvec), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dvec(arma::diagmat(Uvec), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the Us interaction matrix." << std::endl;
    }
}

/// @brief setter - nearest neighbout tunneling interaction energy
/// @param tH1 tunneling interaction energy
void                    FHSystem::setTunneling_BetweenNeighbours(double tH1, bool add)
{
    std::cout << "Tunn function, nrow=" << sitePositions.n_rows << " ncol=" << sitePositions.n_cols << std::endl;
    if (!add)
    {
        tMatrix = arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
    }
    for (unsigned int i = 0; i< sitePositions.n_rows; i++)
    {
        for (unsigned int j = 0; j< sitePositions.n_rows; j++)
        {
            if (neighbourMatrix[i][j])
            {
                tMatrix(i,j) += tH1;
            }
        }
    }
    //printf("J-Matrix : \n");
    //tMatrix.print();
}

/// @brief setter - Same site tunneling interaction energy
/// @param tH0 tunneling interaction energy
void                    FHSystem::setTunneling_SamePosition(double tH0, bool add)
{
    if (!add)
    {
        tMatrix = arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
    }
    for (unsigned int i = 0; i< sitePositions.n_rows; i++)
    {
        for (unsigned int j = 0; j< sitePositions.n_rows; j++)
        {
            if (sameSiteMatrix[i][j] && (i!=j))
            {
                tMatrix(i,j) += tH0;
            }
        }
    }
}

/// @brief setter - Tunneling interaction energy on sites at distance 'dist'
/// @param tHd tunneling interaction energy
/// @param dist distance between two sites
void                    FHSystem::setTunneling_atDistance(double tHd, double dist, bool add)
{
    if (!add)
    {
        tMatrix = arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
    }
    for (unsigned int i = 0; i< sitePositions.n_rows; i++)
    {
        for (unsigned int j = 0; j< sitePositions.n_rows; j++)
        {
            if (abs(distanceMatrix(i,j) - dist) < 0.0001)
            {
                tMatrix(i,j) += tHd;
            }
        }
    }
}

/// @brief setter - matrix tunneling interaction energy
/// @param tH tunneling interaction energy by pair of sites
void                    FHSystem::setTunneling_GivenMatrix(arma::dmat tH, bool add)
{
    if ((tH.n_rows == sitePositions.n_rows) && (tH.n_cols == sitePositions.n_rows))
    {
        if (add)
        {
            tMatrix += arma::cx_dmat(tH, arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
        else
        {
            tMatrix = arma::cx_dmat(tH, arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the tunneling matrix." << std::endl;
    }
}

/// @brief setter - matrix tunneling interaction energy
/// @param tHs1s2 tunneling interaction energy by pair of site types
void                    FHSystem::setTunneling_BySiteType(arma::dmat& tHs1s2, double d, bool add)
{
    if ((tHs1s2.n_rows == siteTypeNames.size()) && (tHs1s2.n_cols == siteTypeNames.size()))
    {
          arma::dvec tHMat = arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
        for (uint iSite = 0; iSite < sitePositions.n_rows; iSite++)
        {
            for (uint jSite = 0; jSite < sitePositions.n_rows; jSite++)
            {
                if (std::abs(distanceMatrix(iSite, jSite) - d) < 0.0001)
                {
                    tHMat(iSite, jSite) = tHs1s2(siteTypes[iSite], siteTypes[jSite]);
                }
            }
        }
        if (add)
        {
            tMatrix += arma::cx_dvec(tHMat, arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
        else
        {
            tMatrix = arma::cx_dvec(tHMat, arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
    }
}

void                    FHSystem::setExchangeCoupling(double J, bool add)
{
    JMatrix = arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
}


void                    FHSystem::setName(std::string name)    
{
    systemName = name;
}


void                    FHSystem::makePotential_pulse()
{
    double intensityA0 = 1;
    double centert0 = 1.;
    double widthtd = 0.1;
    double frequencyw0 = 100;
    arma::dvec polar = arma::dvec(sitePositions.n_cols, arma::fill::zeros);
    polar(0) = 1/std::sqrt(2);
    polar(1) = 1/std::sqrt(2);
    externalPotential = FHExternalPotentialPulse(&sitePositions, intensityA0, centert0, widthtd, frequencyw0, polar, &tMatrix);
}

void                    FHSystem::readPositions(std::string xyzFile)
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
        std::cout << "Position file read, " << nPart << " particles." << std::endl;
        setPositions(armaMatfromVectorsVector(positions3));
    }
    else
    {
        std::cout << "Problem in the xyz file : " << nPart << " particles announced, " << positions3.size() << " found." << std::endl;
    }

}



void                 FHSystem::writePositions(std::string xyzFile, std::string folderName)
{
        // Writing of everything on the file :
    std::string fileAdress =  xyzFile + ".xyz";
    std::cout << "Writing positions in: "<< fileAdress << std::endl;
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

void                    FHSystem::readPotential(std::string fileName, std::string folderName)
{

    std::string fileAdressName = fileName; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + fileName;
    }



    std::ifstream systFile(fileAdressName);
    if (systFile.is_open())
    {
        std::string line;
        std::string coord;
        std::getline(systFile, line);
        
        if (line == "#EXTERNALPOTENTIAL_NULL")                                                //<- ready to read the positions
        {
            externalPotentialType = ExternalPotentialType::ExtPotNone;
        }
        else if (line == "#EXTERNALPOTENTIAL_CONSTANT")   
        {
            externalPotentialType = ExternalPotentialType::ExtPotConstant;
        }
        else if (line == "#EXTERNALPOTENTIAL_PULSE")
        {
            externalPotentialType = ExternalPotentialType::ExtPotPulse;
            std::getline(systFile, line);
            double intensityA0 = std::stod(line);
            std::getline(systFile, line);
            double centert0 = std::stod(line);
            double widthtd = std::stod(line);
            std::getline(systFile, line);
            double frequencyw0 = std::stod(line);
            std::getline(systFile, line);
            arma::dvec polar = arma::dvec(sitePositions.n_cols, arma::fill::zeros);
            std::stringstream ss(line);
            uint coordCountor = 0;
            while (std::getline(ss, coord, ' ')) 
            {
                if (coordCountor < sitePositions.n_cols)
                {
                    polar(coordCountor) = std::stod(coord);
                }
                coordCountor++;
            }
            externalPotential = FHExternalPotentialPulse(&sitePositions, intensityA0, centert0, widthtd, frequencyw0, polar, &tMatrix);
        } 
    }
    else
    {
        std::cout << "External Potential " << fileName << " cannot be read" << std::endl;
    }
}


/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void                    FHSystem::computeMatrices()
{

    unsigned int nSites = sitePositions.n_rows;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    sameSiteMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
    for (int i=0; i<nSites; i++)
    {
        for (int j=0; j<nSites; j++)
        {
            distanceMatrix(i, j) = computeDistance(i,j);
            if (abs(distanceMatrix(i, j))< 0.0001)
            {
                sameSiteMatrix[i][j] = true;
            }
            if (abs(distanceMatrix(i, j) - 1)< 0.0001)
            {
                neighbourMatrix[i][j] = true;
            }
        }
    }
}


void                    FHSystem::initiateEmptyEnergyMatrices()
{
    UInteractionType = FHOppositeSpinInteractionType::Unone;
    tunnelingType = FHTunnelingType::tHnone;
    VInteractionType = FHSameSpinInteractionType::Vnone;
    tJExCType = tJExchangeCouplingType::tJJnone;

    uint nSites = sitePositions.n_rows;

    UMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    VMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    JMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    tMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    externalPotential = FHExternalPotential(&sitePositions);
}



/// @brief Write a JSON version of the Fermi-Hubbard System
/// @param folderName emplacement of the folder
Json::Value                    FHSystem::writeFHSystemJson()
{


    Json::Value systemRoot;
//Positions
    systemRoot["positions"] = systemName + ".xyz";

//Initial condition
    switch (initialCondition)
    {
    case FHInitialConditionType::spinWave:
        systemRoot["IC"] = "spinWave"; 
        break;
    case FHInitialConditionType::chargeWave:
        systemRoot["IC"] = "chargeWave"; 
        break;
    case FHInitialConditionType::quenched:
        systemRoot["IC"] = "quenched"; 
        break;
    case FHInitialConditionType::groundState:
        systemRoot["IC"] = "groundState"; 
        break;
    case FHInitialConditionType::givenDiagonal:
        systemRoot["IC"] = "givenDiagonal"; 
        break;
    default:
        break;
    } 

// Energy
    switch (tunnelingType)
    {
    case FHTunnelingType::tHnone:
        systemRoot["EnergyT"] = "None";
        break;   
    case FHTunnelingType::nearestNeighboursConstant:
    {
        systemRoot["EnergyT"] = "nearestNeighbours";
        arma::uvec q1 = find(tMatrix,1);
        if (q1.n_rows !=0)
        {
            double th = std::real( tMatrix(q1(0)) );
            systemRoot["parametersT"] = th;
        }
        else
        {
            systemRoot["parametersT"] = 0;
        }
        break;
    }
    default:
        break;
    }

    switch (UInteractionType)
    {
    case FHOppositeSpinInteractionType::Unone:
        systemRoot["EnergyU"] = "None";
        break;   
    case FHOppositeSpinInteractionType::UonSiteConstant:
    {
        systemRoot["EnergyU"] = "onSite";
        systemRoot["parametersU"] = std::real( UMatrix(0,0));
        break;
    }
    default:
        break;
    }

    switch (VInteractionType)
    {
    case FHSameSpinInteractionType::Vnone:
        systemRoot["EnergyV"] = "None";
        break;   
    case FHSameSpinInteractionType::VonSiteConstant:
    {
        systemRoot["EnergyV"] = "onSite";
        systemRoot["parametersV"] = std::real( VMatrix(0,0));
        break;
    }
    default:
        break;
    }



    return systemRoot;
}


/// @brief Write an ASCII version of the Fermi-Hubbard System
/// @param folderName emplacement of the folder
void                    FHSystem::writeFHSystemBis(std::string folderName)
{
    // Writing of everything on the file :
    std::string fileAdress =  systemName;

    std::cout << "In system, system name : "<< systemName << std::endl;
    if (folderName != "")
    {
        fileAdress = folderName + "/" + systemName;
    }
   std::ofstream streamSystem(fileAdress);

   std::cout << "In system, file adress : "<< fileAdress << std::endl;
   
    if (streamSystem.is_open())
    {
// I - POSITIONS
    // I - 1 POSITION OF SITES
        bool sitesHaveTags = sitesHaveNames();

        streamSystem << "#POSITIONS" << std::endl;
        for (uint itpos = 0; itpos < sitePositions.n_rows; itpos++)
        {
            if (sitesHaveTags)
            {
                streamSystem << siteTypeNames[siteTypes[itpos]] << "\t" ;
            }
            for(uint itdim = 0 ; itdim < sitePositions.n_cols; itdim++)
            {
                streamSystem << sitePositions(itpos, itdim) << "\t";
            }
            streamSystem << std::endl;          
        }
        streamSystem << "#POSITIONS_END" << std::endl;

    // I - 2 SYMMETRIES
        if (isSymmetric())
        {
            streamSystem << "#SYMMETRIES" << std::endl;
            for (uint itSym = 0; itSym < symmetryVectors.n_rows; itSym++)
            {
                for(uint itdim = 0 ; itdim < symmetryVectors.n_cols; itdim++)
                {
                    streamSystem << sitePositions(itSym, itdim) << "\t";
                }
                streamSystem << std::endl;
            }    
            streamSystem << "#SYMMETRIES_END" << std::endl;       
        }
        
// II - ENERGY PARAMETERS
    // II-1 OPPOSITE SPIN INTERACTION
        switch(UInteractionType)
        {
            case FHOppositeSpinInteractionType::Unone :
            {
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteConstant:
            {
                streamSystem << "#ENERGYPARAMETER_U_CONSTANT" << std::endl;
                streamSystem << UMatrix(0,0).real() << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteGivenDiagonal:
            {
                streamSystem << "#ENERGYPARAMETER_U_BYSITE" << std::endl;
                for(uint iSite = 0; iSite<sitePositions.n_rows; iSite++)
                {
                    streamSystem << UMatrix(iSite, iSite).real() << "\t";
                }
                streamSystem << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteBySiteType:
            {
                streamSystem << "#ENERGYPARAMETER_U_BYSITETYPE" << std::endl;
                for(uint iSiteType = 0; iSiteType < siteTypeNames.size(); iSiteType++)
                {
                    arma::uvec indSite = find(siteTypes = iSiteType, 1);
                    if (indSite.size() !=0)
                    {
                       streamSystem << siteTypeNames[iSiteType]<< "\t" << UMatrix(indSite(0)).real() << std::endl;
                    }            
                }
                streamSystem << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UGivenMatrix:
            {
                streamSystem << "#ENERGYPARAMETER_U_GIVENMATRIX" << std::endl;
                for(uint itU1 = 0; itU1 < UMatrix.n_rows; itU1++)
                {
                    for(uint itU2 = 0; itU2 < UMatrix.n_cols; itU2++)
                    {       
                        streamSystem << std::real(UMatrix(itU1, itU2)) << "\t";
                    }
                    streamSystem << std::endl;
                }
                break;
            }
            default :
            {
                break;
            }
        }



        // II-2 TUNNELING
        switch(tunnelingType)
        {
            case FHTunnelingType::tHnone:
            {
                break;
            }
            case FHTunnelingType::nearestNeighboursConstant:
            {
                streamSystem << "#ENERGYPARAMETER_tH_NEARESTNEIGHBORS" << std::endl;
                arma::uvec q1 = find(tMatrix,1);
                //std::cout << q1 << std::endl;
                if (q1.n_rows !=0)
                {
                    double th = std::real( tMatrix(q1(0)) );
                    streamSystem << th << std::endl;
                }
                else
                {
                     streamSystem << 1 << std::endl;
                }
                break;
            }
            case FHTunnelingType::tunnelingGivenMatrix:
            {
                streamSystem << "#ENERGYPARAMETER_tH_GIVENMATRIX" << std::endl;
                for(uint itt1 = 0; itt1 < tMatrix.n_rows; itt1++)
                {
                    for(uint itt2 = 0; itt2 < tMatrix.n_cols; itt2++)
                    {       
                        streamSystem << std::real(tMatrix(itt1, itt2)) << "\t";
                    }
                    streamSystem << std::endl;
                }
                break;
            }
            default :
            {
                break;
            }
        }



// III - INITIAL CONDITIONS
        switch (initialCondition)
        {
            case FHInitialConditionType::quenched:
                streamSystem << "#INITIALCONDITIONS_QUENCHED" << std::endl;
                break;
            case FHInitialConditionType::givenDiagonal:
            {
                streamSystem << "#INITIALCONDITIONS_GIVENDIAGONAL" << std::endl;
                for (uint itposUp = 0; itposUp < initConditions_upSpin.size(); itposUp++)
                {
                    if (initConditions_upSpin[itposUp])
                    {
                        streamSystem << itposUp << "\t";
                    }
                }
                streamSystem << std::endl;
                for (uint itposDown = 0; itposDown < initConditions_downSpin.size(); itposDown++)
                {
                    if (initConditions_downSpin[itposDown])
                    {
                        streamSystem << itposDown << "\t";
                    }
                }
                streamSystem << std::endl;       
                break;
            }   

        }


        streamSystem.close();
    }  
    else
    {
        std::cout << "Stream not open." << std::endl;
    }
}


FHSystem::~FHSystem()
{

}






// Other functions





std::complex<double> U_os_nij(double U0, double alpha, double d)
{
    double u = 0;
    if (d==0.)
        u = U0;
    else if (alpha != -1)
    {      
        u = U0/pow(d,alpha);      
    }
    return std::complex<double>(u);
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