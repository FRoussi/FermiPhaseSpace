#include "FHSystem.h"


FHSystem::FHSystem()
{

    arma::dmat pos = arma::dmat(1,  1, arma::fill::zeros);

    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>();
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>();

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;

    setPositions(pos);
    setInitialConditions(initCond);
    computeMatrices();

    initiateEmptyEnergyMatrices();

    UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType = FHTunnelingType::nearestNeighbourConstant;

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
    setInitialConditions(initCond);
    computeMatrices();

    initiateEmptyEnergyMatrices();

    UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType = FHTunnelingType::nearestNeighbourConstant;


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
    setInitialConditions(initCond);
    computeMatrices();


    initiateEmptyEnergyMatrices();
    UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
    tunnelingType = FHTunnelingType::nearestNeighbourConstant;


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

    std::ifstream systFile (fileAdressName);
    if (systFile.is_open())
    {
        printf("== Reading %s \n", fileName.c_str());

        WhatToRead whatToR =  WhatToRead::nothing;  

        int nDim = -1;
        troolean orbitalHaveNames = troolean::dontknow;
        bool positionsReady = false;

        std::vector<unsigned int>           orbitalTypes    = std::vector<unsigned int>();
        std::vector<std::string>            orbitalNames    = std::vector<std::string>();
        std::vector<std::vector<double>>    pos             = std::vector<std::vector<double>>();
        std::vector<unsigned int>           initCond_up     = std::vector<unsigned int>();
        std::vector<unsigned int>           initCond_down   = std::vector<unsigned int>();
        std::vector<unsigned int>           particlesNumber(2, 0);

        std::vector<double> Jvector = std::vector<double>();
         std::vector<std::vector<double>> tHMatrix = std::vector<std::vector<double>>();
        std::vector<std::vector<double>> UMatrix = std::vector<std::vector<double>>();


        std::string line;
        while (std::getline(systFile, line))
        {
           if (line == "#POSITIONS")
           {
                whatToR = WhatToRead::positions;
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
                tunnelingType = FHTunnelingType::nearestNeighbourConstant; 
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
                    pos.push_back(coordinates);
                }
                else
                {
                    if (coordinates.size() == nDim)
                    {
                        pos.push_back(coordinates);
                    }
                    else
                    {
                        printf("Dimension problem with line \"%s\"", line.c_str());
                        return;
                    }
                }


                arma::dmat sitePos = armaMatfromVectorsVector(pos);

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
                break;
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
                if (UMatrix.size() == pos.size())
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
                if (tHMatrix.size() == pos.size())
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

            setInitialConditions(initCond);
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




/// @brief getter - number of sites
/// @return number of sites
unsigned int                                FHSystem::getnSites()
{
    return sitePositions.n_rows;
}

bool                                        FHSystem::sitesHaveNames()
{
        return (siteTypes.n_rows == sitePositions.n_rows);
}

arma::drowvec                               FHSystem::getSitePosition(int site)
{
    return sitePositions.row(site);
}

std::string                                 FHSystem::getSiteName(int site)
{
    return siteTypeNames[siteTypes[site]];
}

/// @brief getter - coordinate of a site
/// @param dim coordinate dimension
/// @param pos site number
/// @return coordinate of site 'pos'in dimension 'dim'
double                                      FHSystem::getCoordinate(int dim, int pos)
{ 
    if ((dim < sitePositions.n_cols) && (pos < sitePositions.n_rows))
    {
        return sitePositions(pos, dim);
    } 
    else
    {
        return 0;
    }
}

/// @brief  getter - initial condition type
/// @return initial condition type
FHSystem::FHInitialConditionType            FHSystem::getInitialCondition()
{
    return initialCondition;
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

double                                      FHSystem::getDistance(int i, int j)
{
    return distanceMatrix(i,j);
}

bool                                        FHSystem::areNeighbour(int i, int j)
{
    return neighbourMatrix[i][j];
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
    return std::complex<double>(0);
}

arma::cx_dmat*                              FHSystem::getUMatrixPointer()
{
    return &UMatrix;
}

arma::cx_dmat*                              FHSystem::getJMatrixPointer()
{
    return &JMatrix;
}

arma::cx_dmat*                              FHSystem::gettMatrixPointer()
{
    return &tMatrix;
}

arma::cx_dvec*                              FHSystem::getExtPotPointer()
{
    return &externalPotential;
}

std::string                                 FHSystem::getName()
{
    return systemName;
}

std::string                                 FHSystem::getFolder()
{
    return systemFolder;
}

std::vector<bool>                           FHSystem::getInitialCondition_upSpin()
{
    return initConditions_upSpin;
}

std::vector<bool>                           FHSystem::getInitialCondition_downSpin()
{
    return initConditions_downSpin;
}


// SETTERS

/// @brief setter - on-site interaction energy
/// @param U0 on-site interaction energy
void FHSystem::setUInteraction_OnSite(double U0, bool add)
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
void                    FHSystem::setUInteraction_GivenDiagonal(arma::dvec& U, bool add)
{
    if (U.n_rows == sitePositions.n_rows) 
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(arma::diagmat(U), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(arma::diagmat(U), arma::dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros));
        }      
    }
    else
    {
        std::cout << "Wrong size for the U interaction matrix." << std::endl;
    }
}

/// @brief setter - Matrix interaction energy
/// @param U Matrix interaction energy
void FHSystem::setUInteraction_GivenMatrix(arma::dmat& U, bool add)
{
    if ((U.n_rows == sitePositions.n_rows) && (U.n_cols == sitePositions.n_rows))
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(U, arma::dmat(U.n_rows, U.n_cols, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(U, arma::dmat(U.n_rows, U.n_cols, arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the U interaction matrix." << std::endl;
    }

}

/// @brief setter - nearest neighbout tunneling interaction energy
/// @param t tunneling interaction energy
void FHSystem::setTunneling_BetweenNeighbours(double t, bool add)
{
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
                tMatrix(i,j) += t;
            }
        }
    }
    printf("J-Matrix : \n");
    tMatrix.print();
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

void                    FHSystem::setExchangeCoupling(double J, bool add)
{
    JMatrix = arma::cx_dmat(sitePositions.n_rows, sitePositions.n_rows, arma::fill::zeros);
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

void                    FHSystem::setName(std::string name)    
{
    systemName = name;
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

void                    FHSystem::setPositions(const arma::dmat& positions)
{
    sitePositions = positions;
    computeMatrices();
}

/// @brief Set up the Fermi-Hubbard system from site positions and vectors of site number as initial condition
/// @param diagonalInitConditions vectors of sites on which there is a up-spin/down-spin particle
/// @return void
void FHSystem::setInitialConditions(const std::vector<std::vector<unsigned int>>& diagonalInitConditions )
{

    initialCondition =  FHSystem::FHInitialConditionType::givenDiagonal;

    upSpinParticleNumber = diagonalInitConditions[0].size();
    downSpinParticleNumber = diagonalInitConditions[1].size();

    initConditions_upSpin = std::vector<bool>(sitePositions.n_rows);
    initConditions_downSpin = std::vector<bool>(sitePositions.n_rows);



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
/// @param partNumner vector of the numbers of up-spin, down-spin particle
/// @return void
void FHSystem::setInitialConditions(FHInitialConditionType initCondType, std::vector<unsigned int> partNumber )
{

    initialCondition = initCondType;

    upSpinParticleNumber = partNumber[0];
    downSpinParticleNumber = partNumber[1];

}

/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void FHSystem::computeMatrices()
{

    unsigned int nSites = sitePositions.n_rows;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    sameSiteMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
    for (int i=0; i<nSites; i++)
    {
        for (int j=0; j<nSites; j++)
        {
            double dist = 0;
            for (int d = 0; d<sitePositions.n_cols; d++ )
            {
                dist = dist + pow(sitePositions(j, d) - sitePositions(i, d), 2);
            }
            dist = sqrt(dist);
            distanceMatrix(i, j) = dist;
            if (abs(dist)< 0.0001)
            {
                sameSiteMatrix[i][j] = true;
            }
            if (abs(dist - 1)< 0.0001)
            {
                neighbourMatrix[i][j] = true;
            }
        }
    }
}


void FHSystem::initiateEmptyEnergyMatrices()
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
}

/// @brief Write an ASCII version of the Fermi-Hubbard System
/// @param folderName emplacement of the folder
void     FHSystem::writeFHSystem(std::string folderName)
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
// I- POSITIONS OF SITES
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
        case FHTunnelingType::nearestNeighbourConstant:
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
