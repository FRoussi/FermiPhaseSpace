#include "FHSystem.h"
#include <cmath>
#include <iostream>
#include <fstream>

FHSystem::FHSystem()
{

    arma::dmat pos = arma::dmat(1,  1, arma::fill::zeros);

    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>(0,0);
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>(0,0);

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;

    initiateSystem(pos, initCond);
    computeMatrices();
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

    initiateSystem(pos, initCond);
    computeMatrices();

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

    initiateSystem(pos, initCond);
    computeMatrices();

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

    std::ifstream systFile (fileAdressName);
    if (systFile.is_open())
    {
        printf("== Reading %s \n", fileName.c_str());


        WhatToRead whatToR =  WhatToRead::nothing;  

        int nDim = -1;
        troolean orbitalHaveNames = troolean::dontknow;

        std::vector<unsigned int>       orbitalTypes = std::vector<unsigned int>();
        std::vector<std::string>       orbitalNames = std::vector<std::string>();
        std::vector<std::vector<double>> pos =   std::vector<std::vector<double>>();
        std::vector<unsigned int> initCond_up =   std::vector<unsigned int>();
        std::vector<unsigned int> initCond_down =   std::vector<unsigned int>();
        std::vector<unsigned int> particlesNumber(2, 0);
        

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
                whatToR = WhatToRead::diagonalIC_DownSpin;                      //<- At the end, the up-spin particles positions have been read.
                break;
            }
            case WhatToRead::diagonalIC_DownSpin :                   //<- Else, if the down-spin particles positions have not been read yet
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

        arma::dmat sitePos = arma::dmat(nDim, pos.size());
        for(unsigned int itPos = 0; itPos < pos.size(); itPos ++)
        {
            for (unsigned int itCoor = 0; itCoor<nDim; itCoor ++)
            {
                sitePos(itCoor, itPos) = pos[itPos][itCoor];
            }
        }

    if (orbitalHaveNames == troolean::yes)
    {
        
    }


    if (initialCondition == FHInitialConditionType::givenDiagonal)
    {
            std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
            initCond[0] = initCond_up;
            initCond[1] = initCond_down;

            initiateSystem(sitePos, initCond);
    }
    else if (initialCondition == FHInitialConditionType::quenched)
    {
        initiateSystem(sitePos, initialCondition, particlesNumber);
    }


    computeMatrices();

    }
    else
    {
        printf("Can't open %s", fileName.c_str());
    }


    
}


/// @brief setter - on-site interaction energy
/// @param U0 on-site interaction energy
void FHSystem::setInteraction_OnSite(double U0)
{
    UMatrix = U0 * arma::cx_dmat(sitePositions.n_cols, sitePositions.n_cols, arma::fill::eye);
}

/// @brief setter - nearest neighbout tunneling interaction energy
/// @param t tunneling interaction energy
void FHSystem::setTunneling_BetweenNeighbours(double t)
{
    tMatrix = arma::cx_dmat(sitePositions.n_cols, sitePositions.n_cols, arma::fill::zeros);
    for (unsigned int i = 0; i< sitePositions.n_cols; i++)
    {
        for (unsigned int j = 0; j< sitePositions.n_cols; j++)
        {
            if (neighbourMatrix[i][j])
            {
                tMatrix(i,j) = t;
            }
        }
    }
    printf("J-Matrix : \n");
    tMatrix.print();
}

/// @brief getter - number of sites
/// @return number of sites
unsigned int FHSystem::getnSites()
{
    return sitePositions.n_cols;
}

/// @brief getter - coordinate of a site
/// @param dim coordinate dimension
/// @param pos site number
/// @return coordinate of site 'pos'in dimension 'dim'
double FHSystem::getCoordinate(int dim, int pos)
{ 
    if ((dim < sitePositions.n_rows) && (pos < sitePositions.n_cols))
    {
        return sitePositions(dim, pos);
    } 
    else
    {
        return 0;
    }
}

/// @brief  getter - initial condition type
/// @return initial condition type
FHSystem::FHInitialConditionType     FHSystem::getInitialCondition()
{
    return initialCondition;
}

/// @brief getter - number of spin-up particles
/// @return number of spin-up particles
unsigned int            FHSystem::getUpSpinParticleNumber()
{
    return upSpinParticleNumber;
}

/// @brief getter - number of spin-down particles
/// @return number of spin-down particle
unsigned int            FHSystem::getDownSpinParticleNumber()
{
    return downSpinParticleNumber;
}

double FHSystem::getDistance(int i, int j)
{
    return distanceMatrix(i,j);
}

bool FHSystem::areNeighbour(int i, int j)
{
    return neighbourMatrix[i][j];
}

std::complex<double> FHSystem::getInteractionEnergy(int i, int j)
{
    return UMatrix(i,j);
}

std::complex<double> FHSystem::getTunnelingEnergy(int i, int j)
{
    return tMatrix(i,j);
}

std::complex<double> FHSystem::getExchangeCouplingEnergy(int i, int j)
{
    return JMatrix(i,j);
}

arma::cx_dmat* FHSystem::getUMatrixPointer()
{
    return &UMatrix;
}

arma::cx_dmat* FHSystem::getJMatrixPointer()
{
    return &JMatrix;
}

arma::cx_dmat* FHSystem::gettMatrixPointer()
{
    return &tMatrix;
}

std::string FHSystem::getName()
{
    return systemName;
}

std::string FHSystem::getFolder()
{
    return systemFolder;
}

std::vector<bool> FHSystem::getInitialCondition_upSpin()
{
    return initConditions_upSpin;
}

std::vector<bool> FHSystem::getInitialCondition_downSpin()
{
    return initConditions_downSpin;
}

/// @brief Set up the Fermi-Hubbard system from site positions and vectors of site number as initial condition
/// @param positions Matrix of site position
/// @param diagonalInitConditions vectors of sites on which there is a up-spin/down-spin particle
/// @return void
void FHSystem::initiateSystem(const arma::dmat& positions, const std::vector<std::vector<unsigned int>>& diagonalInitConditions )
{

    sitePositions = positions;

    initialCondition =  FHSystem::FHInitialConditionType::givenDiagonal;

    upSpinParticleNumber = diagonalInitConditions[0].size();
    downSpinParticleNumber = diagonalInitConditions[1].size();

    initConditions_upSpin = std::vector<bool>(sitePositions.n_cols);
    initConditions_downSpin = std::vector<bool>(sitePositions.n_cols);



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
/// @param positions Matrix of site position
/// @param initCondType a kind of initial condition (quenched, charge-wave, spin-wave ...)
/// @param partNumner vector of the numbers of up-spin, down-spin particle
/// @return void
void FHSystem::initiateSystem(const arma::dmat& positions, FHInitialConditionType initCondType, std::vector<unsigned int> partNumber )
{

    sitePositions = positions;

    initialCondition = FHSystem::FHInitialConditionType::quenched;

    upSpinParticleNumber = partNumber[0];
    downSpinParticleNumber = partNumber[1];


}

/// @brief Computes usefull matrices: distanceMatrix and neighbourMatrix.
/// @return void
void FHSystem::computeMatrices()
{

    unsigned int nSites = sitePositions.n_cols;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    sameSiteMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
    for (int i=0; i<nSites; i++)
    {
        for (int j=0; j<nSites; j++)
        {
            double dist = 0;
            for (int d = 0; d<sitePositions.n_rows; d++ )
            {
                dist = dist + pow(sitePositions(d, j) - sitePositions(d, i), 2);
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


/*bool FHSystem::isIn(T word, vector<T>& wordList)
{
    bool isInB = false;
    vector<string>::iterator itWord = wordList.begin();
    while ( (!isInB) && (itWord != wordList.end()) )
    {
        if (*itWord == word)  isInB = true;
        itWord ++;
    }

    return isInB;
}*/



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
