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

        // 0 : nothing
        // 1 : positions
        // 2 : diagonal initial conditions, up spin
        // 3 : diagonal initial conditions, down spin
        // 4 : number of particles, up spin
        // 5 : number of particles, down spin
        unsigned int whatToRead = 0;    


        int nDim = -1;

        std::string         initialConditionType;
        std::vector<std::vector<double>> pos =   std::vector<std::vector<double>>();
        std::vector<unsigned int> initCond_up =   std::vector<unsigned int>();
        std::vector<unsigned int> initCond_down =   std::vector<unsigned int>();
        std::vector<unsigned int> particlesNumber(2, 0);
        

        std::string line;
        while (std::getline(systFile, line))
        {
           if (line == "#POSITIONS")
           {
                whatToRead = 1;
           }
           else if (line == "#INITIALCONDITIONS_GIVENDIAGONAL")
           {
                initialConditionType = "givenDiagonal";
                whatToRead = 2;
           }
           else if (line == "#INITIALCONDITIONS_QUENCHED")
           {
                initialConditionType = "quenched";
                whatToRead = 4;
           }
           else switch (whatToRead)
           {
           case 1:  // Read positions
            {
                std::stringstream ss(line);
                std::string coord;
                std::vector<double> coordinates = std::vector<double>();

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
           case 2:  // Read  up spin particles
           {
                std::stringstream ss(line);
                std::string siteN;
                
                while (std::getline(ss, siteN, ' '))            //<- We read the line word by word
                {
                        //std::cout << "Reading initial positions of up-spin particles, siteN = " << siteN << std::endl;
                    initCond_up.push_back(std::stoi(siteN));    //<- And fill the inital condition vector with the site numbers
                    particlesNumber[0]++;
                }
                whatToRead = 3;                      //<- At the end, the up-spin particles positions have been read.
                break;
            }
            case 3:                   //<- Else, if the down-spin particles positions have not been read yet
            {
                std::stringstream ss(line);
                std::string siteN;
                while (std::getline(ss, siteN, ' ')) 
                {
                       // std::cout << "Reading initial positions of down-spin particles, siteN = " << siteN << std::endl;
                    initCond_down.push_back(std::stoi(siteN)); 
                    particlesNumber[1]++;
                }
                whatToRead = 0; 
                break;
            }
            case 4:
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[0] = std::stoi(nPart);
                whatToRead = 5;
            }
            break;
            case 5:
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[1] = std::stoi(nPart);
                whatToRead = 0;
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


    if (initialConditionType == std::string("givenDiagonal"))
    {
            std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
            initCond[0] = initCond_up;
            initCond[1] = initCond_down;

            initiateSystem(sitePos, initCond);
    }
    else if (initialConditionType == std::string("quenched"))
    {
        initiateSystem(sitePos, initialConditionType, particlesNumber);
    }


    computeMatrices();

    }
    else
    {
        printf("Can't open %s", fileName.c_str());
    }


    
}



void FHSystem::setInteraction_OnSite(double U0)
{
    UMatrix = U0 * arma::cx_dmat(sitePositions.n_cols, sitePositions.n_cols, arma::fill::eye);
}

void FHSystem::setTunneling_BetweenNeighbours(double J)
{
    JMatrix = arma::cx_dmat(sitePositions.n_cols, sitePositions.n_cols, arma::fill::zeros);
    for (unsigned int i = 0; i< sitePositions.n_cols; i++)
    {
        for (unsigned int j = 0; j< sitePositions.n_cols; j++)
        {
            if (neighbourMatrix[i][j])
            {
                JMatrix(i,j) = J;
            }
        }
    }
    printf("J-Matrix : \n");
    JMatrix.print();
}

unsigned int FHSystem::getnSites()
{
    return sitePositions.n_cols;
}

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

std::string             FHSystem::getInitialCondition()
{
    return initialCondition;
}

    unsigned int            FHSystem::getUpSpinParticleNumber()
    {
        return upSpinParticleNumber;
    }
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


void FHSystem::initiateSystem(const arma::dmat& positions, const std::vector<std::vector<unsigned int>>& diagonalInitConditions )
{

    sitePositions = positions;

    initialCondition = "givenDiagonal";

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

void FHSystem::initiateSystem(const arma::dmat& positions, const std::string initCondition, std::vector<unsigned int> partNumber )
{

    sitePositions = positions;

    initialCondition = initCondition;

    upSpinParticleNumber = partNumber[0];
    downSpinParticleNumber = partNumber[1];


}

void FHSystem::computeMatrices()
{

    unsigned int nSites = sitePositions.n_cols;

    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
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
            if (abs(dist - 1)< 0.0001)
            {
                neighbourMatrix[i][j] = true;
            }
            else
            {
                neighbourMatrix[i][j] = false;
            }
        }
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
