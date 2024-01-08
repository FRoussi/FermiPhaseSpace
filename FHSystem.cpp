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
    name = fileName;
    std::string fileAdressName =  fileName; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + fileName;
    }

    std::ifstream systFile (fileAdressName);
    if (systFile.is_open())
    {
        printf("Reading %s \n", fileName.c_str());

        bool readPosition = false;
        bool readInitCond = false;
        bool initCondUpSpinRead = false;
        bool initCondDownSpinRead = false;
        int nDim = -1;

        std::vector<std::vector<double>> pos =   std::vector<std::vector<double>>();
        std::vector<unsigned int> initCond_up =   std::vector<unsigned int>();
        std::vector<unsigned int> initCond_down =   std::vector<unsigned int>();

        std::string line;
        while (std::getline(systFile, line))
        {
           if (line == "#POSITIONS")
           {
                readPosition = true;
                readInitCond = false;
           }
           else if (line == "#INITIALCONDITIONS")
           {
                readPosition = false;
                readInitCond = true;
           }
           else if(readPosition)
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
           }
           else if(readInitCond)
           {
                std::stringstream ss(line);
                std::string siteN;
                
                if (!initCondUpSpinRead)
                {
                    while (std::getline (ss, siteN, ' ')) 
                    {
                        initCond_up.push_back(std::stoi(siteN)); 
                    }
                    initCondUpSpinRead = true;
                }
                else if (!initCondDownSpinRead)
                {
                    while (std::getline (ss, siteN, ' ')) 
                    {
                        initCond_down.push_back(std::stoi(siteN)); 
                    }
                    initCondDownSpinRead = true;
                }
                else
                {
                    printf("Too many initial conditions");
                }


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


        std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
        initCond[0] = initCond_up;
        initCond[1] = initCond_down;

        initiateSystem(sitePos, initCond);
        computeMatrices();

    }
    else
    {
        printf("Can't open %s", fileName.c_str());
    }


    
}



void FHSystem::setInteraction_OnSite(double U0)
{
    UMatrix = U0 * arma::cx_dmat(nSites, nSites, arma::fill::eye);
}

void FHSystem::setTunneling_BetweenNeighbours(double J)
{
    JMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    for (unsigned int i = 0; i< nSites; i++)
    {
        for (unsigned int j = 0; j< nSites; j++)
        {
            if (neighbourMatrix[i][j])
            {
                JMatrix(i,j) = J;
            }
        }
    }

    JMatrix.print();
}

unsigned int FHSystem::getnSites()
{
    return nSites;
}

double FHSystem::getCoordinate(int dim, int pos)
{ 
    if ((dim < nDimensions) && (pos < nSites))
    {
        return sitePositions(dim, pos);
    } 
    else
    {
        return 0;
    }
}

double FHSystem::getDistance(int i, int j)
{
    return distanceMatrix(i,j);
}

bool FHSystem::areNeighbour(int i, int j)
{
    return neighbourMatrix[i][j];
}

arma::cx_double FHSystem::getInteractionEnergy(int i, int j)
{
    return UMatrix(i,j);
}

arma::cx_double FHSystem::getTunnelingEnergy(int i, int j)
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
    return name;
}

void FHSystem::initiateSystem(arma::dmat positions, std::vector<std::vector<unsigned int>> initConditions )
{

    nDimensions = positions.n_rows;
    nSites = positions.n_cols;

    sitePositions = positions;

    initConditions_upSpin = std::vector<bool>(nSites);
    initConditions_downSpin = std::vector<bool>(nSites);

    for (unsigned int itSite = 0; itSite < initConditions[0].size(); itSite++)
    {
        initConditions_upSpin[itSite] = true;
    }
    for (unsigned int itSite = 0; itSite < initConditions[1].size(); itSite++)
    {
        initConditions_downSpin[itSite] = true;
    }

}

void FHSystem::computeMatrices()
{
    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    
    for (int i=0; i<nSites; i++)
    {
        for (int j=0; j<nSites; j++)
        {
            double dist = 0;
            for (int d = 0; d<nDimensions; d++ )
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
