#include "FHSystem.h"
#include <cmath>
#include <iostream>
#include <fstream>

FHSystem::FHSystem()
{
    nSites = 1;
    nDimentions = 1;
    sitePositions = arma::dmat(nDimentions,  nSites, arma::fill::zeros);
    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));

    computeMatrices();
}

FHSystem::FHSystem(int x)
{
    nSites = x;
    nDimentions = 1;
    sitePositions = arma::dmat(nDimentions,  nSites, arma::fill::zeros);
    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));
    for (int d=0; d<2; d++)
    {
        sitePositions[d] = (double)d;
    }

    computeMatrices();
}

FHSystem::FHSystem(int x, int y, int z)
{
    nDimentions = 3;
    nSites = x * y * z;

    sitePositions = arma::dmat(nDimentions,  nSites, arma::fill::zeros);
    distanceMatrix =  arma::dmat(nSites,  nSites, arma::fill::zeros);
    neighbourMatrix = std::vector<std::vector<bool>>(nSites,  std::vector<bool>(nSites, false));

    int tabIndex = 0;
    for (int i_z = 0; i_z < z; i_z++)
    {
        for (int i_y = 0; i_y < y; i_y++)
        {
            for (int i_x = 0; i_x < x; i_x++)
            {   
                sitePositions[tabIndex] = double(i_x);
                sitePositions[tabIndex+1] = double(i_y);
                sitePositions[tabIndex+2] = double(i_z); 
                tabIndex +=3;
            }
        }       
    }

    computeMatrices();
}

FHSystem::FHSystem(std::string fileName)
{
    std::ifstream systFile (fileName);
    if (systFile.is_open())
    {
        printf("Reading %s", fileName);

        std::string line;
        std::getline(systFile, line);
        while ()
        {
            systFile.getline(line);
        }




        systFile.close();
    }
    else
    {
        printf("Can't open %s", fileName);
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
    if ((dim < nDimentions) && (pos < nSites))
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

void FHSystem::initiateSystem(arma::dmat positions, std::vector<std::vector<unsigned int>> initConditions )
{

    nSites = positions.n_rows;
    nDimentions = positions.n_cols;

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

    
    for (int i=0; i<nSites; i++)
    {
        for (int j=0; j<nSites; j++)
        {
            double dist = 0;
            for (int d = 0; d<nDimentions; d++ )
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
