// EXAMPLE : ./makeSitesSystem <systemName> <systemType> <parameter1>
//           ./makeSitesSystem bigGraphene graphene 4

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "include/FHLattice.h"
#include "include/FHSystem.h"

#include "include/json/json.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;

using namespace std;

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

/// @brief Makes a simple 3D-grid of sites 
/// @param nx number of sites in the X dimension
/// @param ny number of sites in the Y dimension
/// @param nz number of sites in the Z dimension
/// @param a distance between nearest neighbours
/// @return a vector of site positions
vector<vector<double>> makeGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a)
{
    
    unsigned int nPoints = nx*ny*nz;
    vector<vector<double>> pos = vector<vector<double>>(nPoints);

    unsigned int counter = 0;
    for (unsigned int z = 0; z < nz; z++)
    {
        for (unsigned int y = 0; y < ny; y++)
        {
            for (unsigned int x = 0; x < nx; x++)
            {
                vector<double> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
            }
        }
    }

    return pos;
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

/// @brief 
/// @param word 
/// @param vector 
/// @return 
bool isIn(string word, vector<string>& wordList)
{
    bool isInB = false;
    vector<string>::iterator itWord = wordList.begin();
    while ( (!isInB) && (itWord != wordList.end()) )
    {
        if (*itWord == word)  isInB = true;
        itWord ++;
    }

    return isInB;
}




// MAIN

int main(int argc, char* argv[])
{





    string fileName = "newSystem";
    string fileFolder = "systemsFolder";

    // Possible values
    vector<string> possibleGeometry = {"grid", "graphene", "bcc", "fcc"};
    enum geometry{
        grid,
        graphene,
        bcc,
        fcc
    };

    vector<string> possibleIC = {"ic_sw", "ic_cw", "ic_quenched"};
    enum IC{
        spinWave,
        chargeWave,
        quenched
    };

    vector<string> possibleNRJParameters = {"HF", "tJ"};
    enum NRJ{
        HF,
        tJ
    };

// Default values
    geometry geom = geometry::grid;                     
    vector<vector<double>> pos = makeGrid(2,1,1,1);     //<- default site system

    IC initC = IC::spinWave;                            //<- default initial condition

    NRJ en = NRJ::HF;
    vector<double> energyParam = {-1,1};                 //<- default energy parameters



    enum whatToRead{                                //<- what are we reading next
        nothing,
        geometryParameters,
        energyParameters
    };                  
    vector<double> geomParam = vector<double>();    //<- list of geometry parameters stored in double
    double argDouble = 0;                           //<- the current entry in double
    whatToRead wtr = whatToRead::nothing;
    bool isANumber = true;                          //<- is the current entry a number

    if (argc < 2)
    {
        std::cout << "To make a new system :" << std::endl;
        std::cout << "./makeSitesSystem <nameOfTheSystem> <geometryOfTheSites> <param1> <param2> ..." << std::endl;
    }
    else
    {

    cout << "== Writing a new FH system." << endl;
    for (int iArg = 1; iArg < argc ; iArg++)        //<- for each entry
    {

        isANumber = true;
        try                                         //<- We check if it's a number with a try/catch block
        {
            argDouble = stod(argv[iArg]);
        }  catch(const std::invalid_argument& ia) 
        {
            isANumber = false;
        }

        if (isANumber && (wtr == whatToRead::geometryParameters))   //<- if it's a number and we are reading the geometry parameters
        {
            geomParam.push_back(argDouble);                         //<- we store the number in the goemetry parameter vector
        }
        else if (isIn(string(argv[iArg]), possibleGeometry))        //<- else, if the word is among the geometry types
        {
            if (string(argv[iArg]) == possibleGeometry[0])  geom = geometry::grid;       //<- we note the geometry
            if (string(argv[iArg]) == possibleGeometry[1])  geom = geometry::graphene; 
            if (string(argv[iArg]) == possibleGeometry[2])  geom = geometry::bcc;  
            if (string(argv[iArg]) == possibleGeometry[3])  geom = geometry::fcc;
            wtr = whatToRead::geometryParameters;                   //<- We stand ready to read the geometry parameters
        } 
        else if (isIn(string(argv[iArg]), possibleIC)) //((argv[iArg]==string("ic_sw")) || (argv[iArg]==string("ic_cw")))    //<- else, if the word is among the Initial conditions types
        {
            if (string(argv[iArg]) == possibleIC[0])  initC = IC::spinWave; 
            if (string(argv[iArg]) == possibleIC[1])  initC = IC::chargeWave;           
            if (string(argv[iArg]) == possibleIC[2])  initC = IC::quenched;       //<- we note the initial condition
            wtr = whatToRead::nothing;                  //<- We stop reading geometry parameters
        }
        else 
        {
             wtr = whatToRead::nothing;                 //<- Else we don't stand ready for anything
             if (iArg==1)                           //<- and if it's the first argumentm that's the name of the file
             {
                fileName = argv[1];
             }
        }

    }

// Computation of the geometry
    for(unsigned int iGeom = geomParam.size(); iGeom<4; iGeom++)
    {
        geomParam.push_back(1.);
    }
    switch (geom)
    {
    case geometry::grid:
        pos = makeGrid( static_cast<int>(geomParam[0]), static_cast<int>(geomParam[1]), static_cast<int>(geomParam[2]), geomParam[3]);
        break;
    case geometry::graphene :
        pos = makeHoneyCombs(static_cast<int>(geomParam[0]),geomParam[1]);
        break;
    case geometry::bcc :
        pos = makeBCCGrid(static_cast<int>(geomParam[0]),geomParam[1], geomParam[2], 1.);
        break;
    case geometry::fcc :
        pos = makeFCCGrid(static_cast<int>(geomParam[0]),geomParam[1], geomParam[2], 1);
        break;
    }



// MAKING OF THE FHSystem
    arma::dmat sitePosition = armaMatfromVectorsVector(pos);
    //cout << "Position : " << endl << sitePosition << endl; 



    FHSystem newSystem = FHSystem();
    newSystem.setName(fileName);                                                                        //<- Name
    newSystem.setPositions(sitePosition);                                                               //<- Position
    newSystem.setTunnelingType(FHSystem::FHTunnelingType::nearestNeighboursConstant);                    //<- Energy
    newSystem.setTunneling_BetweenNeighbours(energyParam[0]);
    newSystem.setOppositeSpinInteractionType(FHSystem::FHOppositeSpinInteractionType::UonSiteConstant);
    newSystem.setUInteraction_OnSite(energyParam[1]);
    
    switch (initC)
    {
    case IC::spinWave:
        newSystem.setInitialConditions(FHSystem::FHInitialConditionType::spinWave);
        break;
    case IC::chargeWave:
        newSystem.setInitialConditions(FHSystem::FHInitialConditionType::chargeWave);
        break;
    case IC::quenched:
        newSystem.setInitialConditions(FHSystem::FHInitialConditionType::quenched);
        break; 
    default:
        newSystem.setInitialConditions(FHSystem::FHInitialConditionType::spinWave);
        break;
    }
                                                       //<- Initial conditions
    //cout << "file folder : " <<  fileFolder << endl; 


    Json::Value simRoot;
    simRoot["sim_name"] = fileName;
    simRoot["method"] = "HF";

    Json::Value param;
    param["t_max"] = 2.;
    param["dt"] = 0.01;
    param["stepsPerTimePoints"] = 1;
    param["n_traj"] = 1000;
    param["int_order"] = 4;
    simRoot["parameters"] = param;

    simRoot["system"] = newSystem.writeFHSystemJson();

    Json::Value res;
    res["diagDMatrix"] = true;
    res["diag2DMatrix"] = false;
    res["wholeDMatrix"] = false;
    res["momenta"] = false;
    res["energies"] = true;
    simRoot["results"] = res;

    
    //std::cout << simRoot << std::endl;

    std::string  fileAdress = fileFolder + "/" + fileName;
    std::string fileAdressJSON = fileAdress + ".json";

    std::ofstream streamSystem(fileAdressJSON);
    if (streamSystem.is_open())
    {
        streamSystem << simRoot ;
        streamSystem.close();
    }

    newSystem.writePositions(fileName, fileFolder);

   //std::cout << "In system, file adress : "<< fileAdress << std::endl;

//    newSystem.writeFHSystemBis(fileFolder);


    cout << "== " << fileName << " system created with " << pos.size() << " sites and " << newSystem.getUpSpinParticleNumber() + newSystem.getDownSpinParticleNumber() << " particles." << endl;  
    cout << "== Particle positions in " << fileFolder << "/" << fileName << ".xyz" << endl; 
    }
    return 0;
}