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

    vector<string> possibleIC = {"ic_sw", "ic_cw", "ic_tbgs", "ic_mfgs"};
    enum IC{
        spinWave,
        chargeWave,
        tbGroundState,
        mfGroundState
    };

    vector<string> possibleNRJParameters = {"HF", "tJ"};
    enum NRJ{
        HF,
        tJ
    };

// Default values
    std::string latticeTypeName = "grid";             
    IC initC = IC::spinWave;                            //<- default initial condition
    //NRJ en = NRJ::HF;
    vector<double> energyParam = {-1,1};                 //<- default energy parameters
            
    vector<double> geomParam = vector<double>();    //<- list of geometry parameters stored in double
    double argDouble = 0;                           //<- the current entry in double
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
            if (iArg == 1)
            {
                fileName = argv[1];
            }
            else if (iArg == 2)
            {
                latticeTypeName = argv[2];
            }
            else if (isIn(string(argv[iArg]), possibleIC)) //((argv[iArg]==string("ic_sw")) || (argv[iArg]==string("ic_cw")))    //<- else, if the word is among the Initial conditions types
            {
                if (string(argv[iArg]) == possibleIC[0])  initC = IC::spinWave; 
                if (string(argv[iArg]) == possibleIC[1])  initC = IC::chargeWave;           
                if (string(argv[iArg]) == possibleIC[2])  initC = IC::tbGroundState;       //<- we note the initial condition
                if (string(argv[iArg]) == possibleIC[3])  initC = IC::mfGroundState;       //<- we note the initial condition
            }
            else
            {
                try                                         //<- We check if it's a number with a try/catch block
                {
                    argDouble = stod(argv[iArg]);
                }  catch(const std::invalid_argument& ia) 
                {
                    isANumber = false;
                }

                if (isANumber)   //<- if it's a number and we are reading the geometry parameters
                {
                    geomParam.push_back(argDouble);                         //<- we store the number in the goemetry parameter vector
                }
            }
        }


// Computation of the geometry
        for(unsigned int iGeom = geomParam.size(); iGeom<4; iGeom++)
        {
            geomParam.push_back(1.);
        }

        //FHLattice  lat0 = FHLattice(geom, geomParam);
        FHLattice  lat0 = FHLattice(latticeTypeName, geomParam);
        lat0.setLatticeName(fileName);

        //std::cout << "Lattice done" << std::endl;


// MAKING OF THE FHSystem

        FHSystem newSystem = FHSystem(lat0);
        newSystem.setName(fileName);    

        newSystem.setBasisInteractions();


        switch (initC)
        {
        case IC::spinWave:
            newSystem.setInitialConditions(FHSystem::FHInitialConditionType::spinWave);
            break;
        case IC::chargeWave:
            newSystem.setInitialConditions(FHSystem::FHInitialConditionType::chargeWave);
            break;
        case IC::tbGroundState:
            newSystem.setInitialConditions(FHSystem::FHInitialConditionType::tbGroundState);
            break; 
        case IC::mfGroundState:
            newSystem.setInitialConditions(FHSystem::FHInitialConditionType::mfGroundState);
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
        res["tbStates"] = false;
        res["mfStates"] = false;
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

        lat0.writePositions(fileName, fileFolder);

   //std::cout << "In system, file adress : "<< fileAdress << std::endl;

//    newSystem.writeFHSystemBis(fileFolder);


        cout << "== " << fileName << " system created with " << lat0.getnSites() << " sites and " << newSystem.getUpSpinParticleNumber() + newSystem.getDownSpinParticleNumber() << " particles." << endl;  
        cout << "== Particle positions in " << fileFolder << "/" << fileName << ".xyz" << endl; 
    }
    return 0;
}