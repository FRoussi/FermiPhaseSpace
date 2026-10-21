#pragma once

#include <bitset>
//#include <climits>
#include <iostream>
#include <string>
#include <vector>

//#include <cmath>
#include "json/json.h"

#include "FHSystem.h"
#include "ResultWriterTxt.h"
#include "ResultWriterJson.h"




class FHSolver  	
{
public:

    FHSolver(FHSystem* fhSystem);

    void                    setResultFile(std::string workName, std::string folderNameName);
    void                    setResultsWriter(std::string how, bool yn = true);
    void                    setTimeParameters(double t_max, double dt, uint nPerTS = 1);
    void                    setTalkative(bool yn);
    void                    setJsonFile(Json::Value file);

    void                    setComputeTBStates(bool yn);
    void                    setComputeMFStates(bool yn);

    void                    setResultsToSave(Json::Value res);
    void                    setSaveWholeDMatrix(bool yn);
    void                    setSaveDiagDMatrix(bool yn);
    void                    setSaveDiag2DMatrix(bool yn);
    void                    setSaveTBStates(bool yn);
    void                    setSaveMFStates(bool yn);
    void                    setSaveEnergies(bool yn);


    void        initiateSolver();


    void        step(double d_t);
    void        runSimulation();


    const arma::cx_double       jj = arma::cx_double(0,1);

    ~FHSolver();

protected:
// Fermi Hubbard System  :
    FHSystem*               fhSystem;                     //<- Pointer to the Fermi-Hubbard model
    Json::Value             jsonFile{Json::objectValue};                 //<- Json file containing the parameters of the simulation

// results in :
    bool                    writeResultsTxt{false};
    bool                    writeResultsJson{false};
    bool                    writeResultsHDF5{false};


    bool                    saveWholeDMatrix{false};
    bool                    saveDiagDMatrix{true};
    bool                    saveDiag2DMatrix{true};
    bool                    saveNiuNid{false};
    bool                    saveTBStates{false};
    bool                    saveMFStates{false};
    bool                    saveEnergies{true};

    uint                    nEnergies{4};   

    std::string             workName;                     //<- A name for the work, for the results file
    std::string             folderName;                   //<- A name for the work folder, for the results file

    //Momentum spaces :
    //bool                    computeTBStates;
    //bool                    computeMFStates;

// Work parameters
    bool                    talkative{false};

// Time parameters
    double                  tMax{1.};                         //<- Maximum time extend of the simulation 
    double                  dt{0.01};                           //<- Time step of calculated observables        
    unsigned int            nStepsPerTimePoint{1};           //<- Number of simulated time steps between the calculation of observable
    unsigned int            nTimePoints{100};                  //<- Number of time steps (=tMax/dt)
              
    
    std::vector<double> timeLine;
    arma::cx_dcube Nu;
    arma::cx_dcube Nd;
    arma::cx_dmat Niiu;
    arma::cx_dmat Niid;
    arma::cx_dmat Mkku;
    arma::cx_dmat Mkkd;
    arma::cx_dmat Mukku;
    arma::cx_dmat Mukkd;
    arma::cx_dmat Energies;
    arma::cx_dmat N2iiu, N2iid;
    arma::cx_dmat NiuNid;   
    arma::cx_dmat V_1p;
    arma::cx_dmat NiCH;

    
    void                        writeResults();


    


};

