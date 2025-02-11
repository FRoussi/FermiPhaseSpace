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
    Json::Value             jsonFile;           

// results in :
    bool                    writeResultsTxt;
    bool                    writeResultsJson;
    bool                    writeResultsHDF5;


    bool                    saveWholeDMatrix;
    bool                    saveDiagDMatrix;
    bool                    saveDiag2DMatrix;
    bool                    saveTBStates;
    bool                    saveMFStates;
    bool                    saveEnergies;

    std::string             workName;                     //<- A name for the work, for the results file
    std::string             folderName;                   //<- A name for the work folder, for the results file

    //Momentum spaces :
    //bool                    computeTBStates;
    //bool                    computeMFStates;

// Work parameters
    bool                    talkative;

// Time parameters
    double                  tMax;                         //<- Maximum time extend of the simulation 
    double                  dt;                           //<- Time step of calculated observables        
    unsigned int            nStepsPerTimePoint;           //<- Number of simulated time steps between the calculation of observable
    unsigned int            nTimePoints;                  //<- Number of time steps (=tMax/dt)
              
    

    
    void                        writeResults(arma::cx_dcube &Nu_t, arma::cx_dcube &Nd_t, 
        arma::cx_dmat &Niiu_t, arma::cx_dmat &Niid_t, 
        arma::cx_dmat &Mkku, arma::cx_dmat &Mkkd, 
        arma::cx_dmat &Mukku, arma::cx_dmat &Mukkd,                                             
        arma::cx_dmat &N2iiu, arma::cx_dmat &N2iid, 
        arma::cx_dmat &E,std::vector<double> &timeLine);


};

