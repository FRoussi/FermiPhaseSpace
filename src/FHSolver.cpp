#include "../include/FHSolver.h"

FHSolver::FHSolver(FHSystem *fhS)
{
    fhSystem = fhS;

    computeMomenta = false;
    talkative = false;
    jsonFile = Json::Value();

    writeResultsTxt = false;
    writeResultsJson = false;
    writeResultsHDF5 = false;

    saveWholeDMatrix = false;
    saveDiagDMatrix = true;
    saveDiag2DMatrix = true;
    saveMomenta = false;
    saveEnergies = true;

}

void FHSolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
    talkative = false;
}

void                        FHSolver::setResultsWriter(std::string how, bool yn)
{
    if (how  == "txt")
    {
        writeResultsTxt = yn;
    }
    if (how  == "json")
    {
        writeResultsJson = yn;
    }
    if (how  == "hdf5")
    {
        writeResultsHDF5 = yn;
    }
}

void                        FHSolver::setResultsSaved(Json::Value res)
{
    saveDiagDMatrix = res.get("diagDMatrix", true).asBool();
    saveDiag2DMatrix = res.get("diag2DMatrix", false).asBool();
    saveWholeDMatrix = res.get("wholeDMatrix", false).asBool();
    saveMomenta = res.get("momenta", false).asBool();
    saveEnergies = res.get("energies", true).asBool();
}

void FHSolver::setSaveWholeDMatrix(bool yn)     {saveWholeDMatrix = yn;}
void FHSolver::setSaveDiagDMatrix(bool yn)      {saveDiagDMatrix = yn;}
void FHSolver::setSaveDiag2DMatrix(bool yn)     {saveDiag2DMatrix = yn;}
void FHSolver::setSaveMomenta(bool yn)          {saveMomenta = yn;}
void FHSolver::setSaveEnergies(bool yn)         {saveEnergies = yn;}

void FHSolver::setTimeParameters(double t_max, double d_t, uint nPerTS)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nPerTS;
    this->nTimePoints = floor(tMax / dt + 1);
}

void FHSolver::setTalkative(bool yn)
{
    talkative = yn;
}

void FHSolver::setComputeMomenta(bool yn)
{
    computeMomenta = yn;

}

void FHSolver::setJsonFile(Json::Value file)
{
    jsonFile = file;
}



/// @brief Initiate the ED solver, compute the initial wave-function, then the usefull operators.
///@return void
void FHSolver::initiateSolver()
{
 
}


/// @brief Simulation step, update the wave-function psi.
/// @param d_t Time step
/// @return void
void FHSolver::step(double d_t)
{
    
  
}


/// @brief Launch simulation and write the results
void FHSolver::runSimulation()
{

  

    

}


/// @brief Destructor
FHSolver::~FHSolver()
{
}
