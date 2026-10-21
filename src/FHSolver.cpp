#include "../include/FHSolver.h"



FHSolver::FHSolver(FHSystem *fhS) :
    fhSystem(fhS),
    jsonFile(Json::Value())
{

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

void                        FHSolver::setResultsToSave(Json::Value res)
{
    saveDiagDMatrix     = res.get("diagDMatrix", true).asBool();
    saveDiag2DMatrix    = res.get("diag2DMatrix", false).asBool();
    saveNiuNid          = res.get("onSiteDoublons", false).asBool();
    saveWholeDMatrix    = res.get("wholeDMatrix", false).asBool();
    saveEnergies        = res.get("energies", true).asBool();
    saveTBStates        = res.get("tbStates", false).asBool();
    saveMFStates        = res.get("mfStates", false).asBool();


}


void                    FHSolver::setComputeTBStates(bool yn)
{
    saveTBStates = yn;
    if (saveTBStates && (!fhSystem->isTBGroundStateComputed()))
    {
        fhSystem->computeTBGroundState();
    }
}

void                    FHSolver::setComputeMFStates(bool yn)
{
    saveMFStates = yn;
    if (saveMFStates && (!fhSystem->isMFGroundStateComputed()))
    {
        fhSystem->computeMFGroundState();
    }
}

void FHSolver::setSaveWholeDMatrix(bool yn)     {saveWholeDMatrix = yn;}
void FHSolver::setSaveDiagDMatrix(bool yn)      {saveDiagDMatrix = yn;}
void FHSolver::setSaveDiag2DMatrix(bool yn)     {saveDiag2DMatrix = yn;}
void FHSolver::setSaveEnergies(bool yn)         {saveEnergies = yn;}

void FHSolver::setSaveTBStates(bool yn)         
{
    saveTBStates = yn;
    if (saveTBStates)
    {
        if (!fhSystem->isTBGroundStateComputed())
        {
            fhSystem->computeTBGroundState();
        }
    }
}

void FHSolver::setSaveMFStates(bool yn)   
{
    saveMFStates = yn;
    if (saveMFStates)
    {
        if (!fhSystem->isMFGroundStateComputed())
        {
            fhSystem->computeMFGroundState();
        }
    }
}

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


void FHSolver::setJsonFile(Json::Value file)
{
    jsonFile = file;
}



void FHSolver:: writeResults()
{

    // Writing the results :
    if (writeResultsTxt)
    {
        ResultWriterTxt results = ResultWriterTxt(workName, folderName);
        results.newResults();
        results.writeLine("t", timeLine);
        if (saveDiagDMatrix)
        {
            for (unsigned int iSite=0; Niiu.n_rows; iSite++)
            {
                results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Niiu.row(iSite));
            }
            for (unsigned int iSite=0; Niid.n_rows; iSite++)
            {
                results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_d", Niid.row(iSite));
            }
        }

        if (saveTBStates)
        {
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("M_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_u", Mkku.row(kMom));
            }
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("M_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_d", Mkkd.row(kMom));
            }
        }
        if (saveMFStates)
        {
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("Mu_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_u", Mukku.row(kMom));
            }
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("Mu_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_d", Mukkd.row(kMom));
            }
        }
        if (saveEnergies)
        {
            results.writeLineInReal("EU", Energies.row(0));
            results.writeLineInReal("EtH", Energies.row(1));
            results.writeLineInReal("Epe", Energies.row(2));
            results.writeLineInReal("E", Energies.row(3));
        }


         std::cout << "   -- Results written in " << folderName  << workName << ".txt" << std::endl;
    }

    if (writeResultsJson)
    {

        ResultWriterJson rJson; 
        Json::Value resJson = jsonFile;

        resJson["fpsVersion"] = fhSystem->getFPSVersion(); 
        if (fhSystem->isTBGroundStateComputed())
        {
            resJson["tbGroundState"] = fhSystem->writeTBGroundStateJson();
        }
        if (fhSystem->isMFGroundStateComputed())
        {
            resJson["mfGroundState"] = fhSystem->writeMFGroundStateJson();
        }

        Json::Value simJson = Json::Value();

        simJson["t"] = rJson.writeLine(timeLine);

        if (saveWholeDMatrix)
        {
            //Json::Value n12uJS(Json::arrayValue);
            //n12uJS.append(rJson.writeLineInReal(Nu.tube(1,2)));

            simJson["N1hu"] = rJson.writeLine(Nu.tube(0,Nu.n_cols-1));
            simJson["N1hd"] = rJson.writeLine(Nd.tube(0,Nd.n_cols-1));
        }

        if (saveDiagDMatrix)
        {
            Json::Value niiuJson(Json::arrayValue);
            for (uint iSite = 0 ;  iSite < Niiu.n_rows ; iSite++)
            {
                niiuJson.append(rJson.writeLineInReal(Niiu.row(iSite)));
            }
            simJson["Niiu"] = niiuJson;
            Json::Value niidJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Niid.n_rows ; iSite++)
            {
                niidJson.append(rJson.writeLineInReal(Niid.row(iSite)));
            }
            simJson["Niid"] = niidJson;
        }

        if (saveDiag2DMatrix)
        {
            Json::Value n2iiuJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < N2iiu.n_rows ; iSite++)
            {
                n2iiuJson.append(rJson.writeLineInReal(N2iiu.row(iSite)));
            }
            simJson["N2iiu"] = n2iiuJson;
            Json::Value n2iidJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < N2iid.n_rows ; iSite++)
            {
                n2iidJson.append(rJson.writeLineInReal(N2iid.row(iSite)));
            }
            simJson["N2iid"] = n2iidJson;
        }

        if (saveNiuNid)
        {
            Json::Value nundJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < NiuNid.n_rows ; iSite++)
            {
                nundJson.append(rJson.writeLineInReal(NiuNid.row(iSite)));
            }
            simJson["NuNd"] = nundJson;
        }

        if (saveTBStates)
        {
            Json::Value mkkuJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Mkku.n_rows ; iSite++)
            {
                mkkuJson.append(rJson.writeLineInReal(Mkku.row(iSite)));
            }
            simJson["Mkku"] = mkkuJson;
            Json::Value mkkdJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Mkkd.n_rows ; iSite++)
            {
                mkkdJson.append(rJson.writeLineInReal(Mkkd.row(iSite)));
            }
            simJson["Mkkd"] = mkkdJson;
        }
        if (saveMFStates)
        {
            Json::Value mukkuJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Mukku.n_rows ; iSite++)
            {
                mukkuJson.append(rJson.writeLineInReal(Mukku.row(iSite)));
            }
            simJson["Mukku"] = mukkuJson;
            Json::Value mukkdJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Mukkd.n_rows ; iSite++)
            {
                mukkdJson.append(rJson.writeLineInReal(Mukkd.row(iSite)));
            }
            simJson["Mukkd"] = mukkdJson;
        }

        if (saveEnergies)
        {

            Json::Value EJson(Json::objectValue);
            EJson["tH"] = rJson.writeLineInReal(Energies.row(0));
            EJson["U"] = rJson.writeLineInReal(Energies.row(1));
            EJson["V"] = rJson.writeLineInReal(Energies.row(2));
            EJson["pe"] = rJson.writeLineInReal(Energies.row(3));
            EJson["tot"] = rJson.writeLineInReal(Energies.row(nEnergies));
            simJson["E"] = EJson;
        }

        if (fhSystem->getSystemType() == FHSystemType::SysCoreHole)
        {
            Json::Value NiCHJson(Json::objectValue);
            NiCHJson["u"] = rJson.writeLine(NiCH.row(0));
            NiCHJson["d"] = rJson.writeLine(NiCH.row(1));
            simJson["NiCh"] = NiCHJson;
        }

        resJson["simulation"] = simJson;


        Json::StreamWriterBuilder wbuilder;
        wbuilder["indentation"] = "\t";
        std::string document = Json::writeString(wbuilder, resJson);
        std::string fileAdressName = folderName + workName + ".json";
        std::ofstream resultFile(fileAdressName);
        resultFile << document;
        resultFile.close();



        std::cout << "   -- Results written in " << folderName << workName << ".json" << std::endl;
    }

    if (writeResultsHDF5)
    {
        
    }
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
