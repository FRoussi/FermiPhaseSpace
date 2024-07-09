#include "../include/FHHFSolver.h"


FHHFSolver::FHHFSolver(FHSystem* fhS)
{
    fhSystem = fhS;
    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
    J = fhSystem->getJMatrixPointer();
    V = fhSystem->getVMatrixPointer();
    extPot = fhSystem->getExtPotPointer();


    computeMomenta = false;
    talkative = false;

    writeResultsTxt = false;
    writeResultsJson = false;
    writeResultsHDF5 = false;

}


void FHHFSolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
}

void                        FHHFSolver::writeResults(std::string how, bool yn)
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


/// @brief Set time simulation parameters.
/// @param t_max the maximum simulation time
/// @param d_t the time step between stored values in the simulation 
/// @param nSPTP number of simulation steps between two time steps
/// @return void
void FHHFSolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
}

void FHHFSolver::setTalkative(bool yn)
{
    talkative = yn;
}

/// @brief Initiate the HF solver, compute the initial density matrix.
void FHHFSolver::initiateSolver()
{

    FHSystem::FHInitialConditionType initCondition = fhSystem->getInitialCondition();
    integrator = HFIntegrator::RungeKutta4_38;

    currentSimulationTime = 0;

    unsigned int nSites = fhSystem->getnSites();


    if ((initCondition == FHSystem::FHInitialConditionType::givenDiagonal) || (initCondition == FHSystem::FHInitialConditionType::spinWave) || (initCondition == FHSystem::FHInitialConditionType::chargeWave))
    {

        Nu0 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        Nd0 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int it= 0; it < nSites; it ++)
        {
            Nu0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_upSpin()[it],   0) ;
            Nd0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_downSpin()[it], 0) ;
        }

        nPartU = std::floor(std::real(arma::accu(Nu0)));
        nPartD = std::floor(std::real(arma::accu(Nd0)));


        if (talkative)
        {
            std::cout << "E matrices" << std::endl;
            tM->print();
            U->print();
        }

        if (talkative)
        {
            printf("Test 4 - Initial Density matrices : \n");
            Nu0.print();
            Nd0.print();
        }


    }
    else if (initCondition == FHSystem::FHInitialConditionType::quenched)
    {

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();

        setComputeMomenta(true);

        arma::cx_dmat Mku0 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        arma::cx_dmat Mkd0 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

        for (unsigned int k=0; k<nPartU ; k++)
        {
            Mku0(k,k) = 1; 
        }

        for (unsigned int k=0; k<nPartD ; k++)
        {
            Mkd0(k,k) = 1; 
        }

        Nu0 =  V_1p * Mku0 * V_1p_1;
        Nd0 =  V_1p * Mkd0 * V_1p_1;

        // printf("Test 4 - Initial Density matrix : \n");
        // Nu0.print();
    }
    else
    {
       std::cout << "No initial condition found." << std::endl;
    }




}

void FHHFSolver::setComputeMomenta(bool yn)
{
    computeMomenta = yn;
    if (yn)
    {
        arma::dvec eigVal;
        arma::dmat Jreal = arma::real(*tM);
        arma::eig_sym(eigVal, V_1p, Jreal);
        V_1p_1 = V_1p.t();
    }
}

/// @brief Time step of a trajectory, computed with the 1st order Runge-Kutta method (Explicit Euler).
/// @param dt the time step
/// @return void
void FHHFSolver::step_RK1(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(currentSimulationTime, Nut, Ndt, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * dt ;
    Ndt = Ndt + dNd_dt * dt ;
}


/// @brief Time step of a trajectory, computed with the 2nd order Runge-Kutta method.
/// @param dt the time step
/// @return void
void FHHFSolver::step_RK2(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(currentSimulationTime, Nut, Ndt, dNu_dt, dNd_dt);

    arma::cx_dmat Nu_mid = Nut + dNu_dt * dt/2 ;
    arma::cx_dmat Nd_mid = Ndt + dNd_dt * dt/2 ;

    derivative(currentSimulationTime + dt/2, Nu_mid, Nd_mid, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * dt ;
    Ndt = Ndt + dNd_dt * dt ;
}


/// @brief Time step of a trajectory, computed with "original" 4th order Runge-Kutta method.
/// @param dt the time step
/// @return void
void FHHFSolver::step_RK4(double dt)
{
    arma::cx_dmat dNu_dt_k1 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k1 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k2 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k2 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k3 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k3 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k4 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k4 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

// k1 : 
    derivative(currentSimulationTime, Nut, Ndt, dNu_dt_k1, dNd_dt_k1);
// k2 : 0.5 | 0.5
    arma::cx_dmat Nu_k1 = Nut + dNu_dt_k1 * dt/2 ;
    arma::cx_dmat Nd_k1 = Ndt + dNd_dt_k1 * dt/2 ;
    derivative(currentSimulationTime + dt/2, Nu_k1 , Nd_k1, dNu_dt_k2, dNd_dt_k2);
// k3 : 0.5 |  0   0.5
    arma::cx_dmat Nu_k2 = Nut + dNu_dt_k2 * dt/2 ;
    arma::cx_dmat Nd_k2 = Ndt + dNd_dt_k2 * dt/2 ;
    derivative(currentSimulationTime + dt/2, Nu_k2 , Nd_k2, dNu_dt_k3, dNd_dt_k3);
// k3 : 1   |  0    0    1
    arma::cx_dmat Nu_k3 = Nut + dNu_dt_k3 * dt ;
    arma::cx_dmat Nd_k3 = Ndt + dNd_dt_k3 * dt ;
    derivative(currentSimulationTime + dt, Nu_k3 , Nd_k3, dNu_dt_k4, dNd_dt_k4);
// k4 :     |  1/6  2/6  2/6  1/6  
    Nut = Nut + (dNu_dt_k1 + 2* dNu_dt_k2 + 2*dNu_dt_k3 + dNu_dt_k4) * dt/6;
    Ndt = Ndt + (dNd_dt_k1 + 2* dNd_dt_k2 + 2*dNd_dt_k3 + dNd_dt_k4) * dt/6;
}


/// @brief Time step of a trajectory, computed with the "3/8 rule" 4th order Runge-Kutta method.
/// @param dt the time step
/// @return void
void FHHFSolver::step_RK4_38(double dt)
{
    arma::cx_dmat dNu_dt_k1 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k1 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k2 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k2 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k3 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k3 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat dNu_dt_k4 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt_k4 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

// k1 :
    derivative(currentSimulationTime, Nut, Ndt, dNu_dt_k1, dNd_dt_k1);
// k2 : 1/3  |  1/3
    arma::cx_dmat Nu_k = Nut + dNu_dt_k1 * dt/3 ;
    arma::cx_dmat Nd_k = Ndt + dNd_dt_k1 * dt/3 ;
    derivative(currentSimulationTime + dt/3, Nu_k , Nd_k, dNu_dt_k2, dNd_dt_k2);
// k3 : 2/3  |  -1/3  1
    Nu_k = Nut + (- dNu_dt_k1/3 + dNu_dt_k2) * dt ;
    Nd_k = Ndt + (- dNd_dt_k1/3 + dNd_dt_k2) * dt ;
    derivative(currentSimulationTime + 2*dt/3, Nu_k , Nd_k, dNu_dt_k3, dNd_dt_k3);
// k3 :  1   |  1     -1    1
    Nu_k = Nut + (dNu_dt_k1 - dNu_dt_k2 + dNu_dt_k3) * dt;
    Nd_k = Ndt + (dNd_dt_k1 - dNd_dt_k2 + dNd_dt_k3) * dt;
    derivative(currentSimulationTime + dt, Nu_k , Nd_k, dNu_dt_k4, dNd_dt_k4);
//    :      | 1/8    3/8   3/8   1
    Nut = Nut + (dNu_dt_k1 + dNu_dt_k2*3 + dNu_dt_k3*3 + dNu_dt_k4) * dt/8;
    Ndt = Ndt + (dNd_dt_k1 + dNd_dt_k2*3 + dNd_dt_k3*3 + dNd_dt_k4) * dt/8;
}


/// @brief ==NOT IMPLEMENTED== Time step of a trajectory, computed with the adaptive 4th/5th order Runge-Kutta method.
/// @param dt the time step
/// @return void
void FHHFSolver::step_RK45(double dt)
{
}


/// @brief Time derivative of the density matrix.
/// @param nut the up-spin density matrix.
/// @param ndt the down-spin density matrix.
/// @param dnut_dt the time derivative of the up-spin density matrix, to be updated.
/// @param dndt_dt the time derivative of the down-spin density matrix, to be updated.
/// @return void
void FHHFSolver::derivative(double currentTime, const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt)
{

    if (fhSystem->isSystemTimeDependant())
    {
        fhSystem->updateExternalPotential(currentTime);
    }

    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * ndt.diag()) + *extPot + arma::diagmat(*V * nut.diag());
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * nut.diag()) + *extPot + arma::diagmat(*V * ndt.diag());

    dnut_dt = - jj * ( nut * DeltaUp - DeltaUp * nut) ;
    dndt_dt = - jj * ( ndt * DeltaDo - DeltaDo * ndt) ;
}


/// @brief Launch simulation and write the results
void FHHFSolver::runSimulation( )
{


    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    arma::cx_dcube Nu = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube Nd = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    Nut = Nu0;
    Ndt = Nd0;

// Simulation :
    for (unsigned int t = 0;  t<nTimePoints; t++)
    {
        timeLine[t] = t*dt;
        Nu.slice(t) = Nut;
        Nd.slice(t) = Ndt;




        for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
        {
            if      (integrator==HFIntegrator::RungeKutta1)     step_RK1(dt/nStepsPerTimePoint);
            else if (integrator==HFIntegrator::RungeKutta2)     step_RK2(dt/nStepsPerTimePoint);
            else if (integrator==HFIntegrator::RungeKutta4)     step_RK4(dt/nStepsPerTimePoint);
            else if (integrator==HFIntegrator::RungeKutta4_38)  step_RK4_38(dt/nStepsPerTimePoint);
            else if (integrator==HFIntegrator::RungeKutta45)    step_RK45(dt/nStepsPerTimePoint);

            currentSimulationTime += dt/nStepsPerTimePoint;
        }
 
    }


// Computing other values :
    arma::cx_dmat Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat Mkku = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat Mkkd = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EU =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EtH =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec nPartU =  arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec nPartD =  arma::drowvec(nTimePoints, arma::fill::zeros);
    for (unsigned int t = 0;  t<nTimePoints; t++)
    {

        Nut = Nu.slice(t);
        Ndt = Nd.slice(t);
        Niiu.col(t) = Nut.diag();
        Niid.col(t) = Ndt.diag();
        nPartU(t) = std::real(arma::trace(Nut));
        nPartD(t) = std::real(arma::trace(Ndt));

        EtH(t) =  arma::accu(*tM % Nut) + arma::accu(*tM % Ndt);
        EU(t) = (arma::cdot(Nut.diag(), *U * Ndt.diag()) + arma::cdot(Ndt.diag(),  *U * Nut.diag())) / arma::cx_double(2,0);

        if (computeMomenta)
        {
            arma::cx_dmat Mu = V_1p_1 * Nut * V_1p;
            arma::cx_dmat Md = V_1p_1 * Ndt * V_1p;
            Mkku.col(t) = Mu.diag();
            Mkkd.col(t) = Md.diag();
        }

    }

    arma::cx_drowvec E = EU + EtH;

    if (talkative)
    {
        printf("--- Mean field dynamic, n1u : \n");
        (arma::real(Niiu.row(0))).print();
        printf("--- Mean field dynamic, E : \n");
        (arma::real(E)).print();
    }


    std::cout << "- Computation done, writing results." << std::endl;

       // Writing the results :
    if (writeResultsTxt)
    {
        ResultWriterTxt results = ResultWriterTxt(workName, folderName);
        results.newResults();

        results.writeLine("t", timeLine);
        for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
        {
            results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Niiu.row(iSite));
        }
        for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
        {
            results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_d", Niid.row(iSite));
        }
        if (computeMomenta)
        {
            for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
            {
                results.writeLineInReal("M_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Mkku.row(iSite));
            }
            for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
            {
                results.writeLineInReal("M_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_d", Mkkd.row(iSite));
            }

        }
        results.writeLineInReal("EU", EU);
        results.writeLineInReal("EtH", EtH);
        results.writeLineInReal("E", E);

        std::cout << "- Results writen in " << folderName << "/" << workName << ".txt" << std::endl;
    }

    if (writeResultsJson)
    {

        ResultWriterJson rJson; 
        Json::Value resJson = fhSystem->getOriginalJsonSystem();
        std::cout << "- Writing a json file." << std::endl;
        Json::Value niiuJson;
        for (uint iSite = 0 ; iSite < Niiu.n_rows ; iSite++)
        {
            std::string iSiteStr = std::to_string(iSite + 1);
            niiuJson[iSiteStr] = rJson.writeLineInReal(Niiu.row(iSite));
        }
        resJson["Niiu"] = niiuJson;
std::cout << "- Writing niiu." << std::endl;
        Json::Value niidJson;
        for (uint iSite = 0 ; iSite < Niid.n_rows ; iSite++)
        {
            std::string iSiteStr = std::to_string(iSite + 1);
            niidJson[iSiteStr] = rJson.writeLineInReal(Niid.row(iSite));
        }
        resJson["Niid"] = niidJson;


        Json::StreamWriterBuilder wbuilder;
        wbuilder["indentation"] = "";
        std::string document = Json::writeString(wbuilder, resJson);
std::cout << "- Writing the doc." << std::endl;
        std::string fileAdressName = folderName + "/" + workName + ".json";
        std::ofstream resultFile(fileAdressName);
        resultFile << document;
        resultFile.close();

        std::cout << "- Results writen in " << folderName << "/" << workName << ".json" << std::endl;
    }

    if (writeResultsHDF5)
    {
        
    }



    

}




/// @brief Destructor
FHHFSolver::~FHHFSolver()
{

}
 