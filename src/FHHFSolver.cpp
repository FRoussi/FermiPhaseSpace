#include "../include/FHHFSolver.h"


FHHFSolver::FHHFSolver(FHSystem* fhS) : FHSolver(fhS)
{

/*
    fhSystem = fhS;

    computeMomenta = false;
    talkative = false;
    jsonFile = Json::Value();

    writeResultsTxt = false;
    writeResultsJson = false;
    writeResultsHDF5 = false; */

    currentSimulationTime = 0;
    integratorOrder = 2;
    integrator = HFIntegrator::RungeKutta2;
    

    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
    //J = fhSystem->getJMatrixPointer();
    V = fhSystem->getVMatrixPointer();
    extPot = fhSystem->getExtPotPointer();

}




void                        FHHFSolver::setIntegratorOrder(int order)
{
    switch (order)
    {
        case 1:
        {
            integratorOrder = 1;
            integrator = HFIntegrator::RungeKutta1;
            break;
        }
        case 2:
        {
            integratorOrder = 2;
            integrator = HFIntegrator::RungeKutta2;
            break;
        }
        case 4:
        {
            integratorOrder = 4;
            integrator = HFIntegrator::RungeKutta4_38;
            break;
        }
        default:
        {
            std::cout << "Integrator order " << order << " not implemented. Switching to default 2." << std::endl;
            integratorOrder = 2;
            integrator = HFIntegrator::RungeKutta2;
            break;
        }
    }

}

void                        FHHFSolver::setMethodParameters(Json::Value param)
{
    double param_dt = param.get("dt", 0.01).asDouble();
    double param_tm = param.get("t_max", 2).asDouble();
    int param_sptp = param.get("stepsParTimePoints", 1).asInt();
    int param_intOrder = param.get("int_order", 2).asInt();
    setTimeParameters(param_tm, param_dt, param_sptp);
    setIntegratorOrder(param_intOrder);
}

/// @brief Initiate the HF solver, compute the initial density matrix.
void FHHFSolver::initiateSolver()
{

    std::cout << "   -- HF Solver initiation." << std::endl;

    FHSystem::FHInitialConditionType initCondition = fhSystem->getInitialCondition();
    integrator = HFIntegrator::RungeKutta4_38;
    //integrator = HFIntegrator::RungeKutta1_Mem;

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
    else if (initCondition == FHSystem::FHInitialConditionType::tbGroundState)
    {

        if (!fhSystem->isTBGroundStateComputed())
        {
            fhSystem->computeTBGroundState();
        }

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();

        setComputeMomenta(true);

        Nu0 =  arma::cx_dmat(*fhSystem->get_tbgs_Nu_Pointer());
        Nd0 =  arma::cx_dmat(*fhSystem->get_tbgs_Nd_Pointer());

        // printf("Test 4 - Initial Density matrix : \n");
        // Nu0.print();
    }
    else if (initCondition == FHSystem::FHInitialConditionType::mfGroundState)
    {

        if (!fhSystem->isMFGroundStateComputed())
        {
            fhSystem->computeMFGroundState();
        }

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();

        setComputeMomenta(true);

        Nu0 =  arma::cx_dmat(*fhSystem->get_mfgs_Nu_Pointer());
        Nd0 =  arma::cx_dmat(*fhSystem->get_mfgs_Nd_Pointer());

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
        if (!fhSystem->isTBGroundStateComputed())
        {
            fhSystem->computeTBGroundState();
        }
    }

}

/// @brief Time step of a trajectory, computed with the 1st order Runge-Kutta method (Explicit Euler).
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK1(double d_t)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(currentSimulationTime, Nut, Ndt, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * d_t ;
    Ndt = Ndt + dNd_dt * d_t ;
}

/// @brief Time step of a trajectory, computed with the 1st order Runge-Kutta method (Explicit Euler).
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK1_Mem(double d_t, arma::cx_dmat & dNu_dt, arma::cx_dmat & dNd_dt)
{

    derivative(currentSimulationTime, Nut, Ndt, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * d_t ;
    Ndt = Ndt + dNd_dt * d_t ;
}


/// @brief Time step of a trajectory, computed with the 2nd order Runge-Kutta method.
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK2(double d_t)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(currentSimulationTime, Nut, Ndt, dNu_dt, dNd_dt);

    arma::cx_dmat Nu_mid = Nut + dNu_dt * d_t/2 ;
    arma::cx_dmat Nd_mid = Ndt + dNd_dt * d_t/2 ;

    derivative(currentSimulationTime + d_t/2, Nu_mid, Nd_mid, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * d_t ;
    Ndt = Ndt + dNd_dt * d_t ;
}


/// @brief Time step of a trajectory, computed with "original" 4th order Runge-Kutta method.
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK4(double d_t)
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
    arma::cx_dmat Nu_k1 = Nut + dNu_dt_k1 * d_t/2 ;
    arma::cx_dmat Nd_k1 = Ndt + dNd_dt_k1 * d_t/2 ;
    derivative(currentSimulationTime + d_t/2, Nu_k1 , Nd_k1, dNu_dt_k2, dNd_dt_k2);
// k3 : 0.5 |  0   0.5
    arma::cx_dmat Nu_k2 = Nut + dNu_dt_k2 * d_t/2 ;
    arma::cx_dmat Nd_k2 = Ndt + dNd_dt_k2 * d_t/2 ;
    derivative(currentSimulationTime + d_t/2, Nu_k2 , Nd_k2, dNu_dt_k3, dNd_dt_k3);
// k3 : 1   |  0    0    1
    arma::cx_dmat Nu_k3 = Nut + dNu_dt_k3 * d_t ;
    arma::cx_dmat Nd_k3 = Ndt + dNd_dt_k3 * d_t ;
    derivative(currentSimulationTime + d_t, Nu_k3 , Nd_k3, dNu_dt_k4, dNd_dt_k4);
// k4 :     |  1/6  2/6  2/6  1/6  
    Nut = Nut + (dNu_dt_k1 + 2* dNu_dt_k2 + 2*dNu_dt_k3 + dNu_dt_k4) * d_t/6;
    Ndt = Ndt + (dNd_dt_k1 + 2* dNd_dt_k2 + 2*dNd_dt_k3 + dNd_dt_k4) * d_t/6;
}


/// @brief Time step of a trajectory, computed with the "3/8 rule" 4th order Runge-Kutta method.
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK4_38(double d_t)
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
    arma::cx_dmat Nu_k = Nut + dNu_dt_k1 * d_t/3 ;
    arma::cx_dmat Nd_k = Ndt + dNd_dt_k1 * d_t/3 ;
    derivative(currentSimulationTime + d_t/3, Nu_k , Nd_k, dNu_dt_k2, dNd_dt_k2);
// k3 : 2/3  |  -1/3  1
    Nu_k = Nut + (- dNu_dt_k1/3 + dNu_dt_k2) * d_t ;
    Nd_k = Ndt + (- dNd_dt_k1/3 + dNd_dt_k2) * d_t ;
    derivative(currentSimulationTime + 2*dt/3, Nu_k , Nd_k, dNu_dt_k3, dNd_dt_k3);
// k3 :  1   |  1     -1    1
    Nu_k = Nut + (dNu_dt_k1 - dNu_dt_k2 + dNu_dt_k3) * d_t;
    Nd_k = Ndt + (dNd_dt_k1 - dNd_dt_k2 + dNd_dt_k3) * d_t;
    derivative(currentSimulationTime + d_t, Nu_k , Nd_k, dNu_dt_k4, dNd_dt_k4);
//    :      | 1/8    3/8   3/8   1
    Nut = Nut + (dNu_dt_k1 + dNu_dt_k2*3 + dNu_dt_k3*3 + dNu_dt_k4) * d_t/8;
    Ndt = Ndt + (dNd_dt_k1 + dNd_dt_k2*3 + dNd_dt_k3*3 + dNd_dt_k4) * d_t/8;
}


/// @brief ==NOT IMPLEMENTED== Time step of a trajectory, computed with the adaptive 4th/5th order Runge-Kutta method.
/// @param d_t the time step
/// @return void
void FHHFSolver::step_RK45(double d_t)
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

    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * ndt.diag())  + arma::diagmat(*V * nut.diag());
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * nut.diag())  + arma::diagmat(*V * ndt.diag());


    if (!fhSystem->isExternalPotentialNull())
    {

        if (fhSystem->isSystemTimeDependant())
        {
            fhSystem->updateExternalPotential(currentTime);
        }
        DeltaUp += *extPot;
        DeltaDo += *extPot;
    }



    dnut_dt = - jj * ( nut * DeltaUp - DeltaUp * nut) ;
    dndt_dt = - jj * ( ndt * DeltaDo - DeltaDo * ndt) ;
}


/// @brief Launch simulation and write the results
void FHHFSolver::runSimulation( )
{


    std::cout << "   -- Start simulation." << std::endl;


    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);

    Nut = Nu0;
    Ndt = Nd0;


    arma::cx_dcube Nu;
    arma::cx_dcube Nd;
    arma::cx_dmat Niiu;
    arma::cx_dmat Niid;
    arma::cx_dmat Mkku;
    arma::cx_dmat Mkkd;
    arma::cx_drowvec EU;
    arma::cx_drowvec EtH;
    arma::cx_drowvec EpotExt;
    arma::cx_drowvec E;
    arma::drowvec nPartU_t =  arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec nPartD_t =  arma::drowvec(nTimePoints, arma::fill::zeros);
    nPartU_t(0) = std::real(arma::trace(Nut));
    nPartD_t(0) = std::real(arma::trace(Ndt));

    if (saveWholeDMatrix)
    {
        Nu = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        Nd = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        Nu.slice(0) = Nut;
        Nd.slice(0) = Ndt;
    }
    if (saveDiagDMatrix)
    {
        Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        Niiu.col(0) = Nut.diag();
        Niid.col(0) = Ndt.diag();
    }
    if (saveMomenta)
    {
        Mkku = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        Mkkd = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        if (computeMomenta)
        {
            arma::cx_dmat V_1p = arma::cx_dmat(*fhSystem->get_tbgs_V_Pointer());
            arma::cx_dmat Mu = V_1p.t() * Nut * V_1p;
            arma::cx_dmat Md = V_1p.t() * Ndt * V_1p;
            Mkku.col(0) = Mu.diag();
            Mkkd.col(0) = Md.diag();
        }
    }
    if (saveEnergies)
    {
        EU = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
        EtH = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
        EpotExt = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
        EtH(0) =  arma::accu(*tM % Nut) + arma::accu(*tM % Ndt);
        EU(0) = (arma::cdot(Nut.diag(), *U * Ndt.diag()) + arma::cdot(Ndt.diag(),  *U * Nut.diag())) / arma::cx_double(2,0);
        if (fhSystem->getExternalPotentialType() != FHExternalPotentialType::EPnull)
        {
            EpotExt(0) =  arma::accu(*extPot % Nut) + arma::accu(*extPot % Ndt);
        }
    }


    arma::cx_dmat dNu_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


// Simulation :
    for (unsigned int t = 1;  t< nTimePoints; t++)
    {
        
        for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
        {
            switch (integrator)
            {
            case HFIntegrator::RungeKutta1 :
                step_RK1(dt/nStepsPerTimePoint);
                break;
            case HFIntegrator::RungeKutta1_Mem :
                step_RK1_Mem(dt/nStepsPerTimePoint, dNu_dt, dNd_dt);
                break;
            case HFIntegrator::RungeKutta2 :
                step_RK2(dt/nStepsPerTimePoint);
                break;
            case HFIntegrator::RungeKutta4 :
                step_RK4(dt/nStepsPerTimePoint);
                break;
            case HFIntegrator::RungeKutta4_38 :
                step_RK4_38(dt/nStepsPerTimePoint);
                break;
            case HFIntegrator::RungeKutta45 :
                step_RK45(dt/nStepsPerTimePoint);
                break;  
            default:
                break;
            }

            currentSimulationTime += dt/nStepsPerTimePoint;
        }

        // Save Values :
        timeLine[t] = t*dt;
        nPartU_t(t) = std::real(arma::trace(Nut));
        nPartD_t(t) = std::real(arma::trace(Ndt));
        //std::cout << "Save Matrix" << std::endl;
        if (saveWholeDMatrix)
        {
            Nu.slice(t) = Nut;
            Nd.slice(t) = Ndt;
        }
        //std::cout << "Save Diag" << std::endl;
        if (saveDiagDMatrix)
        {
            Niiu.col(t) = Nut.diag();
            Niid.col(t) = Ndt.diag();
        }
        if (saveDiag2DMatrix)
        {
        }
        if (saveMomenta)
        {
            if (computeMomenta)
            {
                arma::cx_dmat V_1p = arma::cx_dmat(*fhSystem->get_tbgs_V_Pointer());
                arma::cx_dmat Mu = V_1p.t() * Nut * V_1p;
                arma::cx_dmat Md = V_1p.t() * Ndt * V_1p;
                Mkku.col(t) = Mu.diag();
                Mkkd.col(t) = Md.diag();
            }
        }
        //std::cout << "Save En, t " << t << std::endl;
        if (saveEnergies)
        {
            EtH(t) =  arma::accu(*tM % Nut) + arma::accu(*tM % Ndt);
            EU(t) = (arma::cdot(Nut.diag(), *U * Ndt.diag()) + arma::cdot(Ndt.diag(),  *U * Nut.diag())) / arma::cx_double(2,0);
            if (fhSystem->getExternalPotentialType() != FHExternalPotentialType::EPnull)
            {
                EpotExt(t) =  arma::accu(*extPot % Nut) + arma::accu(*extPot % Ndt);
            }

        }
        
    }


// Computing other values :
//std::cout << "Save En tot" << std::endl;
    if (saveEnergies)
    {
       E = EU + EtH + EpotExt;
    }

    if (talkative)
    {
        printf("--- Mean field dynamic, n1u : \n");
        (arma::real(Niiu.row(0))).print();
        printf("--- Mean field dynamic, E : \n");
        (arma::real(E)).print();
    }


       // Writing the results :
    if (writeResultsTxt)
    {
        ResultWriterTxt results = ResultWriterTxt(workName, folderName);
        results.newResults();

        results.writeLine("t", timeLine);
        if (saveDiagDMatrix)
        {
            for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
            {
                results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Niiu.row(iSite));
            }
            for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
            {
                results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_d", Niid.row(iSite));
            }
        }
        if (saveMomenta)
        {
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
        }
        if (saveEnergies)
        {
            results.writeLineInReal("EU", EU);
            results.writeLineInReal("EtH", EtH);
            results.writeLineInReal("Epe", EpotExt);
            results.writeLineInReal("E", E);
        }


        std::cout << "   -- Results written in " << folderName << workName << ".txt" << std::endl;
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

        if (saveDiagDMatrix)
        {
        
            Json::Value niiuJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < Niiu.n_rows ; iSite++)
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
        if (saveEnergies)
        {
            Json::Value EJson(Json::objectValue);
            EJson["U"] = rJson.writeLineInReal(EU);
            EJson["tH"] = rJson.writeLineInReal(EtH);
            EJson["pe"] = rJson.writeLineInReal(EpotExt);
            EJson["tot"] = rJson.writeLineInReal(E);
            simJson["E"] = EJson;
        }

        resJson["simulation"] = simJson;


        Json::StreamWriterBuilder wbuilder;
        wbuilder["indentation"] = "\t";
        std::string document = Json::writeString(wbuilder, resJson);
        std::string fileAdressName = folderName  + workName + ".json";
        std::ofstream resultFile(fileAdressName);
        resultFile << document;
        resultFile.close();

        std::cout << "   -- Results written in " << folderName   << workName << ".json" << std::endl;
    }

    if (writeResultsHDF5)
    {
        
    }



    

}




/// @brief Destructor
FHHFSolver::~FHHFSolver()
{

}
 