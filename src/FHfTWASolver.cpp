#include "../include/FHfTWASolver.h"


FHfTWASolver::FHfTWASolver(FHSystem* fhS) : FHSolver(fhS)
{

 
    
    nTraj = 1;
    integrator = fTWAIntegrator::RungeKutta2;

    unsigned int nSites = fhSystem->getnSites();
    arma::dmat zeroMat = arma::dmat(nSites, nSites, arma::fill::zeros);

    halfDiag = 0.5 * arma::cx_dmat(nSites, nSites, arma::fill::eye);
    halfVec = 0.5 * arma::cx_dvec(nSites, arma::fill::eye);

    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
    J = fhSystem->getJMatrixPointer();
    V = fhSystem->getVMatrixPointer();
    extPot = fhSystem->getExtPotPointer();


}


void                        FHfTWASolver::setIntegratorOrder(int order)
{
    switch (order)
    {
        case 2:
        {
            integratorOrder = 2;
            integrator = fTWAIntegrator::RungeKutta2 ;
            break;
        }
        case 4:
        {
            integratorOrder = 4;
            integrator = fTWAIntegrator::RungeKutta4_38 ;
            break;
        }
        default:
        {
            std::cout << "Integrator order " << order << " not implemented. Switching to default 2." << std::endl;
            integratorOrder = 2;
            integrator = fTWAIntegrator::RungeKutta2 ;
            break;
        }
    }

}


void                        FHfTWASolver::setMethodParameters(Json::Value param)
{
    double param_dt = param.get("dt", 0.01).asDouble();
    double param_tm = param.get("t_max", 2).asDouble();
    int param_sptp = param.get("stepsParTimePoints", 1).asInt();
    int param_nTraj = param.get("n_traj", 10).asInt();
    int param_intOrder = param.get("int_order", 2).asInt();

    setTimeParameters(param_tm, param_dt, param_sptp);
    setIntegratorOrder(param_intOrder);
    setfTWAParameters(param_nTraj);
}


/// @brief Set fTWA simulation parameters.
/// @param nT number of trajectories used
/// @return void
void FHfTWASolver::setfTWAParameters(unsigned int nT)
{
    this->nTraj = nT;
}



/// @brief Initiate the fTWA solver, compute the initial density matrix probability and the noise variance to apply.
void FHfTWASolver::initiateSolver()
{

    integrator = fTWAIntegrator::RungeKutta4_38;
    unsigned int nSites = fhSystem->getnSites(); 

    rhou0_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    rhod0_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    diagNu0 = arma::dvec(nSites, arma::fill::zeros);
    diagNd0 = arma::dvec(nSites, arma::fill::zeros);

    if ((fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::givenDiagonal) || (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::spinWave) || (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::chargeWave))
    {

        for (unsigned int it= 0; it < nSites; it ++)
        {
            diagNu0(it) = 1 * fhSystem->getInitialCondition_upSpin()[it];
            diagNd0(it) = 1 * fhSystem->getInitialCondition_downSpin()[it];

            rhou0_mean(it, it) = arma::cx_double(diagNu0(it) - 0.5, 0);
            rhod0_mean(it, it) = arma::cx_double(diagNd0(it) - 0.5, 0);          
        }

        nPartU = static_cast<unsigned int>(std::floor(arma::accu(diagNu0)));
        nPartD = static_cast<unsigned int>(std::floor(arma::accu(diagNd0)));

        sqrtNoiseVarianceMatrixU = arma::dmat(nSites, nSites, arma::fill::zeros);
        sqrtNoiseVarianceMatrixD = arma::dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int i = 0; i < nSites; i++)
        {
            for (unsigned int j = i+1; j < nSites; j++)
            {
                sqrtNoiseVarianceMatrixU(i,j) = (diagNu0(i) + diagNu0(j) - 2 * diagNu0(i) * diagNu0(j))/2;
                sqrtNoiseVarianceMatrixU(j,i) = sqrtNoiseVarianceMatrixU(i,j);
                sqrtNoiseVarianceMatrixD(i,j) = (diagNd0(i) + diagNd0(j) - 2 * diagNd0(i) * diagNd0(j))/2;
                sqrtNoiseVarianceMatrixD(j,i) = sqrtNoiseVarianceMatrixD(i,j);
            }
        }
    }
    else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::tbGroundState)
    {

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();
        
        setComputeMomenta(true);

        diagMu0 = arma::dvec(nSites, arma::fill::zeros);
        diagMd0 = arma::dvec(nSites, arma::fill::zeros);

            for (unsigned int k=0; k<nPartU ; k++)   {
                diagMu0(k) = 1; 
            }
            for (unsigned int k=0; k<nPartD ; k++) {
                diagMd0(k) = 1; 
            }

        sqrtNoiseVarianceMatrixU = arma::dmat(nSites, nSites, arma::fill::zeros);
        sqrtNoiseVarianceMatrixD = arma::dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int i = 0; i < nSites; i++)
        {
            for (unsigned int j = i+1; j < nSites; j++)
            {
                sqrtNoiseVarianceMatrixU(i,j) = (diagMu0(i) + diagMu0(j) - 2 * diagMu0(i) * diagMu0(j))/2;
                sqrtNoiseVarianceMatrixU(j,i) = sqrtNoiseVarianceMatrixU(i,j);
                sqrtNoiseVarianceMatrixD(i,j) = (diagMd0(i) + diagMd0(j) - 2 * diagMd0(i) * diagMd0(j))/2;
                sqrtNoiseVarianceMatrixD(j,i) = sqrtNoiseVarianceMatrixD(i,j);
            }
        }
    }
    
    sqrtNoiseVarianceMatrixU = arma::sqrt(sqrtNoiseVarianceMatrixU);
    sqrtNoiseVarianceMatrixD = arma::sqrt(sqrtNoiseVarianceMatrixD);

    if (talkative)
    {
        std::cout << "--- fTWA Variance matrix : " << std::endl ;
        sqrtNoiseVarianceMatrixU.print();
    }


}


void FHfTWASolver::setComputeMomenta(bool yn)
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





/// @brief Launch simulation and write the results
void FHfTWASolver::runSimulation( )
{

    std::cout << "   -- Start simulation" << std::endl;
    //printf("Test 7 - fTWA dynamic, nTraj=%i \n", nTraj);


    unsigned int nSites = fhSystem->getnSites();
    
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    //std::cout << "nSites=" << nSites << " , nTPoints="<< nTimePoints <<  std::endl;
    arma::cx_dcube rhou_t; 
    arma::cx_dcube rhod_t;
    arma::cx_dmat rhoiiu;
    arma::cx_dmat rhoiid;
    arma::cx_dmat rho2iiu;
    arma::cx_dmat rho2iid;
    arma::cx_drowvec EU;  
    arma::cx_drowvec EtH; 

    if (saveWholeDMatrix)
    {
        rhou_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu    = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid    = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiag2DMatrix)
    {
        rho2iid     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rho2iiu     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        EU          = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
        EtH         = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    }


   // arma::cx_dmat rhou_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    arma::cx_dmat halfDiagMatrix = arma::cx_dmat(0.5 * arma::dmat(nSites, nSites, arma::fill::eye), arma::dmat(nSites, nSites, arma::fill::zeros)); //<- Matrix of diagnoal 1/2





// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel shared(rhou0_mean, rhod0_mean)
{
    // Local (in threads) versions of important variables :

    arma::cx_dcube rhou_t_loc;
    arma::cx_dcube rhod_t_loc; 
    arma::cx_dmat rhoiiu_loc; 
    arma::cx_dmat rhoiid_loc; 
    arma::cx_dmat rho2iiu_loc; 
    arma::cx_dmat rho2iid_loc; 
    arma::cx_drowvec EU_loc;   
    arma::cx_drowvec EtH_loc;  


    if (saveWholeDMatrix)
    {
        rhou_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiag2DMatrix)
    {
        rho2iiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rho2iid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        EU_loc   = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
        EtH_loc  = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    }


    arma::cx_dmat Nu_loc = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat Nd_loc = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhou = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhod = arma::cx_dmat(nSites, nSites);




    #pragma omp for 
    for (unsigned int i_traj=0; i_traj <nTraj; i_traj ++)
    {

        double currentSimulationTime = 0;


        // Complex white noise matrices for initial conditions :
        arma::arma_rng::set_seed_random();
        arma::cx_dmat  wnM = arma::cx_dmat(nSites, nSites, arma::fill::randn)/sqrt(2);  //<- Complex white noise
        arma::cx_dmat  wnM_up = arma::symmatu(wnM);                                     //<- Hermitian complex white noise
        arma::cx_dmat wnM_down = arma::symmatl(wnM);                                    //<- Hermitian complex white noise

        arma::cx_dmat wnu = wnM_up % sqrtNoiseVarianceMatrixU;                          //<- Noise to add to 
        arma::cx_dmat wnd = wnM_down % sqrtNoiseVarianceMatrixD;

        // We miss memeory space so we pay attention to useless matrices :
        wnM.reset();
        wnM_up.reset();
        wnM_up.reset();





        if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::tbGroundState)
        {
            // Initial matrices in momentum space :
            arma::cx_dmat  mU = arma::diagmat(diagMu0) - halfDiag + wnu;
            arma::cx_dmat  mD = arma::diagmat(diagMd0) - halfDiag + wnd;

            // Initial matrices in position space :
            rhou = V_1p * mU * V_1p_1 ;
            rhod = V_1p * mD * V_1p_1;

        }
        else 
        {
           // Initial rho matrices :
           rhou = rhou0_mean + wnu;  
           rhod = rhod0_mean + wnd;

        }


        // We compute the trajectory dynamic and add up the diagonal to the result
        for (unsigned int t = 0;  t<nTimePoints; t++)
        {   
            // We first store the density in the thread-local variables  
            timeLine[t] = t*dt;                                     //<- time line, stored too many times. To Fix ?
            if (saveWholeDMatrix)
            {
                rhou_t_loc.slice(t) +=  rhou;       //<- Up-spin density matrix
                rhod_t_loc.slice(t) +=  rhod;       //<- Down-spin density matrix
            }
            if (saveDiagDMatrix)
            {
                rhoiiu_loc.col(t) += rhou.diag();
                rhoiid_loc.col(t) += rhod.diag();
            }
            if (saveDiag2DMatrix)
            {
                rho2iiu_loc.col(t) += arma::conj(rhou.diag()) % rhou.diag();    //<- Diagonal of the squared norm of up-spin density matrix 
                rho2iid_loc.col(t) += arma::conj(rhod.diag()) % rhod.diag();    //<- Diagonal of the squared norm of down-spin density matrix
            }
            if (saveEnergies)
            {
                Nu_loc =  rhou + halfDiagMatrix;
                Nd_loc =  rhod + halfDiagMatrix;
                EU_loc(t)=  EU_loc(t) + (arma::dot(Nu_loc.diag(), *U * Nd_loc.diag()) + arma::dot(Nd_loc.diag(), *U * Nu_loc.diag())) / arma::cx_double(2,0);
                EtH_loc(t) = EtH_loc(t) + arma::accu((Nu_loc + Nd_loc) % *tM);
            }

            // We do 'nStepsPerTimePoint' simulation steps
            for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
            {
                if (integrator==fTWAIntegrator::RungeKutta2) step_RK2(currentSimulationTime, rhou, rhod, dt/nStepsPerTimePoint);                   //<- RK2 integrator chosen
                else if (integrator==fTWAIntegrator::RungeKutta4_38) step_RK4_38(currentSimulationTime, rhou, rhod, dt/nStepsPerTimePoint);        //<- RK4 integrator chosen
                currentSimulationTime += dt/nStepsPerTimePoint;
            }
        }


        // Little idea of where we are in the computation :
        unsigned int tt = std::floor(nTraj/10);
        if ((i_traj%tt) == 0)
        {
            std::cout << "     - Trajectory " << i_traj << " computed by thread " << omp_get_thread_num() << std::endl;    
        }

    }

    // We add up all thread-local variables into the global ones :
   #pragma omp critical
    {
        if (saveWholeDMatrix)
        {
            rhou_t += rhou_t_loc;
            rhod_t += rhod_t_loc;
        }
        if (saveDiagDMatrix)
        {
            rhoiiu += rhoiiu_loc;
            rhoiid += rhoiid_loc;
        }
        if (saveDiag2DMatrix)
        {
            rho2iiu += rho2iiu_loc;
            rho2iid += rho2iid_loc;
        }
        if (saveEnergies)
        {
            EU += EU_loc;
            EtH += EtH_loc;
        }

    }


}

    //std::cout << "Computation of observables" << std::endl;
// Computation of different observables :
    arma::cx_double nTraj_arma = arma::cx_double(nTraj,0);

    arma::cx_dcube Nu;      
    arma::cx_dcube Nd   ;   
    arma::cx_dmat Niiu  ;   
    arma::cx_dmat Niid  ;   
    arma::cx_dmat N2iiu  ;  
    arma::cx_dmat N2iid ;   
    arma::cx_dmat Mkku  ;   
    arma::cx_dmat Mkkd  ;   
    arma::cx_drowvec E;


    if (saveWholeDMatrix)
    {
        //std::cout << "Computation of N" << std::endl;
        Nu = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        Nd = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            Nu.slice(t) = rhou_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhou_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
            Nd.slice(t) = rhod_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhod_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
        }
    }
    if (saveDiagDMatrix)
    {
        //std::cout << "Computation of Nii" << std::endl;
        Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            Niiu.col(t) = rhoiiu.col(t) / nTraj_arma + arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(2,0);
            Niid.col(t) = rhoiid.col(t) / nTraj_arma + arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(2,0);
        }
    }
    if (saveDiag2DMatrix)
    {
        N2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        N2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            N2iiu.col(t) = rho2iiu.col(t) / nTraj_arma + (rhoiiu.col(t) + arma::conj(rhoiiu.col(t)))/ arma::cx_double(2*nTraj,0) + arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
            N2iid.col(t) = rho2iid.col(t) / nTraj_arma + (rhoiid.col(t) + arma::conj(rhoiid.col(t)))/ arma::cx_double(2*nTraj,0) + arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
        }
    }
    if (saveMomenta)
    {
        Mkku = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        Mkkd = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        arma::cx_dmat V_1p_comp = arma::cx_dmat(V_1p, arma::dmat(arma::size(V_1p), arma::fill::zeros)); 
        if (computeMomenta)
        {
            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                Mkku.col(t) = arma::diagvec(V_1p_comp.t() * Nu.slice(t) * V_1p_comp);
                Mkkd.col(t) = arma::diagvec(V_1p_comp.t() * Nd.slice(t) * V_1p_comp);
            }
        }
    }
    if (saveEnergies)
    {
        EU = EU/nTraj_arma;
        EtH = EtH/nTraj_arma;
        E = EtH+ EU;
    }

    

    writeResults(Nu, Nd, Niiu, Niid, Mkku, Mkkd, N2iiu, N2iid, EU, EtH, E, timeLine);


}







void FHfTWASolver:: writeResults(arma::cx_dcube &Nu, arma::cx_dcube &Nd, arma::cx_dmat &Niiu, arma::cx_dmat &Niid, arma::cx_dmat &Mkku, arma::cx_dmat &Mkkd, arma::cx_dmat &N2iiu, arma::cx_dmat &N2iid, arma::cx_drowvec &EU, arma::cx_drowvec &EtH, arma::cx_drowvec &E,std::vector<double> &timeLine)
{

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

        if (computeMomenta)
        {/*
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("M_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_u", Mkku.row(kMom));
            }
            for (unsigned int kMom=0; kMom<fhSystem->getnSites(); kMom++)
            {
                results.writeLineInReal("M_" + std::to_string(kMom) + "_" + std::to_string(kMom) + "_d", Mkkd.row(kMom));
            }*/
        }
        if (saveEnergies)
        {
            results.writeLineInReal("EU", EU);
            results.writeLineInReal("EtH", EtH);
            results.writeLineInReal("E", E);
        }


         std::cout << "== Results written in " << folderName  << workName << ".txt" << std::endl;
    }

    if (writeResultsJson)
    {

        ResultWriterJson rJson; 
        Json::Value resJson = jsonFile;

        resJson["t"] = rJson.writeLine(timeLine);

        if (saveDiagDMatrix)
        {
            Json::Value niiuJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < fhSystem->getnSites() ; iSite++)
            {
                niiuJson.append(rJson.writeLineInReal(Niiu.row(iSite)));
            }
            resJson["Niiu"] = niiuJson;
            Json::Value niidJson(Json::arrayValue);
            for (uint iSite = 0 ; iSite < fhSystem->getnSites() ; iSite++)
            {
                niidJson.append(rJson.writeLineInReal(Niid.row(iSite)));
            }
            resJson["Niid"] = niidJson;
        }

        if (saveEnergies)
        {
            Json::Value EJson(Json::objectValue);
            EJson["U"] = rJson.writeLineInReal(EU);
            EJson["tH"] = rJson.writeLineInReal(EtH);
            EJson["tot"] = rJson.writeLineInReal(E);
            resJson["E"] = EJson;
        }



        Json::StreamWriterBuilder wbuilder;
        wbuilder["indentation"] = "\t";
        std::string document = Json::writeString(wbuilder, resJson);
        std::string fileAdressName = folderName + workName + ".json";
        std::ofstream resultFile(fileAdressName);
        resultFile << document;
        resultFile.close();

        std::cout << "== Results written in " << folderName << workName << ".json" << std::endl;
    }

    if (writeResultsHDF5)
    {
        
    }
}


/// @brief Time derivative of the trajectories
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param drhou_dt the time derivative of the up-spin density matrix, to be updated
/// @param drhod_dt the time derivative of the down-spin density matrix, to be updated 
/// @return void
void FHfTWASolver::derivative(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{

    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * (rhod.diag() + halfVec)) ;
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * (rhou.diag() + halfVec)) ;
    
    
    if (!fhSystem->isExternalPotentialNull())
    {
        if (fhSystem->isSystemTimeDependant())
        {
            arma::cx_dmat externalPot = fhSystem->getExternalPotential(currentSimulationTime);
            DeltaUp +=externalPot;
            DeltaDo +=externalPot;
        }
        else
        {
            DeltaUp += *extPot;
            DeltaDo += *extPot;
        }
    }


    drhou_dt = - jj * ( rhou * DeltaUp - DeltaUp * rhou) ;
    drhod_dt = - jj * ( rhod * DeltaDo - DeltaDo * rhod) ;
}


/// @brief Time step of a trajectory, computed with the 2nd order Runge-Kutta method. 
/// @param rhou the up-spin density matrix, to be updated
/// @param rhod the down-spin density matrix , to be updated
/// @param d_t the time step
/// @return void
void FHfTWASolver::step_RK2(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double d_t)
{
    arma::cx_dmat drhou_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(currentSimulationTime, rhou, rhod, drhou_dt, drhod_dt);

    arma::cx_dmat rhou_mid = rhou + drhou_dt * d_t/2 ;
    arma::cx_dmat rhod_mid = rhod + drhod_dt * d_t/2 ;

    derivative(currentSimulationTime + d_t/2, rhou_mid, rhod_mid, drhou_dt, drhod_dt);

    rhou = rhou + drhou_dt * d_t ;
    rhod = rhod + drhod_dt * d_t ;
}

/// @brief Time step of a trajectory, computed with the "3/8 rule" 4th order Runge-Kutta method.
/// @param rhou the up-spin density matrix, to be updated
/// @param rhod the down-spin density matrix , to be updated
/// @param d_t the time step
/// @return void
void FHfTWASolver::step_RK4_38(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double d_t)
{
   unsigned int nSites = fhSystem->getnSites();
    arma::cx_dmat drhou_dt_k1 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt_k1 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    arma::cx_dmat drhou_dt_k2 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt_k2 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    arma::cx_dmat drhou_dt_k3 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt_k3 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    arma::cx_dmat drhou_dt_k4 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt_k4 = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

// k1 :
    derivative(currentSimulationTime, rhou, rhod, drhou_dt_k1, drhod_dt_k1);
// k2 : 1/3   |   1/3
    arma::cx_dmat rhou_k = rhou + drhou_dt_k1 * d_t/3 ;
    arma::cx_dmat rhod_k = rhod + drhod_dt_k1 * d_t/3 ;
    derivative(currentSimulationTime + d_t/3, rhou_k, rhod_k, drhou_dt_k2, drhod_dt_k2);
// k3 : 2/3   |   -1/3    1
    rhou_k = rhou + (- drhou_dt_k1 /3 + drhou_dt_k2) * d_t ;
    rhod_k = rhod + (- drhod_dt_k1 /3 + drhod_dt_k2) * d_t ;
    derivative(currentSimulationTime + 2*dt/3, rhou_k, rhod_k, drhou_dt_k3, drhod_dt_k3);
// k4 : 1     |    1      -1     1
    rhou_k = rhou + (drhou_dt_k1 - drhou_dt_k2 + drhou_dt_k3) * d_t ;
    rhod_k = rhod + (drhod_dt_k1 - drhod_dt_k2 + drhod_dt_k3) * d_t ;
    derivative(currentSimulationTime + d_t, rhou_k, rhod_k, drhou_dt_k4, drhod_dt_k4);

    rhou = rhou + (drhou_dt_k1 + drhou_dt_k2*3 + drhou_dt_k3*3 + drhou_dt_k4) * d_t / 8;
    rhod = rhod + (drhod_dt_k1 + drhod_dt_k2*3 + drhod_dt_k3*3 + drhod_dt_k4) * d_t / 8;
}





FHfTWASolver::~FHfTWASolver()
{
}
