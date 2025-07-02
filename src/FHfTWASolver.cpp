#include "../include/FHfTWASolver.h"


FHfTWASolver::FHfTWASolver(FHSystem* fhS) : FHSolver(fhS)
{

 
    
    nTraj = 1;
    integrator = fTWAIntegrator::RungeKutta2;

    unsigned int nSites = fhSystem->getnSites();
    arma::dmat zeroMat = arma::dmat(nSites, nSites, arma::fill::zeros);

    halfDiag = 0.5 * arma::cx_dmat(nSites, nSites, arma::fill::eye);
    halfVec = 0.5 * arma::cx_dvec(nSites, arma::fill::ones);



    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
    //J = fhSystem->getJMatrixPointer();
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

    std::cout << "   -- fTWA Solver initiation." << std::endl;

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
            diagNu0.at(it) = 1 * fhSystem->getInitialCondition_upSpin()[it];
            diagNd0.at(it) = 1 * fhSystem->getInitialCondition_downSpin()[it];

            rhou0_mean.at(it, it) = arma::cx_double(diagNu0(it) - 0.5, 0);
            rhod0_mean.at(it, it) = arma::cx_double(diagNd0(it) - 0.5, 0);          
        }

        nPartU = static_cast<unsigned int>(std::floor(arma::accu(diagNu0)));
        nPartD = static_cast<unsigned int>(std::floor(arma::accu(diagNd0)));

        sqrtNoiseVarianceMatrixU = arma::dmat(nSites, nSites, arma::fill::zeros);
        sqrtNoiseVarianceMatrixD = arma::dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int i = 0; i < nSites; i++)
        {
            for (unsigned int j = i+1; j < nSites; j++)
            {
                sqrtNoiseVarianceMatrixU.at(i,j) = std::sqrt((diagNu0[i] + diagNu0[j] - 2 * diagNu0[i] * diagNu0[j])/2);
                sqrtNoiseVarianceMatrixU.at(j,i) = sqrtNoiseVarianceMatrixU.at(i,j);
                sqrtNoiseVarianceMatrixD.at(i,j) = std::sqrt((diagNd0[i] + diagNd0[j] - 2 * diagNd0[i] * diagNd0[j])/2);
                sqrtNoiseVarianceMatrixD.at(j,i) = sqrtNoiseVarianceMatrixD.at(i,j);
            }
        }
    }
    else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::tbGroundState)
    {

        if (!fhSystem->isTBGroundStateComputed())
        {
            fhSystem->computeTBGroundState();
        }

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();
        

        diagMu0 = arma::dvec(nSites, arma::fill::zeros);
        diagMd0 = arma::dvec(nSites, arma::fill::zeros);

        diagMu0.subvec(0,nPartU-1).ones();
        diagMd0.subvec(0,nPartD-1).ones();

        sqrtNoiseVarianceMatrixU = arma::dmat(nSites, nSites, arma::fill::zeros);
        sqrtNoiseVarianceMatrixD = arma::dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int i = 0; i < nSites; i++)
        {
            for (unsigned int j = i+1; j < nSites; j++)
            {
                sqrtNoiseVarianceMatrixU.at(i,j) = std::sqrt((diagMu0[i] + diagMu0[j] - 2 * diagMu0[i] * diagMu0[j])/2);
                sqrtNoiseVarianceMatrixU.at(j,i) = sqrtNoiseVarianceMatrixU.at(i,j);
                sqrtNoiseVarianceMatrixD.at(i,j) = std::sqrt((diagMd0[i] + diagMd0[j] - 2 * diagMd0[i] * diagMd0[j])/2);
                sqrtNoiseVarianceMatrixD.at(j,i) = sqrtNoiseVarianceMatrixD.at(i,j);
            }
        }
    }
    else if  (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::mfGroundState)
    {

        if (!fhSystem->isMFGroundStateComputed())
        {
            fhSystem->computeMFGroundState();
        }

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();
        

        diagMu0 = arma::dvec(nSites, arma::fill::zeros);
        diagMd0 = arma::dvec(nSites, arma::fill::zeros);

        diagMu0.subvec(0,nPartU-1).ones();
        diagMd0.subvec(0,nPartD-1).ones();

        sqrtNoiseVarianceMatrixU = arma::dmat(nSites, nSites, arma::fill::zeros);
        sqrtNoiseVarianceMatrixD = arma::dmat(nSites, nSites, arma::fill::zeros);
        for (unsigned int i = 0; i < nSites; i++)
        {
            for (unsigned int j = i+1; j < nSites; j++)
            {
                sqrtNoiseVarianceMatrixU.at(i,j) = std::sqrt((diagMu0[i] + diagMu0[j] - 2 * diagMu0[i] * diagMu0[j])/2);
                sqrtNoiseVarianceMatrixU.at(j,i) = sqrtNoiseVarianceMatrixU.at(i,j);
                sqrtNoiseVarianceMatrixD.at(i,j) =  std::sqrt((diagMd0[i] + diagMd0[j] - 2 * diagMd0[i] * diagMd0[j])/2);
                sqrtNoiseVarianceMatrixD.at(j,i) = sqrtNoiseVarianceMatrixD.at(i,j);
            }
        }
    }


    if (saveTBStates)
    {
        if (!fhSystem->isTBGroundStateComputed())
        {
            fhSystem->computeTBGroundState();
        }
    }
    if (saveMFStates)
    {
        if (!fhSystem->isMFGroundStateComputed())
        {
            fhSystem->computeMFGroundState();
        }
    }

    //sqrtNoiseVarianceMatrixU.print();
    
    //sqrtNoiseVarianceMatrixU = arma::sqrt(sqrtNoiseVarianceMatrixU);
    //sqrtNoiseVarianceMatrixD = arma::sqrt(sqrtNoiseVarianceMatrixD);

    if (talkative)
    {
        std::cout << "--- fTWA Variance matrix : " << std::endl ;
        sqrtNoiseVarianceMatrixU.print();
    }

}







/// @brief Launch simulation and write the results
void FHfTWASolver::runSimulation( )
{

    std::cout << "   -- Start simulation." << std::endl;
    //printf("Test 7 - fTWA dynamic, nTraj=%i \n", nTraj);


    unsigned int nSites = fhSystem->getnSites();
    
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    //std::cout << "nSites=" << nSites << " , nTPoints="<< nTimePoints <<  std::endl;
    arma::cx_dcube  rhou_t,     rhod_t;
    arma::cx_dmat   rhoiiu,     rhoiid;
    arma::cx_dmat   rho2iiu,    rho2iid;
    arma::cx_dmat   rhoiurhoid;
    arma::cx_dmat   Energies;
    arma::cx_dmat   mkku, mkkd, mukku, mukkd;

    if (saveWholeDMatrix)
    {
        rhou_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu      = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid      = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiag2DMatrix)
    {
        rho2iid     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rho2iiu     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveNiuNid)
    {
        rhoiurhoid  = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        Energies    = arma::cx_dmat(nEnergies + 1,nTimePoints, arma::fill::zeros);
    }
    if (saveTBStates)
    {
        mkku        = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mkkd        = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveMFStates)
    {
        mukku       = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mukkd       = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }

   // arma::cx_dmat rhou_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    //arma::cx_dmat halfDiagMatrix = arma::cx_dmat(0.5 * arma::dmat(nSites, nSites, arma::fill::eye), arma::dmat(nSites, nSites, arma::fill::zeros)); //<- Matrix of diagnoal 1/2





// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel shared(rhou0_mean, rhod0_mean)
{
    // Local (in threads) versions of important variables :

    arma::cx_dcube rhou_t_loc, rhod_t_loc;
    arma::cx_dmat rhoiiu_loc, rhoiid_loc; 
    arma::cx_dmat rho2iiu_loc, rho2iid_loc; 
    arma::cx_dmat rhoiurhoid_loc;

    arma::cx_dmat En_loc;   
    arma::cx_dmat mkku_loc, mkkd_loc, mukku_loc, mukkd_loc;


    if (saveWholeDMatrix)
    {
        rhou_t_loc      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t_loc      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu_loc      = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid_loc      = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiag2DMatrix)
    {
        rho2iiu_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rho2iid_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveNiuNid)
    {
        rhoiurhoid_loc  = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveTBStates)
    {
        mkku_loc        = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mkkd_loc        = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveMFStates)
    {
        mukku_loc       = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mukkd_loc       = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        En_loc          = arma::cx_dmat(nEnergies,nTimePoints, arma::fill::zeros);
    }

    arma::cx_dmat Nu_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat Nd_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhou      = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhod      = arma::cx_dmat(nSites, nSites);

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

            // Initial matrices in position space 
            rhou = (*fhSystem->get_tbgs_V_Pointer()) * mU * ((*fhSystem->get_tbgs_V_Pointer()).t()) ;
            rhod = (*fhSystem->get_tbgs_V_Pointer()) * mD * ((*fhSystem->get_tbgs_V_Pointer()).t()) ;
    
        }
        else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::mfGroundState)
        {
            // Initial matrices in momentum space :
            arma::cx_dmat  mU = arma::diagmat(diagMu0) - halfDiag + wnu;
            arma::cx_dmat  mD = arma::diagmat(diagMd0) - halfDiag + wnd;

            // Initial matrices in position space :
            rhou = (*fhSystem->get_mfgs_Vu_Pointer()) * mU * ((*fhSystem->get_mfgs_Vu_Pointer()).t()) ;
            rhod = (*fhSystem->get_mfgs_Vd_Pointer()) * mD * ((*fhSystem->get_mfgs_Vd_Pointer()).t()) ;

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
                rho2iiu_loc.col(t)      += arma::conj(rhou.diag()) % rhou.diag();    //<- Diagonal of the squared norm of up-spin density matrix 
                rho2iid_loc.col(t)      += arma::conj(rhod.diag()) % rhod.diag();    //<- Diagonal of the squared norm of down-spin density matrix
            }
            if (saveTBStates)
            {
                mkku_loc.col(t)         += arma::diagvec((*fhSystem->get_tbgs_V_Pointer()).t() * rhou *(*fhSystem->get_tbgs_V_Pointer()));
                mkkd_loc.col(t)         += arma::diagvec((*fhSystem->get_tbgs_V_Pointer()).t() * rhod *(*fhSystem->get_tbgs_V_Pointer()));
            }
            if (saveMFStates)
            {
                mukku_loc.col(t)        += arma::diagvec((*fhSystem->get_mfgs_Vu_Pointer()).t() * rhou *(*fhSystem->get_mfgs_Vu_Pointer()));
                mukkd_loc.col(t)        += arma::diagvec((*fhSystem->get_mfgs_Vd_Pointer()).t() * rhod *(*fhSystem->get_mfgs_Vd_Pointer()));
            }
            if (saveNiuNid)
            {
                rhoiurhoid_loc.col(t)   += rhou.diag() % rhod.diag();
            }
            if (saveEnergies)
            {
                Nu_loc =  rhou + halfDiag;
                Nd_loc =  rhod + halfDiag;
                En_loc.at(0,t) += (arma::dot(Nu_loc.diag(), *U * Nd_loc.diag()) + arma::dot(Nd_loc.diag(), *U * Nu_loc.diag())) / arma::cx_double(2,0);
                En_loc.at(1,t) += arma::accu((Nu_loc + Nd_loc) % (*tM));
                En_loc.at(2,t) += (arma::dot(Nu_loc.diag(), *V * Nu_loc.diag()) + arma::dot(Nd_loc.diag(), *V * Nd_loc.diag())) / arma::cx_double(2,0);
                if (!fhSystem->isExternalPotentialNull())
                {
                    En_loc.at(3,t) +=  arma::accu((Nu_loc + Nd_loc) % (*extPot));
                }
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
        if (saveNiuNid)
        {
            rhoiurhoid += rhoiurhoid_loc;
        }
        if (saveTBStates)
        {
            mkku += mkku_loc;
            mkkd += mkkd_loc;
        }
        if (saveMFStates)
        {
            mukku += mukku_loc;
            mukkd += mukkd_loc;
        }
        if (saveEnergies)
        {
            for (uint ie =0; ie < nEnergies; ie ++)
            {
                Energies.row(ie) +=En_loc.row(ie);
            }
        }

    }


}

    //std::cout << "Computation of observables" << std::endl;
// Computation of different observables :
    arma::cx_double nTraj_arma = arma::cx_double(nTraj,0);

    arma::cx_dcube  Nu, Nd   ;
    arma::cx_dmat   Niiu, Niid ;   
    arma::cx_dmat   N2iiu, N2iid;   
    arma::cx_dmat   Mkku ,  Mkkd ;
    arma::cx_dmat   Mukku,  Mukkd;
    arma::cx_dmat   NiuNid;


    if (saveWholeDMatrix)
    {
        //std::cout << "Computation of N" << std::endl;
        Nu = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        Nd = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            Nu.slice(t) = rhou_t.slice(t) / nTraj_arma + halfDiag;
            Nd.slice(t) = rhod_t.slice(t) / nTraj_arma + halfDiag;
        }

    }
    if (saveDiagDMatrix)
    {
        //std::cout << "Computation of Nii" << std::endl;
        Niiu = rhoiiu / nTraj_arma + arma::cx_dmat(arma::size(rhoiiu), arma::fill::ones)/ arma::cx_double(2,0);
        Niid = rhoiid / nTraj_arma + arma::cx_dmat(arma::size(rhoiid), arma::fill::ones)/ arma::cx_double(2,0);

    }
    if (saveDiag2DMatrix)
    {
        N2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        N2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        arma::cx_dcolvec vec_4 = arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            N2iiu.col(t) = rho2iiu.col(t) / nTraj_arma + (rhoiiu.col(t) + arma::conj(rhoiiu.col(t)))/ arma::cx_double(2*nTraj,0) + vec_4;
            N2iid.col(t) = rho2iid.col(t) / nTraj_arma + (rhoiid.col(t) + arma::conj(rhoiid.col(t)))/ arma::cx_double(2*nTraj,0) + vec_4;
        }
    }
    if (saveNiuNid)
    {
        NiuNid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        arma::cx_dcolvec vec_4 = arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            NiuNid.col(t) = rhoiurhoid.col(t) / nTraj_arma + (rhoiiu.col(t) + rhoiid.col(t))/ arma::cx_double(2*nTraj,0) + vec_4;
        }
    }
    if (saveTBStates)
    {
        Mkku = mkku / nTraj_arma + arma::cx_dmat(arma::size(mkku), arma::fill::ones)/ arma::cx_double(2,0);
        Mkkd = mkkd / nTraj_arma + arma::cx_dmat(arma::size(mkkd), arma::fill::ones)/ arma::cx_double(2,0);
    }
    if (saveMFStates)
    {
        Mukku = mukku / nTraj_arma + arma::cx_dmat(arma::size(mukku), arma::fill::ones)/ arma::cx_double(2,0);
        Mukkd = mukkd / nTraj_arma + arma::cx_dmat(arma::size(mukkd), arma::fill::ones)/ arma::cx_double(2,0);
    }
    if (saveEnergies)
    {
        Energies = Energies/nTraj_arma;
        Energies.row(nEnergies) = arma::sum(Energies,0);
    }

    writeResults(Nu, Nd, Niiu, Niid, Mkku, Mkkd, Mukku, Mukkd, N2iiu, N2iid, NiuNid, Energies, timeLine);

}





/// @brief Time derivative of the trajectories
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param drhou_dt the time derivative of the up-spin density matrix, to be updated
/// @param drhod_dt the time derivative of the down-spin density matrix, to be updated 
/// @return void
void FHfTWASolver::derivative(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{

    arma::cx_dmat DeltaUp = *tM; 
    arma::cx_dmat DeltaDo = *tM;
    
    if (fhSystem->getCoulombInteractionType() != FHCoulombInteractionType::CEnull)
    {
        DeltaUp += arma::diagmat(*U * (rhod.diag() + halfVec)) ;
        DeltaDo += arma::diagmat(*U * (rhou.diag() + halfVec)) ;
        if (fhSystem->getCoulombInteractionType() != FHCoulombInteractionType::CEonSite)
        {
            DeltaUp += arma::diagmat(*V * (rhou.diag() + halfVec)) ;
            DeltaDo += arma::diagmat(*V * (rhod.diag() + halfVec)) ;       
        }
    }

    
    if (!fhSystem->isExternalPotentialNull())
    {
        if (fhSystem->isSystemTimeDependant())
        {
           fhSystem->updateExternalPotential(currentSimulationTime);
        }

        DeltaUp += *extPot;
        DeltaDo += *extPot;

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
