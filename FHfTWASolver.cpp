#include "FHfTWASolver.h"


FHfTWASolver::FHfTWASolver(FHSystem* fhS)
{
    fhSystem = fhS;

    unsigned int nSites = fhSystem->getnSites();
    arma::dmat zeroMat = arma::dmat(nSites, nSites, arma::fill::zeros);
    arma::dmat haldDiagReal = arma::dmat(nSites, nSites, arma::fill::eye)/2;
    halfDiag = arma::cx_dmat(haldDiagReal, zeroMat);

    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
}

void FHfTWASolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
    writeResults = true;
}

/// @brief Set time simulation parameters.
/// @param t_max the maximum simulation time
/// @param d_t the time step between stored values in the simulation 
/// @param nSPTP number of simulation steps between two time steps
/// @return void
void FHfTWASolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
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

    if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::givenDiagonal)
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
    else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::quenched)
    {

        nPartU = fhSystem->getUpSpinParticleNumber();
        nPartD = fhSystem->getDownSpinParticleNumber();
        
        arma::dvec eigVal;
        arma::dmat Jreal = arma::real(*tM);
        arma::eig_sym(eigVal, V_1p, Jreal);
        arma::dmat V_1p_1 = V_1p.t();

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
    printf("Test 6 - Variance matrix : \n");
    sqrtNoiseVarianceMatrixU.print();

}



/// @brief Launch simulation and write the results
void FHfTWASolver::runSimulation( )
{

    printf("Test 7 - fTWA dynamic, nTraj=%i \n", nTraj);
    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    arma::cx_dcube rhou_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

    arma::cx_drowvec EU   = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EtH  = arma::cx_drowvec(nTimePoints, arma::fill::zeros);

   // arma::cx_dmat rhou_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    arma::cx_dmat halfDiagMatrix = arma::cx_dmat(0.5 * arma::dmat(nSites, nSites, arma::fill::eye), arma::dmat(nSites, nSites, arma::fill::zeros)); //<- Matrix of diagnoal 1/2

// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel shared(rhou0_mean, rhod0_mean)
{
    // Local (in threads) versions of important variables :
    arma::cx_dcube rhou_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

    arma::cx_drowvec EU_loc   = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EtH_loc  = arma::cx_drowvec(nTimePoints, arma::fill::zeros);

    #pragma omp for 
    for (unsigned int i_traj=0; i_traj <nTraj; i_traj ++)
    {
        arma::cx_dmat Nu_loc;
        arma::cx_dmat Nd_loc;
        arma::cx_dmat rhou;
        arma::cx_dmat rhod;

        // White Noise Matrices :
        arma::cx_dmat wnu = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
        arma::cx_dmat wnd = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

        arma::arma_rng::set_seed_random();
        arma::cx_dmat wnM = arma::cx_dmat(nSites, nSites, arma::fill::randn)/sqrt(2);

        for (unsigned int i = 0; i<nSites; i++)
        {
            for (unsigned int j = i+1; j<nSites; j++)
            {               
                wnu(i,j) += wnM(i,j) * sqrtNoiseVarianceMatrixU(i,j);
                wnu(j,i) += std::conj(wnM(i,j)) * sqrtNoiseVarianceMatrixU(j,i);
                wnd(i,j) += wnM(j,i) * sqrtNoiseVarianceMatrixD(i,j);
                wnd(j,i) += std::conj(wnM(j,i)) * sqrtNoiseVarianceMatrixD(j,i);
            }
        }



        if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::givenDiagonal)
        {
           // Initial rho matrices :
           rhou = rhou0_mean + wnu;  
           rhod = rhod0_mean + wnd;

        }
        else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::quenched)
        {
            // Initial matrices in momentum space :
            arma::cx_dmat  mU = arma::diagmat(diagMu0) - halfDiag + wnu;
            arma::cx_dmat  mD = arma::diagmat(diagMd0) - halfDiag + wnd;

            // Initial matrices in position space :
            rhou = V_1p * mU * V_1p_1 ;
            rhod = V_1p * mD * V_1p_1;

        }

        // We compute the trajectory dynamic and add up the diagonal to the result
        for (unsigned int t = 0;  t<nTimePoints; t++)
        {   
            // We first store the density in the thread-local variables  
            timeLine[t] = t*dt;                                     //<- time line, stored too many times. To Fix ?
            rhou_t_loc.slice(t) = rhou_t_loc.slice(t) + rhou;       //<- Up-spin density matrix
            rhod_t_loc.slice(t) = rhod_t_loc.slice(t) + rhod;       //<- Down-spin density matrix
            rho2iiu_loc.col(t) = rho2iiu_loc.col(t) + arma::conj(rhou.diag()) % rhou.diag();    //<- Diagonal of the squared norm of up-spin density matrix 
            rho2iid_loc.col(t) = rho2iid_loc.col(t) + arma::conj(rhod.diag()) % rhod.diag();    //<- Diagonal of the squared norm of down-spin density matrix

            Nu_loc =  rhou + halfDiagMatrix;
            Nd_loc =  rhod + halfDiagMatrix;
            EU_loc(t)=  EU_loc(t) + (arma::dot(Nu_loc.diag(), *U * Nd_loc.diag()) + arma::dot(Nd_loc.diag(), *U * Nu_loc.diag())) / arma::cx_double(2,0);
            EtH_loc(t) = EtH_loc(t) + arma::accu((Nu_loc + Nd_loc) % *tM);

            // We do 'nStepsPerTimePoint' simulation steps
            for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
            {
                if (integrator==fTWAIntegrator::RungeKutta2) step_RK2(rhou, rhod, dt/nStepsPerTimePoint);                   //<- RK2 integrator chosen
                else if (integrator==fTWAIntegrator::RungeKutta4_38) step_RK4_38(rhou, rhod, dt/nStepsPerTimePoint);        //<- RK4 integrator chosen
            }
        }


        
unsigned int tt = std::floor(nTraj/10);
    if ((i_traj%tt) == 0)
    {
        printf("Trajectory %i computed by thread %i \n", i_traj, omp_get_thread_num());    
    }

    }

    // We add up all thread-local variables into the global ones :
    #pragma omp critical
    {
        rhou_t +=rhou_t_loc;
        rhod_t +=rhod_t_loc;
        rho2iiu +=rho2iiu_loc;
        rho2iid +=rho2iid_loc;
        EU += EU_loc;
        EtH += EtH_loc;
    }


}

// of 
     arma::cx_double nTraj_arma = arma::cx_double(nTraj,0);

    arma::cx_dcube Nu       = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube Nd       = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat N2iiu     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat N2iid     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

    EU = EU/nTraj_arma;
    EtH = EtH/nTraj_arma;
    

    for (unsigned int t = 0; t < nTimePoints; t++)
    {
        Nu.slice(t) = rhou_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhou_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
        Nd.slice(t) = rhod_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhod_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
        N2iiu.col(t) = rho2iiu.col(t) / nTraj_arma + ((rhou_t.slice(t)).diag() + arma::conj((rhou_t.slice(t)).diag()))/ arma::cx_double(2*nTraj,0) + arma::cx_dmat(nSites, 1, arma::fill::ones) / arma::cx_double(4,0);
        N2iid.col(t) = rho2iid.col(t) / nTraj_arma + ((rhod_t.slice(t)).diag() + arma::conj((rhod_t.slice(t)).diag()))/ arma::cx_double(2*nTraj,0) + arma::cx_dmat(nSites, 1, arma::fill::ones) / arma::cx_double(4,0);
    }

    arma::cx_drowvec E = EtH+ EU;

    // Writing the results :
    if (writeResults)
    {
        ResultWriter results = ResultWriter(workName, folderName);
        results.newResults();
        results.writeLine("t", timeLine);
        for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
        {
            results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Nu.tube(iSite, iSite));
        }
        results.writeLineInReal("EU", EU);
        results.writeLineInReal("EtH", EtH);
        results.writeLineInReal("E", E);
    }
}





/// @brief Time derivative of the trajectories
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param drhou_dt the time derivative of the up-spin density matrix, to be updated
/// @param drhod_dt the time derivative of the down-spin density matrix, to be updated 
/// @return void
void FHfTWASolver::derivative(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{

    arma::cx_dvec half = arma::cx_double(0.5,0) *  arma::cx_dvec(rhou.n_rows, arma::fill::ones); 

    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * (rhod.diag() + half));
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * (rhou.diag() + half));

    drhou_dt = - jj * ( rhou * DeltaUp - DeltaUp * rhou) ;
    drhod_dt = - jj * ( rhod * DeltaDo - DeltaDo * rhod) ;
}


/// @brief Time step of a trajectory, computed with the 2nd order Runge-Kutta method. 
/// @param rhou the up-spin density matrix, to be updated
/// @param rhod the down-spin density matrix , to be updated
/// @param dt the time step
/// @return void
void FHfTWASolver::step_RK2(arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt)
{
    arma::cx_dmat drhou_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(rhou, rhod, drhou_dt, drhod_dt);

    arma::cx_dmat rhou_mid = rhou + drhou_dt * dt/2 ;
    arma::cx_dmat rhod_mid = rhod + drhod_dt * dt/2 ;

    derivative(rhou_mid, rhod_mid, drhou_dt, drhod_dt);

    rhou = rhou + drhou_dt * dt ;
    rhod = rhod + drhod_dt * dt ;
}

/// @brief Time step of a trajectory, computed with the "3/8 rule" 4th order Runge-Kutta method.
/// @param rhou the up-spin density matrix, to be updated
/// @param rhod the down-spin density matrix , to be updated
/// @param dt the time step
/// @return void
void FHfTWASolver::step_RK4_38(arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt)
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
    derivative(rhou, rhod, drhou_dt_k1, drhod_dt_k1);
// k2 : 1/3   |   1/3
    arma::cx_dmat rhou_k = rhou + drhou_dt_k1 * dt/3 ;
    arma::cx_dmat rhod_k = rhod + drhod_dt_k1 * dt/3 ;
    derivative(rhou_k, rhod_k, drhou_dt_k2, drhod_dt_k2);
// k3 : 2/3   |   -1/3    1
    rhou_k = rhou + (- drhou_dt_k1 /3 + drhou_dt_k2) * dt ;
    rhod_k = rhod + (- drhod_dt_k1 /3 + drhod_dt_k2) * dt ;
    derivative(rhou_k, rhod_k, drhou_dt_k3, drhod_dt_k3);
// k4 : 1     |    1      -1     1
    rhou_k = rhou + (drhou_dt_k1 - drhou_dt_k2 + drhou_dt_k3) * dt ;
    rhod_k = rhod + (drhod_dt_k1 - drhod_dt_k2 + drhod_dt_k3) * dt ;
    derivative(rhou_k, rhod_k, drhou_dt_k4, drhod_dt_k4);

    rhou = rhou + (drhou_dt_k1 + drhou_dt_k2*3 + drhou_dt_k3*3 + drhou_dt_k4) * dt / 8;
    rhod = rhod + (drhod_dt_k1 + drhod_dt_k2*3 + drhod_dt_k3*3 + drhod_dt_k4) * dt / 8;
}



FHfTWASolver::~FHfTWASolver()
{
}
