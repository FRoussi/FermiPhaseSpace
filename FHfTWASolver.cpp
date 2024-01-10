#include "FHfTWASolver.h"


FHfTWASolver::FHfTWASolver(FHSystem* fhS)
{
    fhSystem = fhS;

    unsigned int nSites = fhSystem->getnSites();
    arma::dmat zeroMat = arma::dmat(nSites, nSites, arma::fill::zeros);
    arma::dmat haldDiagReal = arma::dmat(nSites, nSites, arma::fill::eye)/2;
    halfDiag = arma::cx_dmat(haldDiagReal, zeroMat);

    U = fhSystem->getUMatrixPointer();
    J = fhSystem->getJMatrixPointer();
}

void FHfTWASolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
}

void FHfTWASolver::setfTWAParameters(unsigned int nT)
{
    this->nTraj = nT;
}

void FHfTWASolver::initiateSolver()
{

    unsigned int nSites = fhSystem->getnSites(); 

    rhou0_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    rhod0_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    diagNu0 = arma::dvec(nSites, arma::fill::zeros);
    diagNd0 = arma::dvec(nSites, arma::fill::zeros);
    for (unsigned int it= 0; it < nSites; it ++)
    {
        diagNu0(it) = 1 * fhSystem->getInitialCondition_upSpin()[it];
        diagNd0(it) = 1 * fhSystem->getInitialCondition_downSpin()[it];

        rhou0_mean(it, it) = arma::cx_double(diagNu0(it) - 0.5, 0);
        rhod0_mean(it, it) = arma::cx_double(diagNd0(it) - 0.5, 0);          
    }


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
    
    sqrtNoiseVarianceMatrixU = arma::sqrt(sqrtNoiseVarianceMatrixU);
    sqrtNoiseVarianceMatrixD = arma::sqrt(sqrtNoiseVarianceMatrixD);
    printf("Test 6 - Variance matrix : \n");
    sqrtNoiseVarianceMatrixU.print();

}

void FHfTWASolver::step_RK2(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt)
{
    arma::cx_dmat drhou_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(rhou, rhod, &drhou_dt, &drhod_dt);

    arma::cx_dmat rhou_mid = *rhou + drhou_dt * dt/2 ;
    arma::cx_dmat rhod_mid = *rhod + drhod_dt * dt/2 ;

    derivative(&rhou_mid, &rhod_mid, &drhou_dt, &drhod_dt);

    *rhou = *rhou + drhou_dt * dt ;
    *rhod = *rhod + drhod_dt * dt ;
}

void FHfTWASolver::step_RK4_38(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt)
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
    derivative(rhou, rhod, &drhou_dt_k1, &drhod_dt_k1);
// k2 : 1/3   |   1/3
    arma::cx_dmat rhou_k = *rhou + drhou_dt_k1 * dt/3 ;
    arma::cx_dmat rhod_k = *rhod + drhod_dt_k1 * dt/3 ;
    derivative(&rhou_k, &rhod_k, &drhou_dt_k2, &drhod_dt_k2);
// k3 : 2/3   |   -1/3    1
    rhou_k = *rhou + (- drhou_dt_k1 /3 + drhou_dt_k2) * dt ;
    rhod_k = *rhod + (- drhod_dt_k1 /3 + drhod_dt_k2) * dt ;
    derivative(&rhou_k, &rhod_k, &drhou_dt_k3, &drhod_dt_k3);
// k4 : 1     |    1      -1     1
    rhou_k = *rhou + (drhou_dt_k1 - drhou_dt_k2 + drhou_dt_k3) * dt ;
    rhod_k = *rhod + (drhod_dt_k1 - drhod_dt_k2 + drhod_dt_k3) * dt ;
    derivative(&rhou_k, &rhod_k, &drhou_dt_k4, &drhod_dt_k4);

    *rhou = *rhou + (drhou_dt_k1 + drhou_dt_k2*3 + drhou_dt_k3*3 + drhou_dt_k4) * dt / 8;
    *rhod = *rhod + (drhod_dt_k1 + drhod_dt_k2*3 + drhod_dt_k3*3 + drhod_dt_k4) * dt / 8;
}

void FHfTWASolver::derivative(const arma::cx_dmat* rhou, const arma::cx_dmat* rhod, arma::cx_dmat* drhou_dt, arma::cx_dmat* drhod_dt)
{

    arma::cx_dvec half = arma::cx_double(0.5,0) *  arma::cx_dvec(rhou->n_rows, arma::fill::ones); 

    arma::cx_dmat DeltaUp = *J + arma::diagmat(*U * (rhod->diag() + half));
    arma::cx_dmat DeltaDo = *J + arma::diagmat(*U * (rhou->diag() + half));

    *drhou_dt = - jj * ( *rhou * DeltaUp - DeltaUp * *rhou) ;
    *drhod_dt = - jj * ( *rhod * DeltaDo - DeltaDo * *rhod) ;
}


void FHfTWASolver::runSimulation( )
{

    printf("Test 7 - fTWA dynamic, nTraj=%i \n", nTraj);
    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    arma::cx_dcube rhou_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

   // arma::cx_dmat rhou_mean = arma::cx_dmat(nSites, nSites, arma::fill::zeros);




// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel shared(rhou0_mean, rhod0_mean)
{
    // Local (in threads) versions of important variables :
    arma::cx_dcube rhou_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

    #pragma omp loop 
    for (unsigned int i_traj=0; i_traj <nTraj; i_traj ++)
    {
        // Noiseless initial rho matrices :
        arma::cx_dmat rhou = rhou0_mean;
        arma::cx_dmat rhod = rhod0_mean;

        // We add up the white noises : 
        arma::arma_rng::set_seed_random();
        arma::cx_dmat wnM = arma::cx_dmat(nSites, nSites, arma::fill::randn);
        for (unsigned int i = 0; i<nSites; i++)
        {
            for (unsigned int j = i+1; j<nSites; j++)
            {
                rhou(i,j) += wnM(i,j) * sqrtNoiseVarianceMatrixU(i,j);
                rhou(j,i) += std::conj(wnM(i,j)) * sqrtNoiseVarianceMatrixU(j,i);
                rhod(i,j) += wnM(j,i) * sqrtNoiseVarianceMatrixD(i,j);
                rhod(j,i) += std::conj(wnM(j,i)) * sqrtNoiseVarianceMatrixD(j,i);
            }
        }

        // We compute the trajectory dynamic and add up the diagnoal to the result
        for (unsigned int t = 0;  t<nTimePoints; t++)
        {
            timeLine[t] = t*dt;          
            rhou_t_loc.slice(t) = rhou_t_loc.slice(t) + rhou;
            rhod_t_loc.slice(t) = rhod_t_loc.slice(t) + rhod;
            rho2iiu_loc.col(t) = rho2iiu_loc.col(t) + arma::conj(rhou.diag()) % rhou.diag();
            rho2iid_loc.col(t) = rho2iid_loc.col(t) + arma::conj(rhod.diag()) % rhod.diag();
            for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
            {
                step_RK4_38(&rhou, &rhod, dt/nStepsPerTimePoint);
            }
        }


        printf("Trajectory %i computed by thread %i \n", i_traj, omp_get_thread_num());       
    }

    #pragma omp critical
    {
        rhou_t +=rhou_t_loc;
        rhod_t +=rhod_t_loc;
        rho2iiu +=rho2iiu_loc;
        rho2iid +=rho2iid_loc;
    }


}

// of 

    arma::cx_dcube Nu       = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube Nd       = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat N2iiu     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat N2iid     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EU     = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EJ     = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    
    for (unsigned int t = 0; t < nTimePoints; t++)
    {
        Nu.slice(t) = rhou_t.slice(t) / arma::cx_double(nTraj,0) + arma::cx_dmat(arma::size(rhou_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
        Nd.slice(t) = rhod_t.slice(t) / arma::cx_double(nTraj,0) + arma::cx_dmat(arma::size(rhod_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
        N2iiu.col(t) = rho2iiu.col(t) / arma::cx_double(nTraj,0) + ((rhou_t.slice(t)).diag() + arma::conj((rhou_t.slice(t)).diag()))/ arma::cx_double(2*nTraj,0) + arma::cx_dmat(arma::size(rho2iiu.col(t)), arma::fill::ones) / arma::cx_double(4,0);
        N2iid.col(t) = rho2iid.col(t) / arma::cx_double(nTraj,0) + ((rhod_t.slice(t)).diag() + arma::conj((rhod_t.slice(t)).diag()))/ arma::cx_double(2*nTraj,0) + arma::cx_dmat(arma::size(rho2iid.col(t)), arma::fill::ones) / arma::cx_double(4,0);
        EU(t)= (arma::cdot(Nu.slice(t).diag(), *U * Nd.slice(t).diag()) + arma::cdot(Nd.slice(t).diag(), *U * Nu.slice(t).diag())) / arma::cx_double(2,0);
        EJ(t) = arma::accu((Nu.slice(t) + Nd.slice(t)) % *J);
    }


    arma::cx_drowvec E = EJ+ EU;

    // Writing the results :
    std::string folderName = fhSystem->getFolder();
    std::string fileName = fhSystem->getName();
    std::string fHFName = fileName.substr(0, fileName.size() - 4);
    std::string fHFAdressName = fHFName + "_fTWA";

    ResultWriter results = ResultWriter(fHFAdressName, folderName);
    results.newResults();
    results.writeLine("t", timeLine);
    for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
    {
        results.writeLineInReal("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", Nu.tube(iSite, iSite));
    }
    results.writeLineInReal("E_U", EU);
    results.writeLineInReal("E_J", EJ);
    results.writeLineInReal("E", E);

}

FHfTWASolver::~FHfTWASolver()
{
}
