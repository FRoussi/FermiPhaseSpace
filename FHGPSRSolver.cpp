#include "FHGPSRSolver.h"


FHGPSRSolver::FHGPSRSolver(FHSystem* fhS)
{
    fhSystem = fhS;


    U = fhSystem->getUMatrixPointer();
    J = fhSystem->getJMatrixPointer();

    sqrtU =  arma::sqrt(*U);
}

void FHGPSRSolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
}

void FHGPSRSolver::setGPSRParameters(unsigned int nT)
{
    this->nTraj = nT;
}

void FHGPSRSolver::initiateSolver()
{

    rhou0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    rhod0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    for (unsigned int it= 0; it < fhSystem->getnSites(); it ++)
    {
           rhou0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_upSpin()[it],   0);
           rhod0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_downSpin()[it], 0);
    }

}

void FHGPSRSolver::step_Drift(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Drift, arma::cx_dmat* drhod_Drift)
{
    arma::cx_dmat drhou_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(rhou, rhod, &drhou_dt, &drhod_dt);

    arma::cx_dmat rhou_mid = *rhou + drhou_dt * dt/2 ;
    arma::cx_dmat rhod_mid = *rhod + drhod_dt * dt/2 ;

    derivative(&rhou_mid, &rhod_mid, &drhou_dt, &drhou_dt);

    *drhou_Drift = drhou_dt * dt;
    *drhod_Drift = drhod_dt * dt;
}

void FHGPSRSolver::step_Diffusion(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt, arma::cx_dmat * drhou_Diff, arma::cx_dmat* drhod_Diff)
{
    unsigned int nSites = fhSystem->getnSites();

    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - *rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - *rhod;

    arma::cx_dvec diagW1 = arma::cx_dvec(nSites, arma::fill::randn) * std::sqrt(dt);
    arma::cx_dvec diagW2 = arma::cx_dvec(nSites, arma::fill::randn) * std::sqrt(dt);


    arma::cx_dmat W1 = arma::diagmat(diagW1);
    arma::cx_dmat W2 = arma::diagmat(diagW2);


    *drhou_Diff = sqrt_mjj * *rhou * sqrtU * W1 * rhou_H                + sqrt_pjj * rhou_H * W2 * *rhou ;
    *drhod_Diff = sqrt_mjj * *rhod * sqrtU * arma::conj(W1) * rhod_H    + sqrt_pjj * rhod_H * arma::conj(W2) * *rhod ;
}

void FHGPSRSolver::step(arma::cx_dmat * rhou, arma::cx_dmat* rhod, double dt)
{
    arma::cx_dmat drhou_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat drhou_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    step_Drift(rhou, rhod, dt, &drhou_Drift, &drhod_Drift);
    step_Diffusion(rhou, rhod, dt, &drhou_Diff, &drhod_Diff);

    *rhou = *rhou + drhou_Drift + drhou_Diff;
    *rhod = *rhod + drhod_Drift + drhod_Diff;

}

void FHGPSRSolver::runSimulation( )
{

    printf("GPSR dynamic, nTraj=%i \n", nTraj);
    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    arma::cx_dcube rhou_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);





// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel
{
    arma::cx_dcube rhou_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dcube rhod_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);

    arma::cx_dmat rho2iiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat rho2iid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

    #pragma omp for 
    for (unsigned int i_traj=0; i_traj <nTraj; i_traj ++)
    {
        // Noiseless initial rho matrices :
        arma::cx_dmat rhou = rhou0;
        arma::cx_dmat rhod = rhod0;

      
        // We compute the trajectory dynamic and add up the diagnoal to the result
        for (unsigned int t = 0;  t<nTimePoints; t++)
        {
            timeLine[t] = t*dt;
            rhou_t_loc.slice(t) = rhou_t_loc.slice(t) + rhou;
            rhod_t_loc.slice(t) = rhod_t_loc.slice(t) + rhod;
            rho2iiu_loc.col(t) = rho2iiu_loc.col(t) + arma::conj(rhou.diag()) % rhou.diag();
            rho2iid_loc.col(t) = rho2iid_loc.col(t) + arma::conj(rhod.diag()) % rhod.diag();
            step(&rhou, &rhod, dt);
            for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s++)
            {
                
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
    

    
    arma::cx_dcube Nu       = rhou_t / arma::cx_double(nTraj,0) ; 
    arma::cx_dcube Nd       = rhod_t / arma::cx_double(nTraj,0) ; 
    arma::cx_dmat N2iiu     = rho2iiu / arma::cx_double(nTraj,0) ;
    arma::cx_dmat N2iid     = rho2iid / arma::cx_double(nTraj,0) ;
    arma::cx_drowvec EU     = arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EJ     = arma::cx_drowvec(nTimePoints, arma::fill::zeros);

    for (unsigned int t = 0; t < nTimePoints; t++)
    {
        EU(t)= (arma::cdot(Nu.slice(t).diag(), *U * Nd.slice(t).diag()) + arma::cdot(Nd.slice(t).diag(), *U * Nu.slice(t).diag())) / arma::cx_double(2,0);
        EJ(t) = arma::accu((Nu.slice(t) + Nd.slice(t)) % *J);
    }

    arma::cx_drowvec E     = EJ + EU;


    // Writing the results :
    std::string folderName = fhSystem->getFolder();
    std::string fileName = fhSystem->getName();
    std::string fHFName = fileName.substr(0, fileName.size() - 4);
    std::string fHFAdressName = fHFName + "_GPSR";

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

void FHGPSRSolver::derivative(arma::cx_dmat* rhou, arma::cx_dmat* rhod, arma::cx_dmat* drhou_dt, arma::cx_dmat* drhod_dt)
{
    arma::cx_dmat DeltaUp = *J + arma::diagmat(*U * rhod->diag());
    arma::cx_dmat DeltaDo = *J + arma::diagmat(*U * rhou->diag());

    *drhou_dt = - jj * ( *rhou * DeltaUp - DeltaUp * *rhou) ;
    *drhod_dt = - jj * ( *rhod * DeltaDo - DeltaDo * *rhod) ;
}


FHGPSRSolver::~FHGPSRSolver()
{
}
