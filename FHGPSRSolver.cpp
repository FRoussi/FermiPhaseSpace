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

void FHGPSRSolver::step_Drift_RK2(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double dt, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift)
{
    unsigned int nSites = fhSystem->getnSites();
    arma::cx_dmat drhou_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    derivative(rhou, rhod, drhou_dt, drhod_dt);

    arma::cx_dmat rhou_mid = rhou + drhou_dt * dt/2 ;
    arma::cx_dmat rhod_mid = rhod + drhod_dt * dt/2 ;

    derivative(rhou_mid, rhod_mid, drhou_dt, drhou_dt);

    drhou_Drift = drhou_dt * dt;
    drhod_Drift = drhod_dt * dt;
}

void FHGPSRSolver::step_Drift_RK4_38(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double dt, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift)
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


    drhou_Drift = (drhou_dt_k1 + drhou_dt_k2*3 + drhou_dt_k3*3 + drhou_dt_k4 )* dt/8;
    drhod_Drift = (drhod_dt_k1 + drhod_dt_k2*3 + drhod_dt_k3*3 + drhod_dt_k4)* dt/8;
}

void FHGPSRSolver::step_Diffusion(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double dt, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff)
{
    unsigned int nSites = fhSystem->getnSites();

    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dvec diagW1 = arma::cx_dvec(nSites, arma::fill::randn) * std::sqrt(dt);
    arma::cx_dvec diagW2 = arma::cx_dvec(nSites, arma::fill::randn) * std::sqrt(dt);

    arma::cx_dmat W1 = arma::diagmat(diagW1);
    arma::cx_dmat W2 = arma::diagmat(diagW2);

    drhou_Diff = sqrt_mjj * rhou * sqrtU * W1 * rhou_H                + sqrt_pjj * rhou_H * W2 * rhou ;
    drhod_Diff = sqrt_mjj * rhod * sqrtU * arma::conj(W1) * rhod_H    + sqrt_pjj * rhod_H * arma::conj(W2) * rhod ;
}

void FHGPSRSolver::step_Diffusion_svd(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double dt, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff)
{

    unsigned int nSites = rhou.n_rows;

    arma::cx_dmat D =diffusionMatrix(rhou, rhod, U0);
    double D_norm = arma::norm(D);

    arma::cx_dmat B =sqrtDecomposition(D);

    arma::cx_dvec W_m = arma::cx_dvec(B.n_cols, arma::fill::randn) * std::sqrt(dt);
    arma::cx_dvec BW = B*W_m;
    arma::cx_dmat BW_u =BW.subvec(0,nSites^2-1);
    arma::cx_dmat BW_d =BW.subvec(nSites^2, 2*nSites^2-1);
    drhou_Diff =  BW_u.reshape(nSites, nSites);
    drhod_Diff =  BW_d.reshape(nSites, nSites);

}


arma::cx_dmat FHGPSRSolver::diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0)
{
   
    unsigned int nSites = rhou.n_rows;
    unsigned int sizeMat = nSites * nSites ;

    std::complex<double> cst = std::sqrt(0.5 * jj * U0 ) ;


    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dmat B11u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B12u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B21u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B22u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B11d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B12d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B21d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B22d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);

    for (unsigned int l = 0 ; l < nSites ; l++)
    {
        for (unsigned int k = 0 ; k < nSites ; k++)
        {
            unsigned int lineIndex = l * nSites + k;
            for (unsigned int p = 0 ; p < nSites ; p++)
            {
                B11u(lineIndex, p) =  (rhou)(p,l)* rhou_H(k,p)  ;
                B12u(lineIndex, p) = jj * (rhou)(p,l)*rhou_H(k,p);
                B21u(lineIndex, p) = jj * rhou_H(p,l)*(rhou)(k,p) ;
                B22u(lineIndex, p) = -  rhou_H(p,l)*(rhou)(k,p)  ;
                        
                B11d(lineIndex, p) =  (rhod)(p,l)*rhod_H(k,p)  ;
                B12d(lineIndex, p) = - jj * (rhod)(p,l)*rhod_H(k,p);
                B21d(lineIndex, p) = jj * rhod_H(p,l)*(rhod)(k,p) ;
                B22d(lineIndex, p) =  rhod_H(p,l)*(rhod)(k,p)  ;
            }
        }

    }
        
    // B = | B11u  B12u  B21u  B22u  |
    //     | B11d  B12d  B21d  B22d  | :
    arma::cx_dmat B = arma::cx_dmat(2 * sizeMat, 4*nSites);
    B.submat(0,         0,          sizeMat-1,      nSites-1)   = B11u;
    B.submat(0,         nSites,     sizeMat-1,      2*nSites-1) = B12u;
    B.submat(0,         2*nSites,   sizeMat-1,      3*nSites-1) = B21u;
    B.submat(0,         3*nSites,   sizeMat-1,      4*nSites-1) = B22u;
    B.submat(sizeMat,   0,          2*sizeMat-1,    nSites-1)   = B11u;
    B.submat(sizeMat,   nSites,     2*sizeMat-1,    2*nSites-1) = B12u;
    B.submat(sizeMat,   2*nSites,   2*sizeMat-1,    3*nSites-1) = B21u;
    B.submat(sizeMat,   3*nSites,   2*sizeMat-1,    4*nSites-1) = B22u;

    B = cst * B;

    return B * B.st();

}

arma::cx_dmat FHGPSRSolver::sqrtDecomposition(arma::cx_dmat& D)
{
    unsigned int sizeD = D.n_rows;

    // D =  | 0     A |
    //      | A^T   0 |
    // D1 = -0.5i A
    arma::cx_dmat D1 = -0.5 * jj * D.submat(0,sizeD/2, sizeD/2-1, sizeD-1);

    arma::cx_dmat U;
    arma::dvec s;
    arma::cx_dmat V;

    // D1 = U s V^T
    arma::svd(U, s, V, D1);

    arma::dmat sqrtS = arma::diagmat(arma::sqrt(s)); 

    // B =  | U.sqrtS           i.U.sqrtS  |
    //      | ii.c(V).sqrtS     c(V).sqrtS |
    arma::cx_dmat B = arma::cx_dmat(sizeD, sizeD, arma::fill::zeros);
    B.submat(0,         0,          sizeD/2-1,    sizeD/2-1)    = U * sqrtS;
    B.submat(0,         sizeD/2,    sizeD/2-1,    sizeD-1)      = FHGPSRSolver::jj * U * sqrtS;
    B.submat(sizeD/2,   0,          sizeD-1,      sizeD/2-1)    = FHGPSRSolver::jj * arma::conj(V) * sqrtS;  
    B.submat(sizeD/2,   sizeD/2,    sizeD-1,      sizeD-1)      = arma::conj(V) * sqrtS; 

    return B;

}

void FHGPSRSolver::step(arma::cx_dmat& rhou, arma::cx_dmat& rhod, double dt)
{
    arma::cx_dmat drhou_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    arma::cx_dmat drhou_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat drhod_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    step_Drift_RK4_38(rhou, rhod, dt, drhou_Drift, drhod_Drift);

    double U0 = std::real( fhSystem->getInteractionEnergy(0,0));
    step_Diffusion_svd(rhou, rhod, U0, dt, drhou_Diff, drhod_Diff);

    rhou = rhou + drhou_Drift + drhou_Diff;
    rhod = rhod + drhod_Drift + drhod_Diff;

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
            step(rhou, rhod, dt);
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

void FHGPSRSolver::derivative(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{
    arma::cx_dmat DeltaUp = *J + arma::diagmat(*U * rhod.diag());
    arma::cx_dmat DeltaDo = *J + arma::diagmat(*U * rhou.diag());

    drhou_dt = - jj * (rhou * DeltaUp - DeltaUp * rhou) ;
    drhod_dt = - jj * (rhod * DeltaDo - DeltaDo * rhod) ;
}


FHGPSRSolver::~FHGPSRSolver()
{
}
