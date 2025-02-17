#include "../include/FHGPSRSolver.h"



FHGPSRSolver::FHGPSRSolver(FHSystem* fhS) : FHSolver(fhS)
{

    nTraj = 1;
    driftIntegrator = GPSRDriftIntegrator::RungeKutta2;
    diffusionIntegrator =  GPSRDiffusionIntegrator::analytical;

    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
   // J = fhSystem->getJMatrixPointer();
    V = fhSystem->getVMatrixPointer();
    extPot = fhSystem->getExtPotPointer();

    sqrtU =  arma::sqrt(*U);




}


void                        FHGPSRSolver::setIntegratorOrder(int order)
{
    switch (order)
    {
        case 2:
        {
            driftIntegratorOrder = 2;
            driftIntegrator = GPSRDriftIntegrator::RungeKutta2 ;
            break;
        }
        case 4:
        {
            driftIntegratorOrder = 4;
            driftIntegrator = GPSRDriftIntegrator::RungeKutta4_38 ;
            break;
        }
        default:
        {
            std::cout << "Integrator order " << order << " not implemented. Switching to default 2." << std::endl;
            driftIntegratorOrder = 2;
            driftIntegrator = GPSRDriftIntegrator::RungeKutta2 ;
            break;
        }
    }

}

void                        FHGPSRSolver::setMethodParameters(Json::Value param)
{
    double param_dt = param.get("dt", 0.01).asDouble();
    double param_tm = param.get("t_max", 2).asDouble();
    int param_sptp = param.get("stepsParTimePoints", 1).asInt();
    int param_nTraj = param.get("n_traj", 10).asInt();
    int param_intOrder = param.get("int_order", 2).asInt();

    setTimeParameters(param_tm, param_dt, param_sptp);
    setIntegratorOrder(param_intOrder);
    setGPSRParameters(param_nTraj);
}


/// @brief Set GPSR simulation parameters.
/// @param nT number of trajectories used
/// @return void
void FHGPSRSolver::setGPSRParameters(unsigned int nT)
{
    this->nTraj = nT;
}



/// @brief Initiate the GPSR solver, compute the initial density matrix probability. 
void FHGPSRSolver::initiateSolver()
{

    std::cout << "   -- GPSR Solver initiation." << std::endl;

    driftIntegrator = GPSRDriftIntegrator::RungeKutta4_38;
    diffusionIntegrator = GPSRDiffusionIntegrator::numerical_rsvd;

    rhou0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    rhod0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    for (unsigned int it= 0; it < fhSystem->getnSites(); it ++)
    {
           rhou0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_upSpin()[it],   0);
           rhod0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_downSpin()[it], 0);
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

}


/// @brief Time derivative of the drift part of the dynamical simulation
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param drhou_dt the time derivative of the up-spin density matrix, to be updated
/// @param drhod_dt the time derivative of the down-spin density matrix, to be updated 
/// @return void
void FHGPSRSolver::derivative(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{
    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * rhod.diag());
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * rhou.diag());

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

    drhou_dt = - jj * (rhou * DeltaUp - DeltaUp * rhou) ;
    drhod_dt = - jj * (rhod * DeltaDo - DeltaDo * rhod) ;
}


/// @brief Drift part a ttime step, computed with Runge-Kutta 2 
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param d_t the time step
/// @param drhou_Drift the up-spin density matrix step, to be updated
/// @param drhod_Drift the down-spin density matrix step, to be updated
/// @return void
void FHGPSRSolver::step_Drift_RK2( double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double d_t, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift)
{
    unsigned int nSites = fhSystem->getnSites();
    arma::cx_dmat drhou_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat drhod_dt = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    derivative(currentSimulationTime, rhou, rhod, drhou_dt, drhod_dt);

    arma::cx_dmat rhou_mid = rhou + drhou_dt * d_t/2 ;
    arma::cx_dmat rhod_mid = rhod + drhod_dt * d_t/2 ;

    derivative(currentSimulationTime + d_t/2, rhou_mid, rhod_mid, drhou_dt, drhou_dt);

    drhou_Drift = drhou_dt * d_t;
    drhod_Drift = drhod_dt * d_t;
}


/// @brief Drift part a ttime step, computed with Runge-Kutta 4 
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param d_t the time step
/// @param drhou_Drift the up-spin density matrix step, to be updated
/// @param drhod_Drift the down-spin density matrix step, to be updated
/// @return void
void FHGPSRSolver::step_Drift_RK4_38(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double d_t, arma::cx_dmat& drhou_Drift, arma::cx_dmat& drhod_Drift)
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
    derivative(currentSimulationTime + 2* d_t/3, rhou_k, rhod_k, drhou_dt_k3, drhod_dt_k3);
// k4 : 1     |    1      -1     1
    rhou_k = rhou + (drhou_dt_k1 - drhou_dt_k2 + drhou_dt_k3) * d_t ;
    rhod_k = rhod + (drhod_dt_k1 - drhod_dt_k2 + drhod_dt_k3) * d_t ;
    derivative(currentSimulationTime + d_t, rhou_k, rhod_k, drhou_dt_k4, drhod_dt_k4);

    drhou_Drift = (drhou_dt_k1 + drhou_dt_k2*3 + drhou_dt_k3*3 + drhou_dt_k4 ) * d_t/8;
    drhod_Drift = (drhod_dt_k1 + drhod_dt_k2*3 + drhod_dt_k3*3 + drhod_dt_k4 ) * d_t/8;

}


/// @brief Diffusion part of a time step, computed with the analytical method 
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param U0 the on-site interaction energy
/// @param d_t the time step
/// @param drhou_Drift the up-spin density matrix step, to be updated
/// @param drhod_Drift the down-spin density matrix step, to be updated
/// @return void
void FHGPSRSolver::step_Diffusion(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double d_t, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff)
{
    unsigned int nSites = fhSystem->getnSites();

    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dvec diagW1 = arma::cx_dvec(nSites, arma::fill::randn)/std::sqrt(2) * std::sqrt(dt);
    arma::cx_dvec diagW2 = arma::cx_dvec(nSites, arma::fill::randn)/std::sqrt(2) * std::sqrt(dt);

    arma::cx_dmat W1 = arma::diagmat(diagW1);
    arma::cx_dmat W2 = arma::diagmat(diagW2);

    std::complex<double> cst = std::sqrt(U0 * jj);

    drhou_Diff = cst * ( rhou_H * W1 * rhou             + jj * rhou * W2 * rhou_H );
    drhod_Diff = cst * ( rhod_H * arma::conj(W1) * rhod + jj * rhod * arma::conj(W2) * rhod_H );
}


/// @brief Diffusion part of a time step, computed with the MaxAndrii method
/// @param rhou the up-spin density matrix 
/// @param rhod the down-spin density matrix 
/// @param U0 the on-site interaction energy 
/// @param d_t the time step
/// @param drhou_Diff the up-spin density matrix step, to be updated
/// @param drhod_Diff the down-spin density matrix step, to be updated
/// @return void
void FHGPSRSolver::step_Diffusion_svd(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0, double d_t, arma::cx_dmat& drhou_Diff, arma::cx_dmat& drhod_Diff)
{

    unsigned int nSites = rhou.n_rows;
    unsigned int sizeMat = std::pow(nSites, 2);

    arma::cx_dmat B = arma::cx_dmat(2*sizeMat, 1, arma::fill::zeros);

    if (fhSystem->getCoulombInteractionType() == FHCoulombInteractionType::CEonSite)
    {
        arma::cx_dmat A =quarterDiffusionMatrix(rhou, rhod, U0);
        double A_norm = arma::norm(A,1);
        if (A_norm > 1.e-7) 
        {
            B = sqrtDecomposition_fromA(A);
        }
    }
    else if (fhSystem->getCoulombInteractionType() == FHCoulombInteractionType::CEgiven)
    {
        arma::cx_dvec Udiag = arma::diagmat(*U);
        arma::cx_dmat D =diffusionMatrix(rhou, rhod, Udiag);
        double D_norm = arma::norm(D,1);
        if (D_norm > 1.e-7) 
        {
            B = sqrtDecomposition(D);
        }
    }
        
    


       

    //printf("Bnorm = %f \n", arma::norm(B,1));

    //std::cout << "Dnorm=" << D_norm << "\t - Decomposition error=" << arma::norm(D - B * B.st())/D_norm << std::endl;


    arma::cx_dvec W_m = arma::cx_dvec( arma::dvec(B.n_cols, arma::fill::randn), arma::dvec(B.n_cols, arma::fill::zeros) ) * std::sqrt(dt); //<- REAL NOISE VECTOR !!!
    arma::cx_dvec BW = B*W_m;
    arma::cx_dmat BW_u =BW.subvec(0,sizeMat-1);
    arma::cx_dmat BW_d =BW.subvec(sizeMat, 2*sizeMat-1);

    //std::cout << "size D=" << arma::size(D) << " - size BW =" << arma::size(BW) << " - size BWu=" << arma::size(BW_u) << " - 2 x nSites^2=" <<  (2*sizeMat) <<  std::endl;

    drhou_Diff =  BW_u.reshape(nSites, nSites);
    drhod_Diff =  BW_d.reshape(nSites, nSites);

}


/// @brief Calculation of a diffusion matrix
/// @param rhou the up-spin density matrix 
/// @param rhod the down-spin density matrix 
/// @param U0 the on-site interaction energy 
/// @return D the diffusion matrix
arma::cx_dmat FHGPSRSolver::diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0)
{
   
    unsigned int nSites = rhou.n_rows;
    unsigned int sizeMat = nSites * nSites ;

    std::complex<double> cst = std::sqrt(0.5 * jj * U0 ) ;

    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dmat B1u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B1d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);

    for (unsigned int l = 0 ; l < nSites ; l++)
    {
        for (unsigned int k = 0 ; k < nSites ; k++)
        {
            unsigned int lineIndex = l * nSites + k;
            for (unsigned int p = 0 ; p < nSites ; p++)
            {
                B1u(lineIndex, p) =  rhou(p,l)* rhou_H(k,p)  ;
                B2u(lineIndex, p) = jj * rhou_H(p,l)*rhou(k,p) ;
                      
                B1d(lineIndex, p) =  rhod(p,l)*rhod_H(k,p)  ;
                B2d(lineIndex, p) = jj * rhod_H(p,l)*rhod(k,p) ;
            }
        }
    }
        

    // B = | B1u  iB1u  B2u   iB2u  |
    //     | B1d -iB1d  B2d  -iB2d  | :
    arma::cx_dmat B = arma::cx_dmat(2 * sizeMat, 4*nSites);
    B.submat(0,         0,          sizeMat-1,      nSites-1)   = B1u;
    B.submat(0,         nSites,     sizeMat-1,      2*nSites-1) = jj * B1u;
    B.submat(0,         2*nSites,   sizeMat-1,      3*nSites-1) = B2u;
    B.submat(0,         3*nSites,   sizeMat-1,      4*nSites-1) = jj*B2u;
    B.submat(sizeMat,   0,          2*sizeMat-1,    nSites-1)   = B1d;
    B.submat(sizeMat,   nSites,     2*sizeMat-1,    2*nSites-1) = - jj * B1d;
    B.submat(sizeMat,   2*nSites,   2*sizeMat-1,    3*nSites-1) = B2d;
    B.submat(sizeMat,   3*nSites,   2*sizeMat-1,    4*nSites-1) = - jj *B2d;

    B = cst * B;

   // std::cout << "size B = " << arma::size(B) << std::endl;

    return B * B.st();

}

arma::cx_dmat       FHGPSRSolver::quarterDiffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, double U0)
{
      unsigned int nSites = rhou.n_rows;
    unsigned int sizeMat = nSites * nSites ;

    std::complex<double> cst = jj * U0  ;

    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dmat B1u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B1d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);

    for (unsigned int l = 0 ; l < nSites ; l++)
    {
        for (unsigned int k = 0 ; k < nSites ; k++)
        {
            unsigned int lineIndex = l * nSites + k;
            for (unsigned int p = 0 ; p < nSites ; p++)
            {
                B1u(lineIndex, p) =  rhou(p,l)* rhou_H(k,p)  ;
                B2u(lineIndex, p) = jj * rhou_H(p,l)*rhou(k,p) ;
                      
                B1d(lineIndex, p) =  rhod(p,l)*rhod_H(k,p)  ;
                B2d(lineIndex, p) = jj * rhod_H(p,l)*rhod(k,p) ;
            }
        }
    }
        

   // // Bu = | B1u  iB1u  B2u   iB2u  | :
   // arma::cx_dmat Bu = arma::cx_dmat(sizeMat, 4*nSites);
   // Bu.submat(0,   0,          sizeMat-1,  nSites-1)   = cst * B1u;
   // Bu.submat(0,   nSites,     sizeMat-1,  2*nSites-1) = cst * jj * B1u;
   // Bu.submat(0,   2*nSites,   sizeMat-1,  3*nSites-1) = cst * B2u;
   // Bu.submat(0,   3*nSites,   sizeMat-1,  4*nSites-1) = cst * jj*B2u;
   // // Bd = | B1d -iB1d  B2d  -iB2d  | :
   // arma::cx_dmat Bd = arma::cx_dmat(sizeMat, 4*nSites);
   // Bd.submat(0,   0,          sizeMat-1,  nSites-1)   = cst * B1d;
   // Bd.submat(0,   nSites,     sizeMat-1,  2*nSites-1) = - cst * jj * B1d;
   // Bd.submat(0,   2*nSites,   sizeMat-1,  3*nSites-1) = cst * B2d;
   // Bd.submat(0,   3*nSites,   sizeMat-1,  4*nSites-1) = - cst * jj *B2d;



   // std::cout << "size B = " << arma::size(B) << std::endl;

    return (cst*(B1u * B1d.st() + B2u * B2d.st()));

}


arma::cx_dmat       FHGPSRSolver::diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dvec Udiag)
{
    uint nSites = rhou.n_cols;
    uint sizeMat = nSites * nSites ;

    std::complex<double> cst = std::sqrt(0.5*jj)  ;


    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dmat B1u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B1d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);

    for (unsigned int l = 0 ; l < nSites ; l++)
    {
        for (unsigned int k = 0 ; k < nSites ; k++)
        {
            unsigned int lineIndex = l * nSites + k;
            for (unsigned int p = 0 ; p < nSites ; p++)
            {
                B1u(lineIndex, p) =  rhou(p,l)* rhou_H(k,p)  ;
                B2u(lineIndex, p) = jj * rhou_H(p,l)*rhou(k,p) ;
                      
                B1d(lineIndex, p) =  rhod(p,l)*rhod_H(k,p)  ;
                B2d(lineIndex, p) = jj * rhod_H(p,l)*rhod(k,p) ;
            }
        }
    }

        // B = | B1u  iB1u  B2u   iB2u  |
    //     | B1d -iB1d  B2d  -iB2d  | :
    arma::cx_dmat B = arma::cx_dmat(2 * sizeMat, 4*nSites);
    B.submat(0,         0,          sizeMat-1,      nSites-1)   = B1u;
    B.submat(0,         nSites,     sizeMat-1,      2*nSites-1) = jj * B1u;
    B.submat(0,         2*nSites,   sizeMat-1,      3*nSites-1) = B2u;
    B.submat(0,         3*nSites,   sizeMat-1,      4*nSites-1) = jj*B2u;
    B.submat(sizeMat,   0,          2*sizeMat-1,    nSites-1)   = B1d;
    B.submat(sizeMat,   nSites,     2*sizeMat-1,    2*nSites-1) = - jj * B1d;
    B.submat(sizeMat,   2*nSites,   2*sizeMat-1,    3*nSites-1) = B2d;
    B.submat(sizeMat,   3*nSites,   2*sizeMat-1,    4*nSites-1) = - jj *B2d;

    B = cst * B;

    
    arma::cx_dvec Udiag4 = arma::cx_dvec(4*nSites);
    Udiag4.subvec(0       ,   nSites-1) = Udiag;        // Udiag4 = | Ud    0    0    0   |
    Udiag4.subvec(nSites  , 2*nSites-1) = Udiag;        //          |  0   Ud    0    0   |  
    Udiag4.subvec(2*nSites, 3*nSites-1) = Udiag;        //          |  0    0   Ud    0   |
    Udiag4.subvec(3*nSites, 4*nSites-1) = Udiag;        //          |  0    0    0   Ud   |

    arma::cx_dmat Ua = arma::diagmat(Udiag4);


    return B * Udiag4 * B.st();
}

arma::cx_dmat       FHGPSRSolver::diffusionMatrix(const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat UMat)
{   
    uint nSites = rhou.n_cols;
    uint sizeMat = nSites * nSites ;

    std::complex<double> cst = std::sqrt(0.5*jj)  ;


    arma::cx_dmat rhou_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhou;
    arma::cx_dmat rhod_H = arma::cx_dmat(nSites, nSites, arma::fill::eye) - rhod;

    arma::cx_dmat B1u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2u = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B1d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);
    arma::cx_dmat B2d = arma::cx_dmat(sizeMat, nSites, arma::fill::zeros);

    for (unsigned int l = 0 ; l < nSites ; l++)
    {
        for (unsigned int k = 0 ; k < nSites ; k++)
        {
            unsigned int lineIndex = l * nSites + k;
            for (unsigned int p = 0 ; p < nSites ; p++)
            {
                B1u(lineIndex, p) =  rhou(p,l)* rhou_H(k,p)  ;
                B2u(lineIndex, p) = jj * rhou_H(p,l)*rhou(k,p) ;
                      
                B1d(lineIndex, p) =  rhod(p,l)*rhod_H(k,p)  ;
                B2d(lineIndex, p) = jj * rhod_H(p,l)*rhod(k,p) ;
            }
        }
    }

        // B = | B1u  iB1u  B2u   iB2u  |
    //     | B1d -iB1d  B2d  -iB2d  | :
    arma::cx_dmat B = arma::cx_dmat(2 * sizeMat, 4*nSites);
    B.submat(0,         0,          sizeMat-1,      nSites-1)   = B1u;
    B.submat(0,         nSites,     sizeMat-1,      2*nSites-1) = jj * B1u;
    B.submat(0,         2*nSites,   sizeMat-1,      3*nSites-1) = B2u;
    B.submat(0,         3*nSites,   sizeMat-1,      4*nSites-1) = jj*B2u;
    B.submat(sizeMat,   0,          2*sizeMat-1,    nSites-1)   = B1d;
    B.submat(sizeMat,   nSites,     2*sizeMat-1,    2*nSites-1) = - jj * B1d;
    B.submat(sizeMat,   2*nSites,   2*sizeMat-1,    3*nSites-1) = B2d;
    B.submat(sizeMat,   3*nSites,   2*sizeMat-1,    4*nSites-1) = - jj *B2d;

    B = cst * B;

    
    arma::cx_dvec UMat4 = arma::cx_dmat(4*nSites, 4*nSites, arma::fill::zeros);
    UMat4.submat(0       , 0        ,  nSites-1     , nSites-1)     = UMat;        // UMat4 =  | Uab   0    0    0   |
    UMat4.submat(nSites  , nSites   ,  2*nSites-1   , 2*nSites-1)   = UMat;        //          |  0   Uab   0    0   |  
    UMat4.submat(2*nSites, 2*nSites , 3*nSites-1    , 3*nSites-1)   = UMat;        //          |  0    0   Uab   0   |
    UMat4.submat(3*nSites, 3*nSites , 4*nSites-1    , 4*nSites-1)   = UMat;        //          |  0    0    0   Uab  |



    return B * UMat4 * B.st();

}

/// @brief Decomposition of a diffusion matrice into its diffusion rate
/// @param D the diffusion matrix to decompose 
/// @return B a square root of D, BB' = D
arma::cx_dmat FHGPSRSolver::sqrtDecomposition(arma::cx_dmat& D)
{

    unsigned int sizeD = D.n_rows;
    unsigned int sizeD_2 = sizeD/2;
    bool successSVD = false;

    if (D.has_nan() || D.has_inf())
    {
        return arma::datum::nan * arma::cx_dmat(sizeD, 1, arma::fill::zeros);
    }
    else
    {
        // D =  | 0     A |
        //      | A^T   0 |
        // D1 = -0.5i A
        arma::cx_dmat D1 = -0.5 * jj * D.submat(0,sizeD_2, sizeD_2-1, sizeD-1);

        arma::cx_dmat P;
        arma::dvec s;
        arma::cx_dmat Q;

        // D1 = U s V^T
        //bool SVDsuccess = arma::svd(U, s, V, D1);

        unsigned int nSites = static_cast<unsigned int>(std::sqrt(sizeD_2));
        if (diffusionIntegrator == GPSRDiffusionIntegrator::numerical_svd)
        {
            successSVD = svd(P,s,Q,D1);
        }
        else if (diffusionIntegrator == GPSRDiffusionIntegrator::numerical_rsvd)
        {
            //successSVD = rsvd(P,s,Q,D1, 2*nSites);
            successSVD = rsvd(P,s,Q,D1, nSites + 5);
        }

        arma::dmat sqrtS = arma::diagmat(arma::sqrt(s)); 

        // B =  | U.sqrtS           i.U.sqrtS  |
        //      | ii.c(V).sqrtS     c(V).sqrtS |
        arma::cx_dmat B = arma::cx_dmat(2*P.n_rows, 2*P.n_cols, arma::fill::zeros);
        B.submat(0,          0,           B.n_rows/2-1, B.n_cols/2-1)    = P * sqrtS;
        B.submat(0,          B.n_cols/2,  B.n_rows/2-1, B.n_cols-1)      = FHGPSRSolver::jj * P * sqrtS;
        B.submat(B.n_rows/2, 0,           B.n_rows-1,   B.n_cols/2-1)    = FHGPSRSolver::jj * arma::conj(Q) * sqrtS;  
        B.submat(B.n_rows/2, B.n_cols/2,  B.n_rows-1,   B.n_cols-1)      = arma::conj(Q) * sqrtS; 

        return B;
    }

}

/// @brief Decomposition of a symmetric and half-empty diffusion matrice into its diffusion rate
/// @param A the upper right quarter of the diffusion matrix to decompose 
/// @return B a square root of D, BB' = D
arma::cx_dmat FHGPSRSolver::sqrtDecomposition_fromA(arma::cx_dmat& A)
{
  unsigned int sizeD_2 = A.n_rows;
  uint sizeD = sizeD_2 * 2;
    bool successSVD = false;

    if (A.has_nan() || A.has_inf())
    {
        return arma::datum::nan * arma::cx_dmat(sizeD, 1, arma::fill::zeros);
    }
    else
    {
        // D =  | 0     A |
        //      | A^T   0 |
        // D1 = -0.5i A
        arma::cx_dmat D1 = -0.5 * jj * A;

        arma::cx_dmat P;
        arma::dvec s;
        arma::cx_dmat Q;

        // D1 = P s Q^T
        //bool SVDsuccess = arma::svd(P, s, Q, D1);

        unsigned int nSites = static_cast<unsigned int>(std::sqrt(sizeD_2));
        if (diffusionIntegrator == GPSRDiffusionIntegrator::numerical_svd)
        {
            successSVD = svd(P,s,Q,D1);
        }
        else if (diffusionIntegrator == GPSRDiffusionIntegrator::numerical_rsvd)
        {
            successSVD = rsvd(P,s,Q,D1, 2*nSites);
        }

        arma::dmat sqrtS = arma::diagmat(arma::sqrt(s)); 

        // B =  | U.sqrtS           i.U.sqrtS  |
        //      | ii.c(V).sqrtS     c(V).sqrtS |
        arma::cx_dmat B = arma::cx_dmat(2*P.n_rows, 2*P.n_cols, arma::fill::zeros);
        B.submat(0,          0,           B.n_rows/2-1, B.n_cols/2-1)    = P * sqrtS;
        B.submat(0,          B.n_cols/2,  B.n_rows/2-1, B.n_cols-1)      = FHGPSRSolver::jj * P * sqrtS;
        B.submat(B.n_rows/2, 0,           B.n_rows-1,   B.n_cols/2-1)    = FHGPSRSolver::jj * arma::conj(Q) * sqrtS;  
        B.submat(B.n_rows/2, B.n_cols/2,  B.n_rows-1,   B.n_cols-1)      = arma::conj(Q) * sqrtS; 

        return B;
    }
}


/// @brief Simulation time step
/// @param rhou the up-spin density matrix, to be updated
/// @param rhod the down-spin density matrix, to be updated 
/// @param d_t the time step
/// @return void
void FHGPSRSolver::step(double currentSimulationTime, arma::cx_dmat& rhou, arma::cx_dmat& rhod, double d_t)
{

    if (rhou.has_nan() || rhou.has_inf() || rhod.has_nan() || rhod.has_inf())
    {
        rhou = arma::datum::nan * rhou;
        rhod = arma::datum::nan * rhod;
    }
    else
    {
        arma::cx_dmat drhou_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
        arma::cx_dmat drhod_Drift = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

        arma::cx_dmat drhou_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
        arma::cx_dmat drhod_Diff = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

        // Step from the drift part : 
        if (driftIntegrator == GPSRDriftIntegrator::RungeKutta2)            step_Drift_RK2(currentSimulationTime, rhou, rhod, d_t, drhou_Drift, drhod_Drift);
        else if (driftIntegrator == GPSRDriftIntegrator::RungeKutta4_38)    step_Drift_RK4_38(currentSimulationTime, rhou, rhod, d_t, drhou_Drift, drhod_Drift);

        // Step from the diffusion part : 
        double U0 = std::real( fhSystem->getCoulombInteraction(0,0));
        if (diffusionIntegrator==GPSRDiffusionIntegrator::analytical)   step_Diffusion(rhou, rhod, U0, d_t, drhou_Diff, drhod_Diff);
        else step_Diffusion_svd(rhou, rhod, U0, d_t, drhou_Diff, drhod_Diff);

        rhou = rhou + drhou_Drift + drhou_Diff;
        rhod = rhod + drhod_Drift + drhod_Diff;
    }

}


/// @brief Launch simulation and write the results
void FHGPSRSolver::runSimulation( )
{

    //printf("GPSR dynamic, nTraj=%i \n", nTraj);
    std::cout << "   -- Start simulation." << std::endl;
    unsigned int nSites = fhSystem->getnSites();
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
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
        Energies    = arma::cx_dmat(4,nTimePoints, arma::fill::zeros);
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



// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel
{

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
        En_loc          = arma::cx_dmat(3,nTimePoints, arma::fill::zeros);
    }

    arma::cx_dmat Nu_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat Nd_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhou      = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhod      = arma::cx_dmat(nSites, nSites);

    double simulationTime = 0;

    #pragma omp for 
    for (unsigned int i_traj=0; i_traj <nTraj; i_traj ++)
    {
        // Noiseless initial rho matrices :
        rhou = rhou0;
        rhod = rhod0;

      
        // We compute the trajectory dynamic and add up the diagnoal to the result
        for (unsigned int t = 0;  t<nTimePoints; t++)
        {
            timeLine[t] = t*dt;
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
                En_loc(0,t)=  En_loc(0,t) + (arma::dot(rhou.diag(), *U * rhod.diag()) + arma::dot(rhod.diag(), *U * rhou.diag())) / arma::cx_double(2,0);
                En_loc(1,t) = En_loc(1,t) + arma::accu((rhou + rhod) % (*tM));
                if (!fhSystem->isExternalPotentialNull())
                {
                    En_loc(2,t) = En_loc(2,t) + arma::accu((rhou + rhod) % (*extPot));
                }
            }

            
            for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s++)
            {
                step(simulationTime, rhou, rhod, dt/nStepsPerTimePoint);
                simulationTime += dt/nStepsPerTimePoint;
            }
        }

        unsigned int tt = std::floor(nTraj/10);
        if ((i_traj%tt) == 0)
        {
            std::cout<< "     - Trajectory " << i_traj << " computed by thread " << omp_get_thread_num() << std::endl;    
        }
   
    }

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
            for (uint ie =0; ie < 3; ie ++)
            {
                Energies.row(ie) +=En_loc.row(ie);
            }
        }
    }
}
    

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
            Nu.slice(t) = rhou_t.slice(t) / nTraj_arma ;
            Nd.slice(t) = rhod_t.slice(t) / nTraj_arma;
        }
    }
    if (saveDiagDMatrix)
    {
        //std::cout << "Computation of Nii" << std::endl;
        Niiu = rhoiiu / nTraj_arma ;
        Niid = rhoiid / nTraj_arma;
    }
    if (saveDiag2DMatrix)
    {
        N2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        N2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            N2iiu.col(t) = rho2iiu.col(t) / nTraj_arma ;
            N2iid.col(t) = rho2iid.col(t) / nTraj_arma ;
        }
    }
    if (saveNiuNid)
    {
        NiuNid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        for (unsigned int t = 0; t < nTimePoints; t++)
        {
            NiuNid.col(t) = rhoiurhoid.col(t) / nTraj_arma;
        }
    }
    if (saveTBStates)
    {
        Mkku = mkku / nTraj_arma ;
        Mkkd = mkkd / nTraj_arma ;
    }
    if (saveMFStates)
    {
        Mukku = mukku / nTraj_arma;
        Mukkd = mukkd / nTraj_arma;
    }
    if (saveEnergies)
    {
        Energies = Energies/nTraj_arma;
        Energies.row(3) = arma::sum(Energies,0);
    }

    writeResults(Nu, Nd, Niiu, Niid, Mkku, Mkkd, Mukku, Mukkd, N2iiu, N2iid, NiuNid, Energies, timeLine);

}





/// @brief Copy of the MATLAB rsvd function.  random SVD
/// Extremely fast computation of the truncated Singular Value Decomposition, using
/// randomized algorithms as described in Halko et al. 'finding structure with randomness'.
/// Antoine Liutkus  (c) Inria 2014
/// @param U 
/// @param S matrix of singular values
/// @param V
/// @param D matrix whose svd we want
/// @param K number of compponents to keep
/// @return false if the decomposition failed
bool rsvd(arma::cx_dmat& U, arma::dvec& s, arma::cx_dmat& V, arma::cx_dmat& D, unsigned int K)
{
    bool successOrth = true;
    bool successSVD = true; 
    
                                                    // MATLAB version
    //unsigned int M = D.n_rows;
    unsigned int N = D.n_cols;                      //<- [M,N] = arma::size(D);

    //unsigned int P = std::min(2*K,N);               //<- P = min(2*K,N);
    unsigned int P = std::min(K+5,N);               //<- P = min(K+5,N);


    arma::dmat X = arma::dmat(N, P, arma::fill::randn); //<- X = randn(N,P);
    arma::cx_dmat Y = D*X;                          //<- Y = A*X;

    arma::cx_dmat W1;
    successOrth = arma::orth(W1, Y);               //<- W1 = orth(Y);

    arma::cx_dmat B = W1.t() * D;                   //<- B = W1'*A;
    
    arma::cx_dmat W2; 

    successSVD = arma::svd_econ(W2, s, V, B);       //<- [W2,S,V] = svd(B,'econ');

    U = W1 * W2;                                    //<- U = W1*W2;

    unsigned int uCols = U.n_cols;

    unsigned int K2 = std::min(K,uCols);            //<- K = min(K, size(U,2))

    U = U.head_cols(K2);                            //<- U = U(:,1:K);
    s =s.head(K2);                                  //<- S = S(1:K,1:K);
    V = V.head_cols(K2);                            //<- V=V(:,1:K);

    return successOrth && successSVD;
} 


FHGPSRSolver::~FHGPSRSolver()
{
}
