#include "FHHFSolver.h"


FHHFSolver::FHHFSolver(FHSystem* fhS)
{
    fhSystem = fhS;
    U = fhSystem->getUMatrixPointer();
    J = fhSystem->getJMatrixPointer();
}



void FHHFSolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
}

void FHHFSolver::initiateSolver(std::vector<unsigned int> pPartU, std::vector<unsigned int> pPartD )
{

    Nu0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    Nd0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    for (std::vector<unsigned int>::iterator it_u = pPartU.begin(); it_u < pPartU.end(); it_u++)
    {
        Nu0(*it_u, *it_u) = arma::cx_double(1,0);
    }

    for (std::vector<unsigned int>::iterator it_d = pPartD.begin(); it_d < pPartD.end(); it_d++)
    {
        Nd0(*it_d, *it_d) = arma::cx_double(1,0);
    }




    printf("Test 4 - Initial Density matrix : \n");
    Nu0.print();


}

void FHHFSolver::step_RK1(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(&Nut, &Ndt, &dNu_dt, &dNd_dt);

    Nut = Nut + dNu_dt * dt ;
    Ndt = Ndt + dNd_dt * dt ;
}

void FHHFSolver::step_RK2(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(&Nut, &Ndt, &dNu_dt, &dNd_dt);

    arma::cx_dmat Nu_mid = Nut + dNu_dt * dt/2 ;
    arma::cx_dmat Nd_mid = Ndt + dNd_dt * dt/2 ;

    derivative(&Nu_mid, &Nd_mid, &dNu_dt, &dNd_dt);

    Nut = Nut + dNu_dt * dt ;
    Ndt = Ndt + dNd_dt * dt ;
}

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
    derivative(&Nut, &Ndt, &dNu_dt_k1, &dNd_dt_k1);
// k2 :
    arma::cx_dmat Nu_k1 = Nut + dNu_dt_k1 * dt/2 ;
    arma::cx_dmat Nd_k1 = Ndt + dNd_dt_k1 * dt/2 ;
    derivative(&Nu_k1 , &Nd_k1 , &dNu_dt_k2, &dNd_dt_k2);
// k3 :
    arma::cx_dmat Nu_k2 = Nut + dNu_dt_k2 * dt/2 ;
    arma::cx_dmat Nd_k2 = Ndt + dNd_dt_k2 * dt/2 ;
    derivative(&Nu_k2 , &Nd_k2, &dNu_dt_k3, &dNd_dt_k3);
// k3 :
    arma::cx_dmat Nu_k3 = Nut + dNu_dt_k2 * dt ;
    arma::cx_dmat Nd_k3 = Ndt + dNd_dt_k2 * dt ;
    derivative(&Nu_k3 , &Nd_k3, &dNu_dt_k4, &dNd_dt_k4);

    Nut = Nut + (dNu_dt_k1 + 2* dNu_dt_k2 + 2*dNu_dt_k3 + dNu_dt_k4) * dt/6;
    Ndt = Ndt + (dNd_dt_k1 + 2* dNd_dt_k2 + 2*dNd_dt_k3 + dNd_dt_k4) * dt/6;
}

void FHHFSolver::step_RK45(double dt)
{
}

void FHHFSolver::derivative(arma::cx_dmat* nut, arma::cx_dmat* ndt, arma::cx_dmat* dnut_dt, arma::cx_dmat* dndt_dt)
{

    arma::cx_dmat DeltaUp = *J + arma::diagmat(*U * ndt->diag());
    arma::cx_dmat DeltaDo = *J + arma::diagmat(*U * nut->diag());

    *dnut_dt = - jj * ( *nut * DeltaUp - DeltaUp * *nut) ;
    *dndt_dt = - jj * ( *ndt * DeltaDo - DeltaDo * *ndt) ;
}


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
            step_RK2(dt/nStepsPerTimePoint);
        }
 
    }


// Computing other values :
    arma::cx_dmat Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EU =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EJ =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec nPartU =  arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec nPartD =  arma::drowvec(nTimePoints, arma::fill::zeros);
    for (unsigned int t = 0;  t<nTimePoints; t++)
    {
        arma::cx_dvec diagNu = Nu.slice(t).diag();
        arma::cx_dvec diagNd = Nd.slice(t).diag();
        Niiu.col(t) = diagNu;
        Niid.col(t) = diagNd;
        nPartU(t) = std::real(arma::accu(diagNu));
        nPartD(t) = std::real(arma::accu(diagNd));

        EJ(t) =  arma::accu(*J % Nu.slice(t)) + arma::accu(*J % Nd.slice(t));
        EU(t) = (arma::cdot(diagNu, *U * diagNd) + arma::cdot(diagNd,  *U * diagNu)) / arma::cx_double(2,0);

    }

    arma::cx_drowvec E = EU + EJ;

    printf("--- Mean field dynamic, n1u : \n");
    (arma::real(Niiu.row(0))).print();
    printf("--- Mean field dynamic, E : \n");
    (arma::real(E)).print();

}

FHHFSolver::~FHHFSolver()
{
}
 