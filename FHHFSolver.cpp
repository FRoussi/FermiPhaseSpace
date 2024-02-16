#include "FHHFSolver.h"


FHHFSolver::FHHFSolver(FHSystem* fhS)
{
    fhSystem = fhS;
    U = fhSystem->getUMatrixPointer();
    J = fhSystem->getJMatrixPointer();
}

void FHHFSolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
    writeResults = true;
}


void FHHFSolver::setTimeParameters(double t_max, double d_t, unsigned int nSPTP)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nSPTP;
    this->nTimePoints = floor(tMax/dt +1);
}

void FHHFSolver::initiateSolver()
{

    Nu0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    Nd0 = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    for (unsigned int it= 0; it < fhSystem->getnSites(); it ++)
    {
        Nu0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_upSpin()[it],   0) ;
        Nd0(it, it) = arma::cx_double(1 * fhSystem->getInitialCondition_downSpin()[it], 0) ;
    }

    printf("Test 4 - Initial Density matrix : \n");
    Nu0.print();


}

void FHHFSolver::step_RK1(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(Nut, Ndt, dNu_dt, dNd_dt);

    Nut = Nut + dNu_dt * dt ;
    Ndt = Ndt + dNd_dt * dt ;
}

void FHHFSolver::step_RK2(double dt)
{
    arma::cx_dmat dNu_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);
    arma::cx_dmat dNd_dt = arma::cx_dmat(fhSystem->getnSites(), fhSystem->getnSites(), arma::fill::zeros);

    derivative(Nut, Ndt, dNu_dt, dNd_dt);

    arma::cx_dmat Nu_mid = Nut + dNu_dt * dt/2 ;
    arma::cx_dmat Nd_mid = Ndt + dNd_dt * dt/2 ;

    derivative(Nu_mid, Nd_mid, dNu_dt, dNd_dt);

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
    derivative(Nut, Ndt, dNu_dt_k1, dNd_dt_k1);
// k2 : 0.5 | 0.5
    arma::cx_dmat Nu_k1 = Nut + dNu_dt_k1 * dt/2 ;
    arma::cx_dmat Nd_k1 = Ndt + dNd_dt_k1 * dt/2 ;
    derivative(Nu_k1 , Nd_k1 , dNu_dt_k2, dNd_dt_k2);
// k3 : 0.5 |  0   0.5
    arma::cx_dmat Nu_k2 = Nut + dNu_dt_k2 * dt/2 ;
    arma::cx_dmat Nd_k2 = Ndt + dNd_dt_k2 * dt/2 ;
    derivative(Nu_k2 , Nd_k2, dNu_dt_k3, dNd_dt_k3);
// k3 : 1   |  0    0    1
    arma::cx_dmat Nu_k3 = Nut + dNu_dt_k3 * dt ;
    arma::cx_dmat Nd_k3 = Ndt + dNd_dt_k3 * dt ;
    derivative(Nu_k3 , Nd_k3, dNu_dt_k4, dNd_dt_k4);
// k4 :     |  1/6  2/6  2/6  1/6  
    Nut = Nut + (dNu_dt_k1 + 2* dNu_dt_k2 + 2*dNu_dt_k3 + dNu_dt_k4) * dt/6;
    Ndt = Ndt + (dNd_dt_k1 + 2* dNd_dt_k2 + 2*dNd_dt_k3 + dNd_dt_k4) * dt/6;
}

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
    derivative(Nut, Ndt, dNu_dt_k1, dNd_dt_k1);
// k2 : 1/3  |  1/3
    arma::cx_dmat Nu_k = Nut + dNu_dt_k1 * dt/3 ;
    arma::cx_dmat Nd_k = Ndt + dNd_dt_k1 * dt/3 ;
    derivative(Nu_k , Nd_k , dNu_dt_k2, dNd_dt_k2);
// k3 : 2/3  |  -1/3  1
    Nu_k = Nut + (- dNu_dt_k1/3 + dNu_dt_k2) * dt ;
    Nd_k = Ndt + (- dNd_dt_k1/3 + dNd_dt_k2) * dt ;
    derivative(Nu_k , Nd_k, dNu_dt_k3, dNd_dt_k3);
// k3 :  1   |  1     -1    1
    Nu_k = Nut + (dNu_dt_k1 - dNu_dt_k2 + dNu_dt_k3) * dt;
    Nd_k = Ndt + (dNd_dt_k1 - dNd_dt_k2 + dNd_dt_k3) * dt;
    derivative(Nu_k , Nd_k, dNu_dt_k4, dNd_dt_k4);
//    :      | 1/8    3/8   3/8   1
    Nut = Nut + (dNu_dt_k1 + dNu_dt_k2*3 + dNu_dt_k3*3 + dNu_dt_k4) * dt/8;
    Ndt = Ndt + (dNd_dt_k1 + dNd_dt_k2*3 + dNd_dt_k3*3 + dNd_dt_k4) * dt/8;
}

void FHHFSolver::step_RK45(double dt)
{
}

void FHHFSolver::derivative(const arma::cx_dmat& nut, const arma::cx_dmat& ndt, arma::cx_dmat& dnut_dt, arma::cx_dmat& dndt_dt)
{

    arma::cx_dmat DeltaUp = *J + arma::diagmat(*U * ndt.diag());
    arma::cx_dmat DeltaDo = *J + arma::diagmat(*U * nut.diag());

    dnut_dt = - jj * ( nut * DeltaUp - DeltaUp * nut) ;
    dndt_dt = - jj * ( ndt * DeltaDo - DeltaDo * ndt) ;
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
            step_RK4_38(dt/nStepsPerTimePoint);
        }
 
    }


// Computing other values :
    arma::cx_dmat Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_dmat Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EU =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_drowvec EtH =  arma::cx_drowvec(nTimePoints, arma::fill::zeros);
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

        EtH(t) =  arma::accu(*J % Nu.slice(t)) + arma::accu(*J % Nd.slice(t));
        EU(t) = (arma::cdot(diagNu, *U * diagNd) + arma::cdot(diagNd,  *U * diagNu)) / arma::cx_double(2,0);

    }

    arma::cx_drowvec E = EU + EtH;

    printf("--- Mean field dynamic, n1u : \n");
    (arma::real(Niiu.row(0))).print();
    printf("--- Mean field dynamic, E : \n");
    (arma::real(E)).print();


       // Writing the results :
    if (writeResults)
    {
        ResultWriter results = ResultWriter(workName, folderName);
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
        results.writeLineInReal("EU", EU);
        results.writeLineInReal("EtH", EtH);
        results.writeLineInReal("E", E);
    }


    

}

FHHFSolver::~FHHFSolver()
{

}
 