#include "../include/FHfTWAStateSearch.h"


FHfTWAStateSearch::FHfTWAStateSearch(FHSystem* fhS) 
{

    fhSystem = fhS;

    nPoints = 0;

    nMaxPoints = 1000;
    nMaxTries = 100000;


    unsigned int nSites = fhSystem->getnSites();
    arma::dmat zeroMat = arma::dmat(nSites, nSites, arma::fill::zeros);

    halfDiag = 0.5 * arma::cx_dmat(nSites, nSites, arma::fill::eye);
    halfVec = 0.5 * arma::cx_dvec(nSites, arma::fill::ones);

    std::cout <<  "---- Init Search 1 "<< std::endl;

    U = fhSystem->getUMatrixPointer();
    tM = fhSystem->gettMatrixPointer();
    //J = fhSystem->getJMatrixPointer();
    V = fhSystem->getVMatrixPointer();
    extPot = fhSystem->getExtPotPointer();

    std::cout <<  "---- Init Search 2"<< std::endl;



}



void                        FHfTWAStateSearch::findDensity()
{

    arma::cx_double EtH, EU, Epe;
    while( (nPoints < nMaxPoints) && (nTries < nMaxTries) )
    {


        makeTry(rhou0, rhod0, rhou, rhod);
        nTries++;


        energy(0, rhou, rhod, EtH, EU, Epe);

        double eTry = (EtH + EU + Epe).real();

        if (eTry < meanE)
        {
            rhou_mean = ((rhou_mean * nPoints) + rhou) / (double)(nPoints+1);
            rhod_mean = ((rhod_mean * nPoints) + rhod) / (double)(nPoints+1);

            meanE = ((meanE * nPoints) + eTry) / (double)(nPoints+1);

            if (eTry < minE)
            {
                minE = eTry;
            }

            std::cout << "   --- try " << nTries << ", new Point, n " << nPoints << ", new E=" << meanE << std::endl ;

            nPoints ++;
        }


    }
}


void                        FHfTWAStateSearch::initiateStateSearch()
{

    fhSystem->computeTBGroundState();

    uint nSites = fhSystem->getnSites();

    rhou0 = arma::cx_dmat( fhSystem->get_tbgs_Nu_Pointer()->memptr(), nSites, nSites);
    rhod0 = arma::cx_dmat( fhSystem->get_tbgs_Nd_Pointer()->memptr(), nSites, nSites);

    rhou_mean = rhou0;
    rhod_mean = rhod0;


    arma::cx_double EtH, EU, Epe;
    energy(0, rhou0, rhod0, EtH, EU, Epe);
    double e = (EtH + EU + Epe).real();

    maxE = e;
    minE = e;
    meanE = e;

    nPoints = 1;
    nTries = 1;



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


/// @brief Time derivative of the trajectories
/// @param rhou the up-spin density matrix
/// @param rhod the down-spin density matrix 
/// @param drhou_dt the time derivative of the up-spin density matrix, to be updated
/// @param drhod_dt the time derivative of the down-spin density matrix, to be updated 
/// @return void
void FHfTWAStateSearch::derivative(double currentSimulationTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_dmat& drhou_dt, arma::cx_dmat& drhod_dt)
{

    const arma::cx_double       jj = arma::cx_double(0,1);

    arma::cx_dmat DeltaUp = *tM + arma::diagmat(*U * (rhod.diag() + halfVec)) ;
    arma::cx_dmat DeltaDo = *tM + arma::diagmat(*U * (rhou.diag() + halfVec)) ;  
    
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


void FHfTWAStateSearch::energy(double currentTime, const arma::cx_dmat& rhou, const arma::cx_dmat& rhod, arma::cx_double& EtH, arma::cx_double& EU, arma::cx_double& Epe)
{
    arma::cx_dmat Nu_loc =  rhou + halfDiag;
    arma::cx_dmat Nd_loc =  rhod + halfDiag;
    EU =  (arma::dot(Nu_loc.diag(), *U * Nd_loc.diag()) + arma::dot(Nd_loc.diag(), *U * Nu_loc.diag())) / arma::cx_double(2,0);
    EtH = arma::accu((Nu_loc + Nd_loc) % (*tM));
    Epe = 0;
    if (!fhSystem->isExternalPotentialNull())
    {
        Epe =  arma::accu((Nu_loc + Nd_loc) % (*extPot));
    }
}




void FHfTWAStateSearch::makeTry( const arma::cx_dmat& rhou0, const arma::cx_dmat& rhod0,  arma::cx_dmat& rhou, arma::cx_dmat& rhod)
{
    rhou = rhou0;
    rhod = rhod0;

    uint nSites = fhSystem->getnSites();

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


    rhou = rhou + (*fhSystem->get_tbgs_V_Pointer()) * wnu * ((*fhSystem->get_tbgs_V_Pointer()).t()) ;
    rhod = rhod + (*fhSystem->get_tbgs_V_Pointer()) * wnd * ((*fhSystem->get_tbgs_V_Pointer()).t()) ;



}


FHfTWAStateSearch::~FHfTWAStateSearch() 
{
}