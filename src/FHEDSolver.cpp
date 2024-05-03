#include "../include/FHEDSolver.h"

FHEDSolver::FHEDSolver(FHSystem *fhS)
{
    fhSystem = fhS;
    writeResults = false;
    currentSimulationTime = 0;
}

void FHEDSolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
    writeResults = true;
}


void FHEDSolver::setTimeParameters(double t_max, double d_t, uint nPerTS)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nPerTS;
    this->nTimePoints = floor(tMax / dt + 1);
}

void FHEDSolver::setComputeMomenta(bool yn)
{
    computeMomenta = yn;
    if (yn)
    {
        if (basisSize > 1000)
        {
            std::cout << "! Caution, The ED basis size is " << basisSize << ", the computation of the momentum space could be long." << std::endl << std::flush;
        }

        nOfMomentumStates = std::min(basisSize, (uint)10);

        arma::sp_dmat H_t_real = arma::real(H_t);
        arma::eigs_sym(momentaNRJ, momentaVec, H_t_real, nOfMomentumStates);  
    }
}


/// @brief Initiate the ED solver, compute the initial wave-function, then the usefull operators.
///@return void
void FHEDSolver::initiateSolver()
{
    this->posPartU = bitsetFromVector(fhSystem->getInitialCondition_upSpin());
    this->posPartD = bitsetFromVector(fhSystem->getInitialCondition_downSpin());
    this->nPartU = posPartU.count();
    this->nPartD = posPartD.count();


    this->basisSizeSpinUp = nChoosep(fhSystem->getnSites(), nPartU);
    this->basisSizeSpinDown = nChoosep(fhSystem->getnSites(), nPartD);
    this->basisSize = basisSizeSpinUp * basisSizeSpinDown;
    this->basisArray = std::vector<std::array<std::bitset<MAX_BITS>, 2>>(basisSize); // 1 basis vector = 1 combination of up-spin position and 1 combination of down-spin position
    
    std::vector<std::vector<bool>> c_puVec = getCombinations(this->fhSystem->getnSites(), this->nPartU, this->fhSystem->getnSites());
    c_puVec.shrink_to_fit();
    this->c_pu = std::vector<std::bitset<MAX_BITS>>(basisSizeSpinUp);
    for (uint it = 0 ; it<basisSizeSpinUp; it++)
    {
        c_pu[it] = bitsetFromVector(c_puVec[it]);
    }
    fusionSortVector(c_pu);
    

    std::vector<std::vector<bool>> c_pdVec = getCombinations(this->fhSystem->getnSites(), this->nPartD, this->fhSystem->getnSites());
    c_pdVec.shrink_to_fit();
    this->c_pd = std::vector<std::bitset<MAX_BITS>>(basisSizeSpinDown);


    for (uint it = 0 ; it<basisSizeSpinDown; it++)
    {
        c_pd[it] = bitsetFromVector(c_pdVec[it]);
    }

    fusionSortVector(c_pd);


    uint countor = 0;
    for (uint l = 0; l < basisSizeSpinUp; l++)
    {
        for (uint k = 0; k < basisSizeSpinDown; k++)
        {
            basisArray[countor][0] = c_pu[l];         //<- up-spin particle combinations of position
            basisArray[countor][1] = c_pd[k];         //<- down-spin particle combinations of position
            countor += 1;
        }
    }

    uint l = getBasisIndex(this->posPartU, this->posPartD);
    this->psi_0 = arma::cx_dvec(this->basisSize, arma::fill::zeros);
    this->psi_0[l] = 1;                         //<- Compute the initial wave-function

    this->psi = this->psi_0;                    //< The current wave-function becomes the initial one.


    buildPairOperators();                       //<- Usefull pair operators
    buildHamiltonian();                         //<- Hamiltonian

}


///@brief Compute the system's pair operators. First the number operators, then the important tunneling operators.
///@return void
void FHEDSolver::buildPairOperators()
{
    uint nSites = fhSystem->getnSites();
    NijuSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);
    NijdSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);

// Filling the 2nd dimention of the vector-vector and computing the number operators
    for (uint i = 0; i < nSites; i++)
    {
        NijuSet[i] = std::vector<arma::sp_cx_dmat>(nSites);
        NijuSet[i][i] = numberOperatorSpinUp(i);

        NijdSet[i] = std::vector<arma::sp_cx_dmat>(nSites);
        NijdSet[i][i] = numberOperatorSpinDown(i);

    }


    for (uint i = 0; i < nSites; i++)
    {
        for (uint j = i + 1; j < nSites; j++)
        {

            if ((fhSystem->getTunnelingEnergy(i, j) != arma::cx_double(0,0)) || (fhSystem->getExchangeCouplingEnergy(i, j) != arma::cx_double(0,0)))
            {

                NijuSet[i][j] = tunnelingOperatorSpinUp(i, j);
                NijdSet[i][j] = tunnelingOperatorSpinDown(i, j);
                NijuSet[j][i] = NijuSet[i][j].t();
                NijdSet[j][i] = NijdSet[i][j].t();

            }
        }
    }

    /*printf("Some pair operators :");
    NijuSet[0][1].print_dense();
    NijuSet[1][0].print_dense();*/
}


///@brief Compute the system's Hamiltonian interaction part
///@return void
void FHEDSolver::buildHamiltonian_U()
{

    H_U = arma::sp_cx_dmat(basisSize, basisSize);
    for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite ++)
    {
        for (uint jSite = 0 ; jSite < fhSystem->getnSites(); jSite ++)
        {
            if(fhSystem->getInteractionEnergy(iSite,jSite) != arma::cx_double(0,0))
            {
                H_U = H_U + fhSystem->getInteractionEnergy(iSite,jSite) * NijuSet[iSite][iSite] * NijdSet[jSite][jSite];
            }
        }
    }

}


///@brief Compute the system's Hamiltonian tunneling part
///@return void
void FHEDSolver::buildHamiltonian_tM()
{
    H_t = arma::sp_cx_dmat(basisSize, basisSize);
    for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite++)
    {
        for (uint jSite = 0; jSite < fhSystem->getnSites(); jSite++)
        {
            if (fhSystem->getTunnelingEnergy(iSite, jSite) != arma::cx_double(0,0))
            {
                H_t = H_t + fhSystem->getTunnelingEnergy(iSite, jSite) * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]) ; //<- Up-spin AND down-spin !!!
            }
        }
    }

}


///@brief Compute the system's Hamiltonian tunneling part
///@return void
void FHEDSolver::buildHamiltonian_J()
{
    H_J = arma::sp_cx_dmat(basisSize, basisSize);
    for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite++)
    {
        for (uint jSite = 0; jSite < fhSystem->getnSites(); jSite++)
        {
            if (fhSystem->getExchangeCouplingEnergy(iSite, jSite) != arma::cx_double(0,0))
            {
                H_J = H_J + fhSystem->getExchangeCouplingEnergy(iSite, jSite) * arma::cx_double(0.5,0) * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]).t() * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]) ;
            }
        }
    }
}


void FHEDSolver::buildHamiltonian_ExternalPotential(double time)
{
   H_externalPotential = arma::sp_cx_dmat(basisSize, basisSize);
   if (!fhSystem->isExternalPotentialNull())
   {
       for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite++)
       {
           for (uint jSite = 0; jSite < fhSystem->getnSites(); jSite++)
           {
               if (fhSystem->getExternalPotential(iSite, jSite, time) != arma::cx_double(0,0))
               {
                   H_externalPotential = H_externalPotential + fhSystem->getExternalPotential(iSite, jSite, time) * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]);
               }
           }
       }
   }
}


///@brief Compute the system's Hamiltonian
///@return void
void FHEDSolver::buildHamiltonian()
{
    buildHamiltonian_U();
    buildHamiltonian_tM();
    buildHamiltonian_J();

    buildHamiltonian_ExternalPotential(0);
    
    H = H_U + H_t + H_J + H_externalPotential;

    //printf("Test 2 - small system Hamiltonian : \n");
    //H.print_dense();
    //Hdense = arma::cx_dmat(H);
}


void FHEDSolver::updateHamiltonian(double time)
{
    buildHamiltonian_ExternalPotential(time);

    H = H_U + H_t + H_J + H_externalPotential;
}


///@brief Compute the operator "$\\hat{n}_{ii\\uparrow}$"
///@param iSite Site number
///@return The operator "$\\hat{n}_{ii\\uparrow}$"
arma::sp_cx_dmat FHEDSolver::numberOperatorSpinUp(uint iSite)
{
    arma::umat locations = arma::umat(2, basisSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(basisSize, arma::fill::zeros);
    uint index = 0;

    if (nPartU != 0)
    {
        for (uint lu = 0; lu < basisSizeSpinUp; lu++)
        {
            if (c_pu[lu][iSite])
            {
                for (uint ld = 0; ld < basisSizeSpinDown; ld++)
                {
                    locations(0, index) = ld + basisSizeSpinDown * lu;
                    locations(1, index) = ld + basisSizeSpinDown * lu;
                    values(index) = 1;
                    index++;
                }
            }
        }
    }
    return arma::sp_cx_dmat(locations, values, basisSize, basisSize);
}


///@brief Compute the operator "$\\hat{n}_{ii\\downarrow}$"
///@param iSite Site number
///@return The operator "$\\hat{n}_{ii\\downarrow}$"
arma::sp_cx_dmat FHEDSolver::numberOperatorSpinDown(uint iSite)
{
    arma::umat locations = arma::umat(2, basisSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(basisSize, arma::fill::zeros);
    uint index = 0;

    if (nPartD != 0)
    {
        for (uint ld = 0; ld < basisSizeSpinDown; ld++)
        {
            if (c_pd[ld][iSite])
            {
                for (uint lu = 0; lu < basisSizeSpinUp; lu++)
                {
                    locations(0, index) = ld + basisSizeSpinDown * lu;
                    locations(1, index) = ld + basisSizeSpinDown * lu;
                    values(index) = 1;
                    index++;
                }
            }
        }
    }
    return arma::sp_cx_dmat(locations, values, basisSize, basisSize);
}


///@brief Compute the operator "$\\hat{n}_{ij\\uparrow}$"
///@param iSite Site number
///@param jSite Site number
///@return The operator "$\\hat{n}_{ij\\uparrow}$"
arma::sp_cx_dmat FHEDSolver::tunnelingOperatorSpinUp(uint iSite, uint jSite)
{
    uint halfMatSizeUp = floor(basisSize * basisSizeSpinUp / 2);
    arma::umat locations = arma::umat(2, halfMatSizeUp, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(halfMatSizeUp, arma::fill::zeros);

    uint index = 0;

    for (uint lu = 0; lu < basisSizeSpinUp; lu++)
    {
        if (c_pu[lu][iSite] && (!c_pu[lu][jSite]))          // <- has the particle the possiblity to tunnel from iSite to jSite ?
        {
            std::bitset<MAX_BITS> c_pu_hoppdone = c_pu[lu];     // <- if yes, we write the new configuration in 'c_pu_hoppdone'
            c_pu_hoppdone[iSite] = false;
            c_pu_hoppdone[jSite] = true;
            for (uint mu = 0; mu < basisSizeSpinUp; mu++)
            {
                if (c_pu_hoppdone == c_pu[mu])              // <-we look for the basis function index of the new configuration, mu
                {
                    for (uint ld = 0; ld < basisSizeSpinDown; ld++) // <- we note it in the whole base, with the down spin particle part equals
                    {
                        locations(0,index) = basisSizeSpinDown * lu  + ld;      
                        locations(1,index) = basisSizeSpinDown * mu  + ld;
                        values(index) = pow(-1, getParticleNumberBetween(c_pu[lu], iSite, jSite));
                        index ++;
                    }

                }
            }
        }
    }

    return arma::sp_cx_dmat(locations, values, basisSize, basisSize);
}


///@brief Compute the operator "$\\hat{n}_{ij\\downarrow}$"
///@param iSite Site number
///@param jSite Site number
///@return The operator "$\\hat{n}_{ij\\downarrow}$"
arma::sp_cx_dmat FHEDSolver::tunnelingOperatorSpinDown(uint iSite, uint jSite)
{
     uint halfMatSizeDown = floor(basisSize * basisSizeSpinDown / 2);
    arma::umat locations = arma::umat(2, halfMatSizeDown, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(halfMatSizeDown, arma::fill::zeros);

    uint index = 0;

    for (uint ld = 0; ld < basisSizeSpinDown; ld++)
    {
        if (c_pd[ld][iSite] && (!c_pd[ld][jSite]))          // <- has the particle the possiblity to tunnel from iSite to jSite ?
        {
            std::bitset<MAX_BITS> c_pd_hoppdone = c_pd[ld];     // <- if yes, we write the new configuration in 'c_pd_hoppdone'
            c_pd_hoppdone.reset(iSite);                         //<- jSite to false / 0
            c_pd_hoppdone.set(jSite);                           //<- iSite to true / 1
            for (uint md = 0; md < basisSizeSpinDown; md++)
            {
                if (c_pd_hoppdone == c_pd[md])               // <-we look for the basis function index of the new configuration, md
                {
                    for (uint lu = 0; lu < basisSizeSpinUp; lu++) // <- we note it in the whole base, with the up spin particle part equals
                    {
                        locations(0,index) = basisSizeSpinDown * lu  + ld;
                        locations(1,index) = basisSizeSpinDown * lu  + md;
                        values(index) = pow(-1, getParticleNumberBetween(c_pd[ld], iSite, jSite));
                        index ++;
                    }

                }
            }
        }
    }

    return arma::sp_cx_dmat(locations, values, basisSize, basisSize);
}


///@brief Find the index of a wave-function basis vector
///@param posPU Positions of up-spin particles in the wave-function basis vector
///@param posPD Positions of down-spin particles in the wave-function basis vector
///@return The index of the corresponding wave-function basis vector, 0 if not found. 
uint FHEDSolver::getBasisIndex(std::bitset<MAX_BITS>& posPU, std::bitset<MAX_BITS>& posPD)
{
    uint lu = 0;
    bool index1Found = false;
    while (!index1Found && lu < this->basisSizeSpinUp)
    {
        if (posPU == c_pu[lu])
            index1Found = true;
        else
            lu++;
    }
    uint ld = 0;
    bool index2Found = false;
    while (!index2Found && ld < this->basisSizeSpinDown)
    {
        if (posPD == c_pd[ld])
            index2Found = true;
        else
            ld++;
    }
    if (index1Found && index2Found)
    {
        //std::cout << "Looking for the basis index, lu = " << lu << ", ld = " << ld << std::endl;
        return lu * this->basisSizeSpinDown + ld;
    }
    else
    {
        printf("The basis vector of a particle configuration cannot be found. \n");
        return 0;
    }
}

///@brief Contribution of the interaction of opposite-spin electrons of basis wave-function
///@param phi_l the index of the wave-function basis
///@return Opposite-spin interaction energy of basis wave-function phi_l 
std::complex<double> FHEDSolver::U_os_psil(uint phi_l)
{
    std::bitset<MAX_BITS> *posPU = &basisArray[phi_l][0];
    std::bitset<MAX_BITS> *posPD = &basisArray[phi_l][1];

    std::complex<double> u = 0;
    for (uint i_u = 0; i_u < fhSystem->getnSites(); i_u++)
    {
        if ((*posPU)[i_u])
        {
            for (uint i_d = 0; i_d < fhSystem->getnSites(); i_d++)
            {
                if ((*posPD)[i_d])
                    u = u + fhSystem->getInteractionEnergy(i_u, i_d);
            }
        }
    }
    return u;
}


/// @brief Simulation step, update the wave-function psi.
/// @param d_t Time step
/// @return void
void FHEDSolver::step(double d_t)
{
    
    if (fhSystem->isSystemTimeDependant())
    {
        updateHamiltonian(currentSimulationTime);
    }

    std::complex<double> alpha = std::complex<double>(0, -d_t); // <- alpha = -i dt / h

    this->psi = expmv(alpha, H, this->psi);
    currentSimulationTime +=dt;
}


/// @brief Launch simulation and write the results
void FHEDSolver::runSimulation()
{

    std::cout << "Start simulation" << std::endl << std::flush;


    std::vector<double> timeLine = std::vector<double>(nTimePoints, 0);
    arma::dmat niiu = arma::dmat(fhSystem->getnSites(), nTimePoints, arma::fill::zeros);
    arma::dmat niid = arma::dmat(fhSystem->getnSites(), nTimePoints, arma::fill::zeros);
    arma::drowvec EU = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec EtH = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec E = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::cx_dmat mk; 

    if (computeMomenta)
    {
        mk = arma::cx_dmat(nOfMomentumStates, nTimePoints, arma::fill::zeros);
    }

    this->psi = this->psi_0;

    // Bad filling thing :
    uint tt = nTimePoints/20;
    std::cout << "|                     |" << std::endl;
    std::cout << "|" << std::flush;

    for (uint t = 0; t < nTimePoints; t++)
    {
        timeLine[t] = t * dt;
        for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite++)
        {
            niiu(iSite,t) = static_cast<double>(std::real(arma::cdot(psi, NijuSet[iSite][iSite] * psi)));
            niid(iSite,t) = static_cast<double>(std::real(arma::cdot(psi, NijdSet[iSite][iSite] * psi)));
        }     
        EU[t] = static_cast<double>(std::real(arma::cdot(psi, H_U * psi)));
        EtH[t] = static_cast<double>(std::real(arma::cdot(psi, H_t * psi)));
        E[t] = EU[t] + EtH[t];

        if (computeMomenta)
        {
           mk.col(t) = arma::dot( momentaVec,  psi);
        }     

        step(this->dt);

        if ((t%tt) == 0) {  std::cout << "=" << std::flush ;   }

    }

    std::cout << "|" << std::endl << std::flush; 

    printf("Test 3 - small system dynamic, n_1u : \n");
    niiu.row(0).print();
    

    // Writing the results :
    if (writeResults)
    {
        ResultWriter results = ResultWriter(workName, folderName);
        results.newResults();
        results.writeLine("t", timeLine);
        for (uint iSite=0; iSite<fhSystem->getnSites(); iSite++)
        {
            results.writeLine("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", niiu.row(iSite));
        }
        for (uint iSite=0; iSite<fhSystem->getnSites(); iSite++)
        {
            results.writeLine("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_d", niid.row(iSite));
        }
        if (computeMomenta)
        {
            for (uint kMom = 0; kMom < nOfMomentumStates ; kMom++)
            {
                results.writeLineInReal("M_" + std::to_string(kMom), mk.row(kMom));
            }
        }
        results.writeLine("EU", EU);
        results.writeLine("EtH", EtH);
        results.writeLine("E", E);
    }


    

}


/// @brief Destructor
FHEDSolver::~FHEDSolver()
{
}

// ---------------------------
// Other functions


/// @brief Binomial coefficient 
/// @param n Number of possible choices - here sites
/// @param p Number of selected choices - here particles
/// @return Number of non-ordered combinations 
uint nChoosep(uint n, uint p)
{
    uint num = 1;
    uint den = 1;
    if (p < (n - p))
    {
        for (uint k = 1; k <= p; k++)
            den *= k;

        for (uint k = n - p + 1; k <= n; k++)
            num *= k;
    }
    else
    {
        for (uint k = 1; k <= n - p; k++)
            den *= k;

        for (uint k = p + 1; k <= n; k++)
            num *= k;
    }

    return (uint)num / den;
}





/// @brief Builds a vector of boolean vectors of all combination of p particles in n sites. 
/// There are no check on the entrance variables. This programm will malfunction if there are not adapted, like if n<p.
/// @param n number of sites
/// @param p number of particles
/// @param max_n is the final size of the vectors, the total number of sites. It does not crash if it's wrong, it will just be less efficient.
/// @return All combination of p particles in n sites.
std::vector<std::vector<bool>> getCombinations(uint n, uint p, uint max_n)
{
    // getCombinations is a recursive function:
    if (p == 0) // if there is 0 particles, we return the vector of the empty vector.
    {
        std::vector<bool> emptyVect = std::vector<bool>(n, false);
        emptyVect.reserve(max_n);
        return std::vector<std::vector<bool>>(1, emptyVect);
    }
    else if (p == n) // if there are as many particles as sites, we return the vector of the filled vector.
    {
        std::vector<bool> fullVect = std::vector<bool>(n, true);
        fullVect.reserve(max_n);
        return std::vector<std::vector<bool>>(1, fullVect);
    }
    else // else, the last site is either empty or filled. We find the combinations of the rest of the vector, of size n-1,
    {    // place the last hole or particle, and concatenate the two parts.
        std::vector<std::vector<bool>> leftHand = getCombinations(n - 1, p, max_n);
        for (std::vector<std::vector<bool>>::iterator it = leftHand.begin(); it != leftHand.end(); it++)
        {
            (*it).push_back(false);
        }

        std::vector<std::vector<bool>> rightHand = getCombinations(n - 1, p - 1, max_n);
        for (std::vector<std::vector<bool>>::iterator it = rightHand.begin(); it != rightHand.end(); it++)
        {
            (*it).push_back(true);
        }

        leftHand.insert(leftHand.end(), rightHand.begin(), rightHand.end());
        return leftHand;
    }

} // There are no check on the entrance variables. This programm will malfunction if there are not adapted, like if n<p.





/// @brief Check if two wave-function basis vectors are one particle tunneling away.
/// @param phiH1 wave function half basis bitset
/// @param phiH2 wave function half basis bitset
/// @param nSites number of sites (usefull !!!)
/// @param indexFrom If there is a difference, update it to be the last difference with a particle of phiH1
/// @param indexTo If there is a difference, update it to be the last difference with a particle of phiH2
/// @return True if exactly one particle can jump.
bool canHopHalf(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, uint nSitesHalf, uint& indexFrom, uint& indexTo)
{
    uint idFrom = 0;
    uint idTo = 0;
    uint numOfSim = countSimilarities(phiH1, phiH2, nSitesHalf, idFrom, idTo); // check if there are differences in the up-spin particle positions
    if ((nSitesHalf - numOfSim) == 2)                                            // if there is one difference, check if the hop is allowed
    {
        indexFrom = idFrom;
        indexTo = idTo;
        return true;
    }
    else
        return false;
}





/// @brief Convert a vector of boolean into a bitset of boolean 
/// @param vec vector of boolean to convert
/// @return the converted bitset 
std::bitset<MAX_BITS> bitsetFromVector(const std::vector<bool>& vec)
{
    std::bitset<MAX_BITS> bs = std::bitset<MAX_BITS>();
    uint maxV = std::min((int)vec.size(), MAX_BITS);
    for (uint it = 0; it<maxV; it++)
    {
        if (vec[it])  bs.set(it);
    }

    return bs;
}


/// @brief Swap two bitset 
/// @param bs1 bitset
/// @param bs2 bitset
/// @return void 
void bitsetswap(std::bitset<MAX_BITS>& bs1, std::bitset<MAX_BITS>& bs2)
{
   bs1 = bs1 ^ bs2;
   bs2 = bs2 ^ bs1;
   bs1 = bs1 ^ bs2;
}


/// @brief Sorting of a vector of bistset with the tri-fusion algorithm 
/// @param vec vector of particle positions, in bitset
/// @return void 
void fusionSortVector(std::vector<std::bitset<MAX_BITS>>& vec)
{
    uint vecSize = vec.size();

    if (vecSize == 1);

    else if (vecSize == 2)
    {
        if (vec[0].to_ulong() > vec[1].to_ulong())
        {
            bitsetswap(vec[0], vec[1]);
        }
        
    }

    else
    {

        uint mid = (uint)((vecSize / 2));

        std::vector<std::bitset<MAX_BITS>> vec1(vec.begin(), vec.begin() + mid);
        fusionSortVector(vec1);

        std::vector<std::bitset<MAX_BITS>> vec2(vec.begin() + mid + 1, vec.end());
        fusionSortVector(vec2);

        uint countor1 = 0;
        uint countor2 = 0;

        while (countor1 + countor2 < vecSize)
        {
            if (countor1 == vec1.size())
            {
                while (countor2 < vec2.size())
                {
                    vec[countor1 + countor2] = vec2[countor2];
                    countor2 ++;
                }
            }
            else if (countor2 == vec2.size())
            {
                while (countor1 < vec1.size())
                {
                    vec[countor1 + countor2] = vec1[countor1];
                    countor1 ++;
                }
            }
            else if (vec1[countor1].to_ulong() < vec2[countor2].to_ulong())
            {
                vec[countor1 + countor2] = vec1[countor1];
                countor1++;
            }
            else
            {
                vec[countor1 + countor2] = vec2[countor2];
                countor2++;
            }
        }

    }
}



/// @brief Convert a position vector into an empty/filled vector.
/// @param pos vector of particle positions
/// @return Boolean vector of empty/filled sites 
std::bitset<MAX_BITS> positionToBinary(std::vector<uint>& pos)
{
    std::bitset<MAX_BITS> bs = std::bitset<MAX_BITS>();
    for (std::vector<uint>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        bs.set(*it);
    }
    return bs;
}


/// @brief Convert an empty/filled vector into a position vector.
/// @param binaryPos a bitset of particle presence
/// @return Position vector 
std::vector<uint> binaryToPosition(std::bitset<MAX_BITS>& binaryPos)
{
    std::vector<uint> vec = std::vector<uint>();
    for (uint it = 0; it < binaryPos.size(); it++)
    {
        if (binaryPos[it])
        {
            vec.push_back(it);
        }
    }
    return vec;
}





/// @brief Count the similarities between 2 boolean half vectors and store the last difference.
/// @param phiH1 wave-function half basis bitset
/// @param phiH2 wave-function half basis bitset
/// @param lastIndexWrong1 If there is a difference, update it to be the last difference with a particle of phiH1
/// @param lastIndexWrong2 If there is a difference, update it to be the last difference with a particle of phiH2
/// @return Number of similarities between the 2 wave-function basis vector.
uint countSimilarities(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, uint nSites, uint& lastIndexWrong1, uint& lastIndexWrong2)
{
    uint simCount = 0;
    uint lastIdW1 = 0;
    uint lastIdW2 = 0;

    for (uint l = 0; l < nSites; l++)
    {
        if (phiH1[l] == phiH2[l])
            simCount += 1;
        else
        {
            if (phiH1[l] == true)
                    lastIdW1 = l;
            else if (phiH2[l] == true)
                    lastIdW2 = l;
        }
    }

    lastIndexWrong1 = lastIdW1;
    lastIndexWrong2 = lastIdW2;
    return simCount;
}


/// @brief Count the doublons in 2 bitsets.
/// @param phiH1 wave-function half basis bitset
/// @param phiH2 wave-function half basis bitset
/// @return Number of doublons
uint countDoublons(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2)
{
    std::bitset<MAX_BITS> andBits = phiH1 & phiH2;
    return andBits.count();
}


/// @brief Count the number of particles in a wave-function basis vector between 2 sites
/// @param posP empty/filled vector
/// @param from index of 1st site
/// @param to index of 2nd site
/// @return Number of particles between 'from' and 'to'
uint getParticleNumberBetween(std::bitset<MAX_BITS>& posP, uint from, uint to)
{
    uint nPart = 0;
    if ((from < MAX_BITS) && (to < MAX_BITS))
    {
        uint f = std::min(from, to);
        uint t = std::max(from, to);
        for (uint p = f + 1; p < t; p++)
        {
            if (posP[p])
                nPart++;
        }
    }

    return nPart;
}





std::bitset<MAX_BITS> expandBitSetUnderConstraint(std::bitset<MAX_BITS> reducedBS, std::bitset<MAX_BITS> constraintBS, uint sizeReducedBS)
{
    std::bitset<MAX_BITS> bs = std::bitset<MAX_BITS>();
    uint countorBS = 0;
    uint countorReducedBS = 0;
    while ((countorBS < MAX_BITS) && (countorReducedBS < sizeReducedBS))
    {
        if (constraintBS[countorBS])
        {
            countorBS ++;
        }
        else
        {
            bs[countorBS] = reducedBS[countorReducedBS];
            countorBS++;
            countorReducedBS++;
        }
    } 

    return bs;
}



/// @brief Copy of the EXPOKIT expmv function
/// @param t double
/// @param A sparse matrix
/// @param v vector
/// @param tol
/// @param m1
/// @return exp(t.A)v
arma::cx_dvec expmv(std::complex<double> t, arma::sp_cx_dmat& A, arma::cx_dvec& v, double tol, unsigned int m1 )
{
    unsigned int n = A.n_rows;
    unsigned int m = std::min(n, m1);
    double anorm = arma::norm(A, "inf");
    unsigned int mxrej = 10;
    double btol = 1.E-7;
    double gamma = 0.9;
    double delta = 1.2;
    unsigned int mb = m;
    unsigned int mx = mb;
    double t_out = std::abs(t);
    int nstep = 0;
    double t_new = 0;
    double t_now = 0;
    double t_step = 0;
    double s = 0;
    double s_error = 0;
    double rndoff = anorm * 1.E-15;

    unsigned int k1 = 2;
    double xm = 1/(double)m;
    double normv = arma::norm(v,2);
    double avnorm = 0;
    double phi1 = 0;
    double phi2 = 0;
    double err_loc = 0;
    double beta = normv;
    double fact = std::pow((m+1)/std::exp(1), m+1) * std::sqrt(2*M_PI*(m+1));
    t_new = (1/anorm) * pow((fact * tol) / 4*beta*anorm, xm);     
    std::complex<double> sgn = t/std::abs(t);
    nstep = 0;

    arma::cx_dvec w = v;
    double hump = normv;

    while(t_now < t_out)
    {
        nstep +=1;
        t_step = std::min(t_out - t_now, t_new);
        arma::cx_dmat V = arma::cx_dmat(n, m+1); 
        arma::cx_dmat H = arma::cx_dmat(m+2, m+2);
        arma::cx_dmat F;

        V.col(0) = (1/beta) * w;
        for (unsigned int j = 0; j<m; j++)
        {
            arma::cx_dvec p = A * V.col(j);
            for (unsigned int i = 0; i < j ; i++)
            {
                H(i,j) = arma::cdot(V.col(i), p);
                p = p - H(i,j)*V.col(i);
            }
            s = arma::norm(p,2);
            if (s<btol)
            {
                k1 = 0;
                mb = j;
                t_step = t_out - t_now;
                break;
            }
            H(j+1, j) = s;
            V.col(j+1) = (1/s) * p;
        }

        if (k1 != 0)
        {
            H(m+1, m) = 1;
            avnorm = arma::norm(A*V.col(m), 2);
        }
        unsigned int ireject = 0;
        while(ireject < mxrej)
        {
            mx = mb + k1;
            F = arma::expmat(sgn * t_step * H.submat(0, 0, mx-1, mx-1));
            if (k1 ==0)
            {
                err_loc = btol;
                break;
            }
            else
            {
                phi1 = std::abs(beta * F(m, 0));
                phi2 = std::abs(beta * F(m+1,0) * avnorm);
                if (phi1 > 10*phi2)
                {
                    err_loc = phi2;
                    xm = 1/(double)m;
                }
                else if (phi1 > phi2)
                {
                    err_loc = (phi1*phi2)/(phi1-phi2);
                    xm = 1/(double)m;
                }
                else
                {
                    err_loc = phi1;
                    xm = 1/(double)(m-1);
                }
            }
            if (err_loc <= delta * t_step * tol)
            {
                break;
            }
            else
            {
                t_step = gamma * t_step * pow((t_step * tol / err_loc), xm);
                s = pow(10, std::floor(std::log10(t_step))-1);
                t_step = std::ceil(t_step/s) * s;
                if (ireject == mxrej)
                {
                    printf("The requested tolerance is too high.");
                    // !!!! THROW AN ERROR HERE
                }
                ireject+=1;
            }
        }
        mx = mb + std::max((unsigned int)0, k1-1);
        w = V.cols(0, mx-1) * (beta * (F.col(0)).subvec(0,mx-1));
        beta = arma::norm(w, 2);
        hump = std::max(hump, beta);

        t_now = t_now + t_step;
        t_new = gamma * t_step * pow((t_step*tol)/err_loc, xm);
        s = pow(10, std::floor(std::log10(t_new))-1);
        t_new = std::ceil(t_new / s) * s;

        err_loc = std::max(err_loc, rndoff);
        s_error = s_error + err_loc;
    }

    return w;

}


