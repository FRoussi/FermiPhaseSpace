#include "tJEDSolver.h"

tJEDSolver::tJEDSolver(FHSystem *fhS)
{
    fhSystem = fhS;
    writeResults = false;
}

void tJEDSolver::setResultFile(std::string wN, std::string fN)
{
    workName = wN;
    folderName = fN;
    writeResults = true;
}


void tJEDSolver::setTimeParameters(double t_max, double d_t, uint nPerTS)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nPerTS;
    this->nTimePoints = floor(tMax / dt + 1);
}


/// @brief Initiate the ED solver, compute the initial wave-function, then the usefull operators.
///@return true if it's a successful initiation
bool tJEDSolver::initiateSolver()
{
    bool succcess = false;

    this->posPartU = bitsetFromVector(fhSystem->getInitialCondition_upSpin());
    this->posPartD = bitsetFromVector(fhSystem->getInitialCondition_downSpin());
    this->nPartU = posPartU.count();
    this->nPartD = posPartD.count();

    if (nPartU + nPartD >= fhSystem->getnSites())
    {
        std::cout << "Initialisation problem : too many particles for a tJ model" << std::endl;
        return;
    }
    if (countDoublons(posPartU, posPartD) > 0)
    {
        std::cout << "Initialisation problem : the inital position have a doublon." << std::endl;
        return;  
    }


    this->basisSizeSpinUp = nChoosep(fhSystem->getnSites(), nPartU);
    this->basisSizeSpinDown = nChoosep(fhSystem->getnSites() - nPartU, nPartD);
    this->basisSize = basisSizeSpinUp * basisSizeSpinDown;
    this->basisArray = std::vector<std::vector<uint>>(2, std::vector<uint>(basisSize, 0)); // 1 basis vector = 1 combination of up-spin position and 1 combination of down-spin position
    
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
        for (uint k = 0; k < basisSizeSpinUp; k++)
        {
            basisArray[0][countor] = l;         //<- up-spin particle combinations of position
            basisArray[1][countor] = k;         //<- down-spin particle combinations of position
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
void tJEDSolver::buildPairOperators()
{
    uint nSites = fhSystem->getnSites();
    NijuSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);
    NijdSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);
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
            if (fhSystem->areNeighbour(i, j))
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
void tJEDSolver::buildHamiltonian_U()
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
/*
    arma::umat locations = arma::umat(2, basisSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(basisSize, arma::fill::zeros);
    for (uint l = 0; l < basisSize; l++)
    {
        values[l] = U_os_psil(l);
        locations(0, l) = l;
        locations(1, l) = l;
    }
    H_U = arma::sp_cx_dmat(locations, values, basisSize, basisSize);
    */
}


///@brief Compute the system's Hamiltonian tunneling part
///@return void
void tJEDSolver::buildHamiltonian_tM()
{
    H_t = arma::sp_cx_dmat(basisSize, basisSize);
    for (uint iSite = 0; iSite < fhSystem->getnSites(); iSite++)
    {
        for (uint jSite = 0; jSite < fhSystem->getnSites(); jSite++)
        {
            if (fhSystem->getTunnelingEnergy(iSite, jSite) != arma::cx_double(0,0))
            {
                H_t = H_t + fhSystem->getTunnelingEnergy(iSite, jSite) * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]);
            }
        }
    }

}


///@brief Compute the system's Hamiltonian
///@return void
void tJEDSolver::buildHamiltonian()
{
    buildHamiltonian_U();
    buildHamiltonian_tM();

    H = H_U + H_t;

    //printf("Test 2 - small system Hamiltonian : \n");
    //H.print_dense();
    //Hdense = arma::cx_dmat(H);
}


///@brief Compute the operator "$\\hat{n}_{ii\\uparrow}$"
///@param iSite Site number
///@return The operator "$\\hat{n}_{ii\\uparrow}$"
arma::sp_cx_dmat tJEDSolver::numberOperatorSpinUp(uint iSite)
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
arma::sp_cx_dmat tJEDSolver::numberOperatorSpinDown(uint iSite)
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
                for (uint lu = 0; lu < basisSizeSpinDown; lu++)
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
arma::sp_cx_dmat tJEDSolver::tunnelingOperatorSpinUp(uint iSite, uint jSite)
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
arma::sp_cx_dmat tJEDSolver::tunnelingOperatorSpinDown(uint iSite, uint jSite)
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
uint tJEDSolver::getBasisIndex(std::bitset<MAX_BITS>& posPU, std::bitset<MAX_BITS>& posPD)
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
std::complex<double> tJEDSolver::U_os_psil(uint phi_l)
{
    std::bitset<MAX_BITS> *posPU = &c_pu[basisArray[0][phi_l]];
    std::bitset<MAX_BITS> *posPD = &c_pd[basisArray[1][phi_l]];

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
void tJEDSolver::step(double d_t)
{
    //std::complex<double> alpha = std::complex<double>(0, -1);

    //this->psi = expmv(this->dt, alpha * H, this->psi);

    std::complex<double> alpha = std::complex<double>(0, -d_t); // <- alpha = -i dt / h

    this->psi = expmv(alpha, H, this->psi);
}


/// @brief Launch simulation and write the results
void tJEDSolver::runSimulation()
{

std::cout << "Start simulation" << std::endl << std::flush;



    std::vector<double> timeLine = std::vector<double>(nTimePoints, 0);
    arma::dmat niiu = arma::dmat(fhSystem->getnSites(), nTimePoints, arma::fill::zeros);
    arma::dmat niid = arma::dmat(fhSystem->getnSites(), nTimePoints, arma::fill::zeros);
    arma::drowvec EU = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec EtH = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec E = arma::drowvec(nTimePoints, arma::fill::zeros);

    this->psi = this->psi_0;

    // Bad filling thing :
    uint tt = nTimePoints/20;
    std::cout << "|                     |" << std::endl;
    std::cout << "|" ;

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

        step(this->dt);

        if ((t%tt) == 0) {
            std::cout << "=" << std::flush ;        }

    }
std::cout << "|" << std::endl; ;
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
        results.writeLine("EU", EU);
        results.writeLine("EtH", EtH);
        results.writeLine("E", E);
    }


    

}


/// @brief Destructor
tJEDSolver::~tJEDSolver()
{
}

// ---------------------------
// Other functions


/// @brief Binomial coefficient or all combinations
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
/// @param phiH1 wave function half basis vector
/// @param phiH2 wave function half basis vector
/// @param nSites number of sites (usefull ???)
/// @param indexFrom If there is a difference, update it to be the last difference with a particle of phiH1
/// @param indexTo If there is a difference, update it to be the last difference with a particle of phiH2
/// @return True if exactly one particle can jump.
bool canHopHalf(std::bitset<MAX_BITS>& phiH1, std::bitset<MAX_BITS>& phiH2, uint nSites, uint& indexFrom, uint& indexTo)
{
    uint idFrom = 0;
    uint idTo = 0;
    uint numOfSim = countSimilarities(phiH1, phiH2, nSites, idFrom, idTo); // check if there are differences in the up-spin particle positions
    if ((nSites - numOfSim) == 2)                                            // if there is one difference, check if the hop is allowed
    {
        indexFrom = idFrom;
        indexTo = idTo;
        return true;
    }
    else
        return false;
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
