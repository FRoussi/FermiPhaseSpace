#include "FHEDSolver.h"

FHEDSolver::FHEDSolver(FHSystem *fhS)
{
    fhSystem = fhS;
}

void FHEDSolver::setTimeParameters(double t_max, double d_t, unsigned int nPerTS)
{
    this->tMax = t_max;
    this->dt = d_t;
    this->nStepsPerTimePoint = nPerTS;
    this->nTimePoints = floor(tMax / dt + 1);
}

void FHEDSolver::initiateSolver()
{
    this->posPartU = fhSystem->getInitialCondition_upSpin();
    this->posPartD = fhSystem->getInitialCondition_upSpin();
    this->nPartU = 0;
    for (std::vector<bool>::iterator itpu = posPartU.begin(); itpu != posPartU.end(); itpu++)
    {
        if (*itpu)
            nPartU +=1;
    }
    this->nPartD = 0;
    for (std::vector<bool>::iterator itpd = posPartD.begin(); itpd != posPartD.end(); itpd++)
    {
        if (*itpd)
            nPartD +=1;
    }
    this->basisSizeSpinUp = nChoosep(fhSystem->getnSites(), nPartU);
    this->basisSizeSpinDown = nChoosep(fhSystem->getnSites(), nPartD);
    this->basisSize = basisSizeSpinUp * basisSizeSpinDown;
    this->basisArray = std::vector<std::vector<unsigned int>>(2, std::vector<unsigned int>(basisSize, 0)); // 1 basis vector = 1 combination of up-spin position and 1 combination of down-spin position
    this->c_pu = getCombinations(this->fhSystem->getnSites(), this->nPartU, this->fhSystem->getnSites());
    this->c_pu.shrink_to_fit();
    this->c_pd = getCombinations(this->fhSystem->getnSites(), this->nPartD, this->fhSystem->getnSites());
    this->c_pd.shrink_to_fit();
    unsigned int countor = 0;
    for (uint l = 0; l < basisSizeSpinUp; l++)
    {
        for (uint k = 0; k < basisSizeSpinUp; k++)
        {
            basisArray[0][countor] = l; // up-spin particle combinations of position
            basisArray[1][countor] = k; // down-spin particle combinations of position
            countor += 1;
        }
    }

    unsigned int l = getBasisIndex(this->posPartU, this->posPartD);
    this->psi_0 = arma::cx_dvec(this->basisSize, arma::fill::zeros);
    this->psi_0[l] = 1;

    this->psi = this->psi_0;

    buildPairOperators();

    // Hamiltonian
    buildHamiltonian();

    // Operators :

    printf("Test 3 - Number operator :\n");
    NijuSet[0][0].print_dense();
}

void FHEDSolver::buildPairOperators()
{
    unsigned int nSites = fhSystem->getnSites();
    NijuSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);
    NijdSet = std::vector<std::vector<arma::sp_cx_dmat>>(nSites);
    for (unsigned int i = 0; i < nSites; i++)
    {
        NijuSet[i] = std::vector<arma::sp_cx_dmat>(nSites);
        NijuSet[i][i] = numberOperatorSpinUp(i);

        NijdSet[i] = std::vector<arma::sp_cx_dmat>(nSites);
        NijdSet[i][i] = numberOperatorSpinDown(i);
    }

    for (unsigned int i = 0; i < nSites; i++)
    {
        for (unsigned int j = i + 1; j < nSites; j++)
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

void FHEDSolver::buildHamiltonian_U()
{

    H_U = arma::sp_cx_dmat(basisSize, basisSize);
    for (unsigned int iSite = 0; iSite < fhSystem->getnSites(); iSite ++)
    {
        for (unsigned int jSite = 0 ; jSite < fhSystem->getnSites(); jSite ++)
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
    for (unsigned int l = 0; l < basisSize; l++)
    {
        values[l] = U_os_psil(l);
        locations(0, l) = l;
        locations(1, l) = l;
    }
    H_U = arma::sp_cx_dmat(locations, values, basisSize, basisSize);
    */
}

void FHEDSolver::buildHamiltonian_J()
{
    H_J = arma::sp_cx_dmat(basisSize, basisSize);
    for (unsigned int iSite = 0; iSite < fhSystem->getnSites(); iSite++)
    {
        for (unsigned int jSite = 0; jSite < fhSystem->getnSites(); jSite++)
        {
            if (fhSystem->getTunnelingEnergy(iSite, jSite) != arma::cx_double(0,0))
            {
                H_J = H_J + fhSystem->getTunnelingEnergy(iSite, jSite) * (NijuSet[iSite][jSite] + NijdSet[iSite][jSite]);
            }
        }
    }

}
/*
    // Up-spin particle hopps :
    unsigned int halhMatSizeUp = floor(pow(basisSizeSpinUp, 2) / 2);
    arma::umat locationsUp = arma::umat(2, halhMatSizeUp, arma::fill::zeros);
    arma::cx_drowvec valuesUp = arma::cx_drowvec(halhMatSizeUp, arma::fill::zeros);

    for (unsigned int l = 0; l < basisSizeSpinUp; l++)
    {
        unsigned int index = 0;
        for (unsigned int m = l + 1; m < basisSizeSpinUp; m++) // (l+1):upSpinBSize: # down spin is the same : m_d=l_d
        {                                                      // m0 = l_u*downSpinBSize + l_d = l + downSpinBSize
            unsigned int iSiteFrom = 0;
            unsigned int iSiteTo = 0;
            bool yn = canHopHalf(c_pu[l], c_pu[m], fhSystem->getnSites(), &iSiteFrom, &iSiteTo);
            if (yn)
            {
                locationsUp(0, index) = iSiteFrom;
                locationsUp(1, index) = iSiteTo;
                valuesUp(index) = fhSystem->getTunnelingEnergy(iSiteFrom, iSiteTo);
                index++;
            }
        }
    }

    arma::sp_cx_dmat Hopp_spinU = arma::sp_cx_dmat(true, locationsUp, valuesUp, basisSizeSpinUp, basisSizeSpinUp);

    // Down-spin particles hopps :
    unsigned int halhMatSizeDown = floor(pow(basisSizeSpinDown, 2) / 2);
    arma::umat locationsDown = arma::umat(2, halhMatSizeDown, arma::fill::zeros);
    arma::cx_drowvec valuesDown = arma::cx_drowvec(halhMatSizeDown, arma::fill::zeros);

    for (unsigned int l = 0; l < basisSizeSpinDown; l++)
    {
        unsigned int index = 0;
        for (unsigned int m = l + 1; m < basisSizeSpinDown; m++)
        {
            unsigned int iSiteFrom = 0;
            unsigned int iSiteTo = 0;
            bool yn = canHopHalf(c_pd[l], c_pd[m], fhSystem->getnSites(), &iSiteFrom, &iSiteTo);
            if (yn)
            {
                locationsDown(0, index) = iSiteFrom;
                locationsDown(1, index) = iSiteTo;
                valuesDown(index) = fhSystem->getTunnelingEnergy(iSiteFrom, iSiteTo);
                index++;
            }
        }
    }
    arma::sp_cx_dmat Hopp_spinD = arma::sp_cx_dmat(true, locationsDown, valuesDown, basisSizeSpinDown, basisSizeSpinDown);

    // Fusion of the two :
    unsigned int MatSize = floor(pow(basisSize, 2) / 2);
    arma::umat locations = arma::umat(2, MatSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(MatSize, arma::fill::zeros);

    unsigned index = 0;
    for (unsigned int l_part_u = 0; l_part_u < basisSizeSpinUp; l_part_u++)
    {
        unsigned int l_u = l_part_u * basisSizeSpinDown;
        for (unsigned int id = 0; id < basisSizeSpinUp; id++)
        {
            locations(0, index) = locationsUp(0, id) + l_u;
            locations(1, index) = locationsUp(1, id) + l_u;
            values(index) = valuesUp(id);
            index++;
        }
    }

    for (unsigned int l_part_d = 0; l_part_d < basisSizeSpinDown; l_part_d++)
    {
        for (unsigned int id = 0; id < basisSizeSpinDown; id++)
        {
            locations(0, index) = locationsDown(0, id) * basisSizeSpinDown + l_part_d;
            locations(1, index) = locationsDown(1, id) * basisSizeSpinDown + l_part_d;
            values(index) = valuesDown(id);
            index++;
        }
    }
    arma::sp_cx_dmat H_triU = arma::sp_cx_dmat(true, locations, values, basisSize, basisSize);


   H_J = H_triU + H_triU.st();
   */



void FHEDSolver::buildHamiltonian()
{
    buildHamiltonian_U();
    buildHamiltonian_J();

    H = H_U + H_J;

    printf("Test 2 - small system Hamiltonian : \n");
    H.print_dense();
    Hdense = arma::cx_dmat(H);
}

arma::sp_cx_dmat FHEDSolver::numberOperatorSpinUp(unsigned int iSite)
{
    arma::umat locations = arma::umat(2, basisSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(basisSize, arma::fill::zeros);
    unsigned int index = 0;

    if (nPartU != 0)
    {
        for (unsigned int lu = 0; lu < basisSizeSpinUp; lu++)
        {
            if (c_pu[lu][iSite])
            {
                for (unsigned int ld = 0; ld < basisSizeSpinDown; ld++)
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

arma::sp_cx_dmat FHEDSolver::numberOperatorSpinDown(unsigned int iSite)
{
    arma::umat locations = arma::umat(2, basisSize, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(basisSize, arma::fill::zeros);
    unsigned int index = 0;

    if (nPartD != 0)
    {
        for (unsigned int ld = 0; ld < basisSizeSpinDown; ld++)
        {
            if (c_pd[ld][iSite])
            {
                for (unsigned int lu = 0; lu < basisSizeSpinDown; lu++)
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

arma::sp_cx_dmat FHEDSolver::tunnelingOperatorSpinUp(unsigned int iSite, unsigned int jSite)
{
    unsigned int halfMatSizeUp = floor(basisSize * basisSizeSpinUp / 2);
    arma::umat locations = arma::umat(2, halfMatSizeUp, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(halfMatSizeUp, arma::fill::zeros);

    unsigned int index = 0;

    for (unsigned int lu = 0; lu < basisSizeSpinUp; lu++)
    {
        if (c_pu[lu][iSite] && (!c_pu[lu][jSite]))          // <- has the particle the possiblity to tunnel from iSite to jSite ?
        {
            std::vector<bool> c_pu_hoppdone = c_pu[lu];     // <- if yes, we write the new configuration in 'c_pu_hoppdone'
            c_pu_hoppdone[iSite] = false;
            c_pu_hoppdone[jSite] = true;
            for (unsigned int mu = 0; mu < basisSizeSpinUp; mu++)
            {
                if (c_pu_hoppdone == c_pu[mu])              // <-we look for the basis function index of the new configuration, mu
                {
                    for (unsigned int ld = 0; ld < basisSizeSpinDown; ld++) // <- we note it in the whole base, with the down spin particle part equals
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
arma::sp_cx_dmat FHEDSolver::tunnelingOperatorSpinDown(unsigned int iSite, unsigned int jSite)
{
     unsigned int halfMatSizeDown = floor(basisSize * basisSizeSpinDown / 2);
    arma::umat locations = arma::umat(2, halfMatSizeDown, arma::fill::zeros);
    arma::cx_drowvec values = arma::cx_drowvec(halfMatSizeDown, arma::fill::zeros);

    unsigned int index = 0;

    for (unsigned int ld = 0; ld < basisSizeSpinDown; ld++)
    {
        if (c_pd[ld][iSite] && (!c_pd[ld][jSite]))          // <- has the particle the possiblity to tunnel from iSite to jSite ?
        {
            std::vector<bool> c_pd_hoppdone = c_pd[ld];     // <- if yes, we write the new configuration in 'c_pd_hoppdone'
            c_pd_hoppdone[iSite] = false;
            c_pd_hoppdone[jSite] = true;
            for (unsigned int md = 0; md < basisSizeSpinDown; md++)
            {
                if (c_pd_hoppdone == c_pd[md])               // <-we look for the basis function index of the new configuration, md
                {
                    for (unsigned int lu = 0; lu < basisSizeSpinUp; lu++) // <- we note it in the whole base, with the up spin particle part equals
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

unsigned int FHEDSolver::getBasisIndex(std::vector<bool> posPU, std::vector<bool> posPD)
{
    unsigned int lu = 0;
    bool index1Found = false;
    while (!index1Found && lu < this->basisSizeSpinUp)
    {
        if (posPU == c_pu[lu])
            index1Found = true;
        else
            lu++;
    }
    unsigned int ld = 0;
    bool index2Found = false;
    while (!index2Found && ld < this->basisSizeSpinDown)
    {
        if (posPD == c_pd[ld])
            index2Found = true;
        else
            ld++;
    }
    if (index1Found && index2Found)
        return lu * this->basisSizeSpinDown + ld;
    else
    {
        printf("The basis vector of a particle configuration cannot be found. \n");
        return 0;
    }
}

// Contribution of the interaction of opposite-spin electrons of basis
std::complex<double> FHEDSolver::U_os_psil(unsigned int phi_l)
{
    std::vector<bool> *posPU = &c_pu[basisArray[0][phi_l]];
    std::vector<bool> *posPD = &c_pd[basisArray[1][phi_l]];

    std::complex<double> u = 0;
    for (unsigned int i_u = 0; i_u < fhSystem->getnSites(); i_u++)
    {
        if ((*posPU)[i_u])
        {
            for (unsigned int i_d = 0; i_d < fhSystem->getnSites(); i_d++)
            {
                if ((*posPD)[i_d])
                    u = u + fhSystem->getInteractionEnergy(i_u, i_d);
            }
        }
    }
    return u;
}

void FHEDSolver::step()
{
    std::complex<double> alpha = std::complex<double>(0, -1);
    // this->psi =  this->psi + alpha *  this->psi;

    // arma::cx_dmat a = arma::expmat(this->Hdense);
    //this->psi = arma::expmat(alpha * Hdense) * this->psi;

    this->psi = expmv(this->dt, alpha * H, this->psi);
}

void FHEDSolver::runSimulation()
{

    std::vector<double> timeLine = std::vector<double>(nTimePoints, 0);
    arma::dmat niiu = arma::dmat(fhSystem->getnSites(), nTimePoints, arma::fill::zeros);
    arma::drowvec EU = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec EJ = arma::drowvec(nTimePoints, arma::fill::zeros);
    arma::drowvec E = arma::drowvec(nTimePoints, arma::fill::zeros);

    this->psi = this->psi_0;

    for (unsigned int t = 0; t < nTimePoints; t++)
    {
        timeLine[t] = t * dt;
        for (unsigned int iSite = 0; iSite < fhSystem->getnSites(); iSite++)
        {
            niiu(iSite,t) = static_cast<double>(std::real(arma::cdot(psi, NijuSet[iSite][iSite] * psi)));
        }     
        EU[t] = static_cast<double>(std::real(arma::cdot(psi, H_U * psi)));
        EJ[t] = static_cast<double>(std::real(arma::cdot(psi, H_J * psi)));
        E[t] = EU[t] + EJ[t];

        step();
    }

    printf("Test 3 - small system dynamic, n_1u : \n");
    niiu.row(0).print();
    

    // Writing the results :
    std::string folderName = fhSystem->getFolder();
    std::string fileName = fhSystem->getName();
    std::string fEDName = fileName.substr(0, fileName.size() - 4);
    std::string fEDAdressName = fEDName + "_ED";

    ResultWriter results = ResultWriter(fEDAdressName, folderName);
    results.newResults();
    results.writeLine("t", timeLine);
    for (unsigned int iSite=0; iSite<fhSystem->getnSites(); iSite++)
    {
        results.writeLine("N_" + std::to_string(iSite) + "_" + std::to_string(iSite) + "_u", niiu.row(iSite));
    }
    results.writeLine("E_U", EU);
    results.writeLine("E_J", EJ);
    results.writeLine("E", E);

    

}

FHEDSolver::~FHEDSolver()
{
}

// ---------------------------
// Other functions

unsigned int nChoosep(unsigned int n, unsigned int p)
{
    unsigned int num = 1;
    unsigned int den = 1;
    if (p < (n - p))
    {
        for (unsigned int k = 1; k <= p; k++)
            den *= k;

        for (unsigned int k = n - p + 1; k <= n; k++)
            num *= k;
    }
    else
    {
        for (unsigned int k = 1; k <= n - p; k++)
            den *= k;

        for (unsigned int k = p + 1; k <= n; k++)
            num *= k;
    }

    return (unsigned int)num / den;
}

std::vector<bool> positionToBinary(std::vector<unsigned int> pos, unsigned int nSites)
{
    std::vector<bool> vec = std::vector<bool>(nSites, false);
    for (std::vector<unsigned int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        vec[*it] = true;
    }
    return vec;
}

std::vector<unsigned int> binaryToPosition(std::vector<bool> binaryPos)
{
    std::vector<unsigned int> vec = std::vector<unsigned int>();
    for (unsigned int it = 0; it < binaryPos.size(); it++)
    {
        if (binaryPos[it])
        {
            vec.push_back(it);
        }
    }
    return vec;
}

// getCombinations(n, p, max_n) builds a vector of boolean vectors of all combination of p particles in n sites.
// max_n is the final size of the vectors, the total number of sites. It does not crash if it's wrong, it will just be less efficient.
std::vector<std::vector<bool>> getCombinations(unsigned int n, unsigned int p, unsigned int max_n)
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

unsigned int countSimilarities(std::vector<bool> phiH1, std::vector<bool> phiH2, unsigned int *lastIndexWrong1, unsigned int *lastIndexWrong2)
{
    unsigned int simCount = 0;
    unsigned int lastIdW1 = 0;
    unsigned int lastIdW2 = 0;
    if (phiH1.size() == phiH1.size())
    {
        for (unsigned int l = 0; l < phiH1.size(); l++)
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
    }

    *lastIndexWrong1 = lastIdW1;
    *lastIndexWrong2 = lastIdW2;
    return simCount;
}

bool canHopHalf(std::vector<bool> phiH1, std::vector<bool> phiH2, unsigned int nSites, unsigned int *indexFrom, unsigned int *indexTo)
{
    unsigned int idFrom = 0;
    unsigned int idTo = 0;
    unsigned int numOfSim = countSimilarities(phiH1, phiH2, &idFrom, &idTo); // check if there are differences in the up-spin particle positions
    if ((nSites - numOfSim) == 2)                                            // if there is one difference, check if the hop is allowed
    {
        *indexFrom = idFrom;
        *indexTo = idTo;
        return true;
    }
    else
        return false;
}

unsigned int getParticleNumberBetween(std::vector<bool> posP, unsigned int from, unsigned int to)
{
    unsigned int nPart = 0;
    if ((from < posP.size()) && (to < posP.size()))
    {
        unsigned int f = std::min(from, to);
        unsigned int t = std::max(from, to);
        for (unsigned int p = f + 1; p < t; p++)
        {
            if (posP[p])
                nPart++;
        }
    }

    return nPart;
}


arma::cx_dvec expmv(double t, arma::sp_cx_dmat A, arma::cx_dvec v, double tol, unsigned int m1 )
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
    double sgn = t/std::abs(t);
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

        V.col(1) = (1/beta) * w;
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
            H(m+2, m+1) = 1;
            avnorm = arma::norm(A*V.col(m+1), 2);
        }
        unsigned int ireject = 0;
        while(ireject < mxrej)
        {
            mx = mb + k1;
            F = arma::expmat(sgn * t_step * H.submat(0, mx-1, 0, mx-1));
            if (k1 ==0)
            {
                err_loc = btol;
                break;
            }
            else
            {
                phi1 = std::abs(beta * F(m+1, 1));
                phi2 = std::abs(beta * F(m+2,1) * avnorm);
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
