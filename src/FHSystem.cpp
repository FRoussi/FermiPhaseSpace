#include "../include/FHSystem.h"


FHSystem::FHSystem(bool talk)
{

    talkative = talk;
    lattice = FHLattice(talk);

    uint nSites = lattice.getnSites();

    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>(nSites);
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>(nSites);

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;

    TBGroundStateComputed = false;
    MFGroundStateComputed = false;

    setInitialConditionsDiagonal(initCond);
    initiateNullInteraction();

  //  tJExCType               = tJExchangeCouplingType::tJJnone;
   // externalPotentialType   = ExternalPotentialType::ExtPotNone;

    
}


FHSystem::FHSystem(int x, bool talk)
{
    talkative = talk;

lattice = FHLattice(x, talk);

    arma::dmat pos = arma::dmat(1,  x, arma::fill::zeros);
    for (int d=0; d<x; d++)
    {
        pos[d] = (double)d;
    }

    int nPart = std::floor(x/2);
    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>(nPart,0);
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>(nPart,0);
    for(int itPart = 0; itPart <nPart; itPart++)
    {
        initCond_up[itPart] = itPart * 2;
        initCond_down[itPart] = itPart * 2 + 1;
    }

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;


    setPositions(pos);
    setInitialConditionsDiagonal(initCond);

    TBGroundStateComputed = false;
    MFGroundStateComputed = false;

    initiateNullInteraction();

  //  UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
  //  tunnelingType = FHTunnelingType::nearestNeighboursConstant;
  //  VInteractionType        = FHSameSpinInteractionType::Vnone;
  //  tJExCType               = tJExchangeCouplingType::tJJnone;
  //  externalPotentialType   = ExternalPotentialType::ExtPotNone;

}


FHSystem::FHSystem(int x, int y, int z, bool talk)
{

    talkative = talk;

    lattice = FHLattice(x,y,z,talk);

    arma::dmat pos = arma::dmat(3,  x * y * z, arma::fill::zeros);
    int tabIndex = 0;
    for (int i_z = 0; i_z < z; i_z++)
    {
        for (int i_y = 0; i_y < y; i_y++)
        {
            for (int i_x = 0; i_x < x; i_x++)
            {   
                pos(0,tabIndex) = double(i_x);
                pos(1,tabIndex) = double(i_y);
                pos(2,tabIndex) = double(i_z); 
                tabIndex +=1;
            }
        }       
    }

    int nPart = std::floor( x * y * z/2);
    std::vector<unsigned int> initCond_up =  std::vector<unsigned int>(nPart,0);
    std::vector<unsigned int> initCond_down =  std::vector<unsigned int>(nPart,0);
    for(int itPart = 0; itPart <nPart; itPart++)
    {
        initCond_up[itPart] = itPart * 2;
        initCond_down[itPart] = itPart * 2 + 1;
    }

    std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
    initCond[0] = initCond_up;
    initCond[1] = initCond_down;

    setPositions(pos);
    externalPotential = FHExternalPotential(&lattice);
    setInitialConditionsDiagonal(initCond);

    TBGroundStateComputed = false;
    MFGroundStateComputed = false;

    initiateNullInteraction();
  //  UInteractionType        = FHOppositeSpinInteractionType::UonSiteConstant;
  //  tunnelingType           = FHTunnelingType::nearestNeighboursConstant;
  //  VInteractionType        = FHSameSpinInteractionType::Vnone;
  //  tJExCType               = tJExchangeCouplingType::tJJnone;
  //  externalPotentialType   = ExternalPotentialType::ExtPotNone;

}



FHSystem::FHSystem(Json::Value rootSystem, std::string folderName, bool talk)
{

    talkative = talk;

    Json::Value rootLattice = rootSystem.get("lattice","" );

    lattice = FHLattice(rootLattice, folderName, talk);
    

    TBGroundStateComputed = false;
    MFGroundStateComputed = false;

    initiateNullInteraction();

// Tunneling Energy
    Json::Value tunnJson = rootSystem["tunnelingInteraction"];
    tunnelingInteraction = FHTunnelingInteraction(&lattice, tunnJson);


// Coulomb interaction
    Json::Value coulJson = rootSystem["coulombInteraction"];
    std::string coulType = coulJson.get("type", "null").asString();
    coulombInteraction = FHCoulombInteraction(&lattice, coulJson);

// External Potential
    Json::Value ExtPotJ = rootSystem["externalPotential"];
    if (ExtPotJ == Json::nullValue)
    {
        externalPotential = FHExternalPotential(&lattice);
    }
    else
    {
        std::string extPotType = ExtPotJ.get("type", "null").asString();
        if (extPotType == "null")
        {
            externalPotential = FHExternalPotential(&lattice);
        }
        else if (extPotType == "constant")
        {
            externalPotential = FHExternalPotentialConstant(&lattice, arma::dvec(lattice.getnSites()));
        }
        else if (extPotType == "coreHole")
        {
            uint site = ExtPotJ.get("site", 0).asUInt();
            double hPot = ExtPotJ.get("holePotential", -1).asDouble();
            externalPotential = FHExternalPotentialCoreHole(&lattice, site, hPot);
        }
        
    }


// Initial condition
    Json::Value jsonIC = rootSystem["IC"];
    if (jsonIC != Json::Value::null)
    {
        std::string typeIC = jsonIC.get("type", "givenDiagonal").asString();
        if (typeIC == "givenDiagonal")
        {
            initialCondition = FHInitialConditionType::givenDiagonal;
            std::string subtypeIC = jsonIC.get("subtype", "spinWave").asString();
            initConditions_upSpin = std::vector<bool>(getnSites(), false);
            initConditions_downSpin = std::vector<bool>(getnSites(), false);

            if (subtypeIC == "spinWave")
            {
               setInitialConditions(FHInitialConditionType::spinWave);
            }  
            else if (subtypeIC == "chargeWave")
            {
                setInitialConditions(FHInitialConditionType::chargeWave);
            }     
            else if (subtypeIC == "spinWave3")
            {

                uint midParticleNumber = std::ceil(getnSites()/2);
                upSpinParticleNumber = 0;
                downSpinParticleNumber = 0;
                for (uint iPart = 0; iPart < midParticleNumber ; iPart++ )
                {
                    initConditions_upSpin[iPart] = true;
                    upSpinParticleNumber++;
                }
                for (uint iPart = midParticleNumber; iPart < getnSites() ; iPart++ )
                {
                    initConditions_downSpin[iPart] = true;
                    downSpinParticleNumber++;
                }

            }  
            else if (subtypeIC == "chargeWave3")
            {
                uint midParticleNumber = std::ceil(getnSites()/2);
                upSpinParticleNumber = 0;
                downSpinParticleNumber = 0;
                for (uint iPart = 0; iPart < midParticleNumber ; iPart++ )
                {
                    initConditions_upSpin[iPart] = true;
                    initConditions_downSpin[iPart] = true;
                    upSpinParticleNumber++;
                    downSpinParticleNumber++;
                }

            }         
            else if (subtypeIC == "givenInJson")
            {
                Json::Value ICNu = jsonIC["IC_NuDiag"];
                if (ICNu.isArray())
                {
                    for (uint iSite = 0; iSite < std::min((uint)ICNu.size(), (uint)getnSites()); iSite++)
                    {
                        initConditions_upSpin[iSite] = ICNu[iSite].asBool();
                    }
                }
                Json::Value ICNd = jsonIC["IC_NdDiag"];
                if (ICNd.isArray())
                {
                    for (uint iSite = 0; iSite < std::min((uint)ICNd.size(), (uint)getnSites()); iSite++)
                    {
                        initConditions_downSpin[iSite] = ICNd[iSite].asBool();
                    }
                }
            }   
            else if (subtypeIC == "givenInFile")
            {
                std::string initConditions_fileName = jsonIC.get("IC_fileName", "").asString();
                if (initConditions_fileName == "")
                {
                    std::cout << "Initial condition file not written in Json entry." << std::endl;
                }
                else 
                {
                    std::cout << "Reading Initial condition in " << initConditions_fileName << "." << std::endl;
                }

            }        

        }
        else if (typeIC == "tbGroundState")
        {
            //std::cout << "  -- Initial Condition : TB ground State." << std::endl;
            setInitialConditions(FHInitialConditionType::tbGroundState);
        }
        else if (typeIC == "mfGroundState")
        {
            setInitialConditions(FHInitialConditionType::mfGroundState);
        }
        else if (typeIC == "groundState")
        {
            setInitialConditions(FHInitialConditionType::groundState);
        }
        else if (typeIC == "givenInFile")
        {

        }
        else
        {

        }

        std::cout << "   -- Initial condtion :            " << typeIC << std::endl;  
        std::cout << "   -- Particle number :             " << upSpinParticleNumber + downSpinParticleNumber << std::endl;  
    }
    else
    {
        std::cout << "No Initial condition in the json file.  " << std::endl;  
    }



}


FHSystem::FHSystem(FHLattice lat)
{
    talkative = true;

    lattice = lat;
    TBGroundStateComputed = false;
    MFGroundStateComputed = false;
    
    initiateNullInteraction();
}



FHLattice                                   FHSystem::getLattice()                                  {    return lattice;}
unsigned int                                FHSystem::getnSites()                                   {    return lattice.getnSites();}
unsigned int                                FHSystem::getnDimension()                               {    return lattice.getnDimension();}
arma::drowvec                               FHSystem::getSitePosition(uint iSite)                   {    return lattice.getSitePosition(iSite);}
double                                      FHSystem::getCoordinate(uint dim, uint iSite)           {    return lattice.getCoordinate(dim, iSite);}
arma::drowvec                               FHSystem::getLatticeSitePosition(uint iSite)            {    return lattice.getLatticeSitePosition(iSite);}
double                                      FHSystem::getLatticeCoordinate(uint dim, uint iSite)    {    return lattice.getLatticeCoordinate(dim, iSite);}

bool                                        FHSystem::sitesHaveNames()                              {    return lattice.sitesHaveNames();}
std::string                                 FHSystem::getSiteName(uint iSite)                       {    return lattice.getSiteName(iSite);}
bool                                        FHSystem::isPeriodic()                                  {    return lattice.isPeriodic();}

arma::dmat                                  FHSystem::getLatticeMatrix()                            {    return lattice.getLatticeMatrix();}
arma::drowvec                               FHSystem::getLatticeZeroVector()                        {    return lattice.getLatticeZeroVector();}

bool                                        FHSystem::isSystemTimeDependant()                       {    return externalPotential.isThePotentialTimeDependant();}
bool                                        FHSystem::isExternalPotentialNull()                     {    return externalPotential.isThePotentialNull();}
bool                                        FHSystem::isTBGroundStateComputed()                     {    return TBGroundStateComputed;}
bool                                        FHSystem::isMFGroundStateComputed()                     {    return MFGroundStateComputed;}

double                                      FHSystem::getDistance(uint i, uint j)                   {    return lattice.getDistance(i,j);}
bool                                        FHSystem::areNeighbour(uint i, uint j)                  {    return lattice.areNeighbour(i,j);}

FHSystem::FHInitialConditionType            FHSystem::getInitialCondition()                         {    return initialCondition;}
std::vector<bool>                           FHSystem::getInitialCondition_upSpin()                  {    return initConditions_upSpin;}
std::vector<bool>                           FHSystem::getInitialCondition_downSpin()                {    return initConditions_downSpin;}
unsigned int                                FHSystem::getUpSpinParticleNumber()                     {    return upSpinParticleNumber;}
unsigned int                                FHSystem::getDownSpinParticleNumber()                   {    return downSpinParticleNumber;}

FHCoulombInteractionType                    FHSystem::getCoulombInteractionType()                   {    return coulombInteraction.getCoulombInteractionType();}
FHTunnelingInteractionType                  FHSystem::getTunnelingInteractionType()                 {    return tunnelingInteraction.getTunnelingInteractionType();}
FHExternalPotentialType                     FHSystem::getExternalPotentialType()                    {    return externalPotential.getExternalPotentialType();}
//FHSystem::FHSameSpinInteractionType         FHSystem::getSameSpinInteractionType()                  {    return VInteractionType;}
//FHSystem::tJExchangeCouplingType            FHSystem::getExchangeCouplingEnergyType()               {    return tJExCType;}
std::complex<double>                        FHSystem::getCoulombInteraction(int i, int j)            {    return coulombInteraction.getCoulombUInteraction(i,j);}
std::complex<double>                        FHSystem::getSameSpinInteractionEnergy(int i, int j)    {    return coulombInteraction.getCoulombVInteraction(i,j);}
std::complex<double>                        FHSystem::getTunnelingInteraction(int i, int j)              {    return tunnelingInteraction.getTunnelingInteraction(i,j);}
//std::complex<double>                        FHSystem::getExchangeCouplingEnergy(int i, int j)       {    return JMatrix(i,j);}

std::complex<double>                        FHSystem::getExternalPotential(int i, double t)         {    return externalPotential.getExternalPotential(t,i,i);}
std::complex<double>                        FHSystem::getExternalPotential(int i, int j, double t)  {    return externalPotential.getExternalPotential(t,i,j);}
arma::cx_dmat                               FHSystem::getExternalPotential(double time)             {    return externalPotential.getExternalPotential(time);}
void                                        FHSystem::updateExternalPotential(double time)          {    externalPotential.updatePotential(time);}

const arma::cx_dmat*                        FHSystem::getUMatrixPointer()                           {   return coulombInteraction.getCoulombUInteractionPointer();}
const arma::cx_dmat*                        FHSystem::getVMatrixPointer()                           {   return coulombInteraction.getCoulombVInteractionPointer();}
//const arma::cx_dmat*                        FHSystem::getJMatrixPointer()                           {   return &JMatrix;}
const arma::cx_dmat*                        FHSystem::gettMatrixPointer()                           {   return tunnelingInteraction.getTunnelingInteractionPointer();}
const arma::cx_dmat*                        FHSystem::getExtPotPointer()                            {   return externalPotential.getExternalPotentialPointer();}
const arma::cx_dmat*                        FHSystem::get_tbgs_Nu_Pointer()                         {   return &tbgs_Nu;}
const arma::cx_dmat*                        FHSystem::get_tbgs_Nd_Pointer()                         {   return &tbgs_Nd;}
const arma::cx_dmat*                        FHSystem::get_tbgs_V_Pointer()                          {   return &tbgs_V;}
const arma::dvec*                           FHSystem::get_tbgs_E_Pointer()                          {   return &tbgs_E;}
const arma::cx_dmat*                        FHSystem::get_mfgs_Nu_Pointer()                         {   return &mfgs_Nu;}
const arma::cx_dmat*                        FHSystem::get_mfgs_Nd_Pointer()                         {   return &mfgs_Nd;}
const arma::cx_dmat*                        FHSystem::get_mfgs_Vu_Pointer()                         {   return &mfgs_Vu;}
const arma::cx_dmat*                        FHSystem::get_mfgs_Vd_Pointer()                         {   return &mfgs_Vd;}
const arma::dvec*                           FHSystem::get_mfgs_Eu_Pointer()                         {   return &mfgs_Eu;}
const arma::dvec*                           FHSystem::get_mfgs_Ed_Pointer()                         {   return &mfgs_Ed;}

std::string                                 FHSystem::getName()                                     { return systemName;}
std::string                                 FHSystem::getFolder()                                   { return systemFolder;}
std::string                                 FHSystem::getFPSVersion()                               { return std::to_string(FPS_VERSION_MAJOR) + "." + std::to_string(FPS_VERSION_MINOR) + "." + std::to_string(FPS_VERSION_PATCH);}



// SETTERS

    void                                 FHSystem::setLattice(FHLattice newLat)
    {
        lattice = newLat;
    }

void                    FHSystem::setPositions(const arma::dmat& positions, std::string sitesName)
{
    lattice.setPositions(positions, sitesName);
}

void                    FHSystem::setPositions(const arma::dmat& positions, std::vector<std::string> sN)
{
    lattice.setPositions(positions, sN);
}


bool                    FHSystem::setSiteTypes(arma::uvec sTypes)
{
    return lattice.setSiteTypes(sTypes);
}

void                    FHSystem::setSiteTypeNames(std::vector<std::string> typeNames)
{
    lattice.setSiteTypeNames(typeNames);
}


void                    FHSystem::setSymmetry(arma::dmat symVectors)
{
     lattice.setSymmetry(symVectors);
}

void                                FHSystem::setLatticeZero(arma::drowvec latticeZ)
{
    lattice.setLatticeZero(latticeZ);
}


void                                FHSystem::setNeighbourDistance(double neigbD)
{
    lattice.setNeighbourDistance(neigbD);
}


void                                FHSystem::setSamePointDistance(double samePointD)
{
    lattice.setSamePointDistance(samePointD);
}




/// @brief Set up the Fermi-Hubbard system from site positions and vectors of site number as initial condition
/// @param diagonalInitConditions vectors of sites on which there is a up-spin/down-spin particle
/// @return void
void FHSystem::setInitialConditionsDiagonal(const std::vector<std::vector<unsigned int>>& diagonalInitConditions )
{

    initialCondition =  FHSystem::FHInitialConditionType::givenDiagonal;

    upSpinParticleNumber = diagonalInitConditions[0].size();
    downSpinParticleNumber = diagonalInitConditions[1].size();

    initConditions_upSpin = std::vector<bool>(getnSites(), false);
    initConditions_downSpin = std::vector<bool>(getnSites(), false);



    for (unsigned int itSite = 0; itSite < diagonalInitConditions[0].size(); itSite++)
    {
        initConditions_upSpin[diagonalInitConditions[0][itSite]] = true;
    }
    for (unsigned int itSite = 0; itSite < diagonalInitConditions[1].size(); itSite++)
    {
        initConditions_downSpin[diagonalInitConditions[1][itSite]] = true;
    }

}


/// @brief Set up the Fermi-Hubbard system from site positions the number of particles as initial condition
/// @param initCondType a kind of initial condition (tbGroundState, charge-wave, spin-wave ...)
/// @return void
void FHSystem::setInitialConditions(FHInitialConditionType initCondType)
{

    initialCondition = initCondType;

    initConditions_upSpin = std::vector<bool>(getnSites(), false);
    initConditions_downSpin = std::vector<bool>(getnSites(), false);
    upSpinParticleNumber = 0;
    downSpinParticleNumber = 0 ; 

    switch (initialCondition)
    {
        case FHInitialConditionType::spinWave:
        {
            for (uint iSite = 0; iSite < getnSites(); iSite ++)
            {
                //std::cout <<  "      ----- part filling : iSite=" << iSite << std::endl;
                if (iSite%2 == 0)
                {
                    initConditions_upSpin[iSite] = true;
                    upSpinParticleNumber++;
                  //  std::cout <<  "      ----- part filling : up part" << iSite << std::endl;
                }
                else
                {
                    initConditions_downSpin[iSite] = true;
                    downSpinParticleNumber++;
                    //std::cout <<  "      ----- part filling : down part" << iSite << std::endl;
                }

            }
            break;
        }   
        case FHInitialConditionType::chargeWave:
        {
            for (uint iSite = 0; iSite < getnSites(); iSite ++)
            {
                if (iSite%2 == 0)
                {
                    initConditions_upSpin[iSite] = true;
                    upSpinParticleNumber++;
                    initConditions_downSpin[iSite] = true;
                    downSpinParticleNumber++;
                }
            }
            break;
        }
        case FHInitialConditionType::tbGroundState:
        {
            //std::cout << "  -- Initial Condition : TB ground State." << std::endl;
            upSpinParticleNumber = std::ceil((double)lattice.getnSites() / 2.); 
            downSpinParticleNumber = std::floor((double)lattice.getnSites() / 2.); 
            computeTBGroundState();
            break;
        }
        case FHInitialConditionType::mfGroundState:
        {
            upSpinParticleNumber = std::ceil((double)lattice.getnSites() / 2.); 
            downSpinParticleNumber = std::floor((double)lattice.getnSites() / 2.); 
            computeMFGroundState();
            break;
        }
        case FHInitialConditionType::groundState:
        {
            upSpinParticleNumber = std::ceil((double)lattice.getnSites() / 2.); 
            downSpinParticleNumber = std::floor((double)lattice.getnSites() / 2.); 
            computeMFGroundState();
            break;
        }
        default:
        {
            upSpinParticleNumber = getnSites()/2;
            downSpinParticleNumber = getnSites()/2;
            break;
        }
    }

}

/// @brief Set up the Fermi-Hubbard system from site positions the number of particles as initial condition
/// @param initCondType a kind of initial condition (tbGroundState, charge-wave, spin-wave ...)
/// @param partNumner vector of the numbers of up-spin, down-spin particle
/// @return void
void FHSystem::setInitialConditions(FHInitialConditionType initCondType, std::vector<unsigned int> partNumber )
{

    initialCondition = initCondType;
    initConditions_upSpin = std::vector<bool>(getnSites(), false);
    initConditions_downSpin = std::vector<bool>(getnSites(), false);

    upSpinParticleNumber = partNumber[0];
    downSpinParticleNumber = partNumber[1];
    switch(initialCondition)
    {
        case FHInitialConditionType::givenDiagonal:
        {
            for (uint iSite = 0; iSite < std::min((uint)getnSites(), upSpinParticleNumber); iSite +=2)
            {
                    initConditions_upSpin[iSite] = true;
            }
            for (uint iSite = 1; iSite < std::min((uint)getnSites(), downSpinParticleNumber); iSite +=2)
            {
                    initConditions_downSpin[iSite] = true;
            }
            break;
        }
        case FHInitialConditionType::spinWave:
        {
            for (uint iSite = 0; iSite < std::min((uint)getnSites(), upSpinParticleNumber); iSite +=2)
            {
                    initConditions_upSpin[iSite] = true;
            }
            for (uint iSite = 0; iSite < std::min((uint)getnSites(), downSpinParticleNumber); iSite +=2)
            {
                    initConditions_downSpin[iSite] = true;
            }
            break;
        }
        case FHInitialConditionType::tbGroundState:
        {
            computeTBGroundState();
            break;
        }
        case FHInitialConditionType::mfGroundState:
        {
            computeMFGroundState();
            break;
        }
    }
}



void                               FHSystem::setBasisInteractions()
{
    tunnelingInteraction = FHTunnelingInteraction(&lattice, 1., -1.);
    coulombInteraction = FHCoulombInteraction(&lattice, 1.);
}


/*
/// @brief setter - on-site interaction energy
/// @param U0 on-site interaction energy
void                    FHSystem::setUInteraction_OnSite(double U0, bool add)
{
    if (add)
    {
        UMatrix += U0 * arma::cx_dmat(getnSites(), getnSites(), arma::fill::eye);
    }
    else
    {
        UMatrix = U0 * arma::cx_dmat(getnSites(), getnSites(), arma::fill::eye);
    }

}

/// @brief setter - On-site interaction energy
/// @param U Vector of on-site interaction energy
void                    FHSystem::setUInteraction_GivenDiagonal(arma::dvec& Ua, bool add)
{
    if (Ua.n_rows == getnSites()) 
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(arma::diagmat(Ua), arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(arma::diagmat(Ua), arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }      
    }
    else
    {
        std::cout << "Wrong size for the Ua interaction vector." << std::endl;
    }
}

/// @brief setter - Matrix interaction energy
/// @param Uab Matrix interaction energy
void                    FHSystem::setUInteraction_GivenMatrix(arma::dmat& Uab, bool add)
{
    if ((Uab.n_rows == getnSites()) && (Uab.n_cols == getnSites()))
    {
        if (add)
        {
            UMatrix += arma::cx_dmat(Uab, arma::dmat(Uab.n_rows, Uab.n_cols, arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dmat(Uab, arma::dmat(Uab.n_rows, Uab.n_cols, arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the Uab interaction matrix." << std::endl;
    }

}


/// @brief setter - Matrix interaction energy
/// @param Us Matrix interaction energy by site type
void                    FHSystem::setUInteraction_BySiteType(arma::dvec& Us, bool add)
{
    if (Us.n_rows == lattice.getSiteNameNumber()) 
    {
        arma::dvec Uvec = arma::dvec(getnSites(), arma::fill::zeros);
        for (uint iSite = 0; iSite < getnSites(); iSite++)
        {
            Uvec[iSite] = Us[lattice.getSiteType(iSite)];
        }
        if (add)
        {
            UMatrix += arma::cx_dvec(arma::diagmat(Uvec), arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
        else
        {
            UMatrix = arma::cx_dvec(arma::diagmat(Uvec), arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the Us interaction matrix." << std::endl;
    }
}

/// @brief setter - nearest neighbout tunneling interaction energy
/// @param tH1 tunneling interaction energy
void                    FHSystem::setTunneling_BetweenNeighbours(double tH1, bool add)
{

    if (!add)
    {
        tMatrix = arma::cx_dmat(getnSites(), getnSites(), arma::fill::zeros);
    }
    for (unsigned int i = 0; i< getnSites(); i++)
    {
        for (unsigned int j = 0; j< getnSites(); j++)
        {
            if (areNeighbour(i,j))
            {
                tMatrix(i,j) += tH1;
            }
        }
    }
    //printf("J-Matrix : \n");
    //tMatrix.print();
}

/// @brief setter - Same site tunneling interaction energy
/// @param tH0 tunneling interaction energy
void                    FHSystem::setTunneling_SamePosition(double tH0, bool add)
{
    if (!add)
    {
        tMatrix = arma::cx_dmat(getnSites(), getnSites(), arma::fill::zeros);
    }
    for (unsigned int i = 0; i< getnSites(); i++)
    {
        for (unsigned int j = 0; j< getnSites(); j++)
        {
            if (lattice.areSameSite(i,j) && (i!=j))
            {
                tMatrix(i,j) += tH0;
            }
        }
    }
}

/// @brief setter - Tunneling interaction energy on sites at distance 'dist'
/// @param tHd tunneling interaction energy
/// @param dist distance between two sites
void                    FHSystem::setTunneling_atDistance(double tHd, double dist, bool add)
{
    if (!add)
    {
        tMatrix = arma::cx_dmat(getnSites(), getnSites(), arma::fill::zeros);
    }
    for (unsigned int i = 0; i< getnSites(); i++)
    {
        for (unsigned int j = 0; j< getnSites(); j++)
        {
            if (abs(getDistance(i,j) - dist) < 0.0001)
            {
                tMatrix(i,j) += tHd;
            }
        }
    }
}

/// @brief setter - matrix tunneling interaction energy
/// @param tH tunneling interaction energy by pair of sites
void                    FHSystem::setTunneling_GivenMatrix(arma::dmat tH, bool add)
{
    if ((tH.n_rows == getnSites()) && (tH.n_cols == getnSites()))
    {
        if (add)
        {
            tMatrix += arma::cx_dmat(tH, arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
        else
        {
            tMatrix = arma::cx_dmat(tH, arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
    }
    else
    {
        std::cout << "Wrong size for the tunneling matrix." << std::endl;
    }
}

/// @brief setter - matrix tunneling interaction energy
/// @param tHs1s2 tunneling interaction energy by pair of site types
void                    FHSystem::setTunneling_BySiteType(arma::dmat& tHs1s2, double d, bool add)
{
    if ((tHs1s2.n_rows == lattice.getSiteNameNumber()) && (tHs1s2.n_cols == lattice.getSiteNameNumber()))
    {
          arma::dvec tHMat = arma::dmat(getnSites(), getnSites(), arma::fill::zeros);
        for (uint iSite = 0; iSite < getnSites(); iSite++)
        {
            for (uint jSite = 0; jSite < getnSites(); jSite++)
            {
                if (std::abs(getDistance(iSite, jSite) - d) < 0.0001)
                {
                    tHMat(iSite, jSite) = tHs1s2(lattice.getSiteType(iSite), lattice.getSiteType(jSite));
                }
            }
        }
        if (add)
        {
            tMatrix += arma::cx_dvec(tHMat, arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
        else
        {
            tMatrix = arma::cx_dvec(tHMat, arma::dmat(getnSites(), getnSites(), arma::fill::zeros));
        }
    }
}

void                    FHSystem::setExchangeCoupling(double J, bool add)
{
    JMatrix = arma::cx_dmat(getnSites(), getnSites(), arma::fill::zeros);
}

*/

void                    FHSystem::setName(std::string name)    
{
    systemName = name;
}


void                    FHSystem::makePotential_pulse()
{
    double intensityA0 = 1;
    double centert0 = 1.;
    double widthtd = 0.1;
    double frequencyw0 = 100;
    arma::dvec polar = arma::dvec(getnDimension(), arma::fill::zeros);
    polar(0) = 1/std::sqrt(2);
    polar(1) = 1/std::sqrt(2);
    externalPotential = FHExternalPotentialPulse(&lattice, intensityA0, centert0, widthtd, frequencyw0, polar, tunnelingInteraction.getTunnelingInteractionPointer());
}



/*
/// @brief Read a potential stored in fileName.
/// @return void
void                    FHSystem::readPotential(std::string fileName, std::string folderName)
{

    std::string fileAdressName = fileName; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + fileName;
    }



    std::ifstream systFile(fileAdressName);
    if (systFile.is_open())
    {
        std::string line;
        std::string coord;
        std::getline(systFile, line);
        
        if (line == "#EXTERNALPOTENTIAL_NULL")                                                //<- ready to read the positions
        {
            externalPotentialType = ExternalPotentialType::ExtPotNone;
        }
        else if (line == "#EXTERNALPOTENTIAL_CONSTANT")   
        {
            externalPotentialType = ExternalPotentialType::ExtPotConstant;
        }
        else if (line == "#EXTERNALPOTENTIAL_PULSE")
        {
            externalPotentialType = ExternalPotentialType::ExtPotPulse;
            std::getline(systFile, line);
            double intensityA0 = std::stod(line);
            std::getline(systFile, line);
            double centert0 = std::stod(line);
            double widthtd = std::stod(line);
            std::getline(systFile, line);
            double frequencyw0 = std::stod(line);
            std::getline(systFile, line);
            arma::dvec polar = arma::dvec(getnDimension(), arma::fill::zeros);
            std::stringstream ss(line);
            uint coordCountor = 0;
            while (std::getline(ss, coord, ' ')) 
            {
                if (coordCountor < getnDimension())
                {
                    polar(coordCountor) = std::stod(coord);
                }
                coordCountor++;
            }
            externalPotential = FHExternalPotentialPulse(&lattice, intensityA0, centert0, widthtd, frequencyw0, polar, &tMatrix);
        } 
    }
    else
    {
        std::cout << "External Potential " << fileName << " cannot be read" << std::endl;
    }
}

*/

void                    FHSystem::initiateNullInteraction()
{

   // tJExCType = tJExchangeCouplingType::tJJnone;

    //uint nSites = getnSites();

    //UMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    //VMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
   // JMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    //tMatrix = arma::cx_dmat(nSites, nSites, arma::fill::zeros);

    coulombInteraction = FHCoulombInteraction(&lattice);
    tunnelingInteraction = FHTunnelingInteraction(&lattice);
    externalPotential = FHExternalPotential(&lattice);
}

void                                FHSystem::computeTBGroundState(bool withExternalPotential)
{
arma::cx_dmat* tMatrix = tunnelingInteraction.getTunnelingInteractionPointer();

    if (withExternalPotential)
    {
        arma::eig_sym(tbgs_E, tbgs_V, *tMatrix + *externalPotential.getExternalPotentialPointer());
    }
    else
    {
        arma::eig_sym(tbgs_E, tbgs_V, *tMatrix);
    }

//std::cout << "     --- Comp TB 1, npu=" << upSpinParticleNumber << std::endl;
    arma::cx_dvec diagMku = arma::cx_dvec(lattice.getnSites(), arma::fill::zeros);
    diagMku.subvec(0,upSpinParticleNumber-1).ones();
    arma::cx_dmat Mku = arma::diagmat(diagMku);

    arma::cx_dvec diagMkd = arma::cx_dvec(lattice.getnSites(), arma::fill::zeros);
    diagMkd.subvec(0,downSpinParticleNumber-1).ones();
    arma::cx_dmat Mkd = arma::diagmat(diagMkd);

//std::cout << "     --- Comp TB 2" << std::endl;
    tbgs_Nu = tbgs_V * Mku * tbgs_V.t();
    tbgs_Nd = tbgs_V * Mkd * tbgs_V.t();
//std::cout << "     --- Comp TB 3" << std::endl;
    TBGroundStateComputed  = true;


}


void                                FHSystem::computeMFGroundState(bool withExternalPotential)
{

arma::cx_dmat* tMatrix = tunnelingInteraction.getTunnelingInteractionPointer();
arma::cx_dmat* UMatrix = coulombInteraction.getCoulombUInteractionPointer();

    uint nSites = lattice.getnSites();

    arma::cx_dvec diagMku = arma::cx_dvec(lattice.getnSites(), arma::fill::zeros);
    diagMku.subvec(0,upSpinParticleNumber-1).ones();
    arma::cx_dmat Mku = arma::diagmat(diagMku);

    arma::cx_dvec diagMkd = arma::cx_dvec(lattice.getnSites(), arma::fill::zeros);
    diagMkd.subvec(0,downSpinParticleNumber-1).ones();
    arma::cx_dmat Mkd = arma::diagmat(diagMkd);

    double thresholdConverged = 0.0001 * nSites;
    bool isConverged = false;
    int maxIt = 100;
    iterationsForMFState = 0;

    arma::cx_dvec diagNu = arma::cx_dvec(nSites) * (double)upSpinParticleNumber / (double)nSites;
    arma::cx_dvec diagNd = arma::cx_dvec(nSites) * (double)downSpinParticleNumber / (double)nSites;

    arma::cx_dmat mfuHuU = arma::cx_dmat(nSites, nSites, arma::fill::zeros);
    arma::cx_dmat mfuHuD = arma::cx_dmat(nSites, nSites, arma::fill::zeros);


    while ((iterationsForMFState < maxIt) && !isConverged)
    {

        for (int i = 0; i < nSites; i++)
        {
            mfuHuU.at(i,i) = arma::dot(UMatrix->row(i), diagNd);
            mfuHuD.at(i,i) = arma::dot(UMatrix->row(i), diagNu);
        } 

        if (withExternalPotential)
        {
            arma::eig_sym(mfgs_Eu, mfgs_Vu, *tMatrix + mfuHuU);
            arma::eig_sym(mfgs_Ed, mfgs_Vd, *tMatrix + mfuHuD);
        }
        else
        {
            arma::eig_sym(mfgs_Eu, mfgs_Vu, *tMatrix + mfuHuU);
            arma::eig_sym(mfgs_Ed, mfgs_Vd, *tMatrix + mfuHuD);
        }
 
        mfgs_Nu = mfgs_Vu * Mku * mfgs_Vu.t();
        mfgs_Nd = mfgs_Vd * Mkd * mfgs_Vd.t();

        iterationsForMFState+=1;
        residueForMFState = arma::norm(mfgs_Nu.diag() - diagNu) + arma::norm(mfgs_Nd.diag() - diagNd);

        if (residueForMFState < thresholdConverged)
        {
            isConverged = true;
        }
        else
        {
            //double r = (double)iterationsForMFState / (double)maxIt;
            //diagNu = ((1-r)*diagNu + r*mfgs_Nu.diag())/2;
            //diagNd = ((1-r)*diagNd + r*mfgs_Nd.diag())/2;
            diagNu = mfgs_Nu.diag();
            diagNd = mfgs_Nd.diag();
        }

    }

    if (iterationsForMFState == maxIt)
    {
        std::cout << "!! Carefull : the MF ground state computation has not converged. Residue=" << residueForMFState << std::endl;
    }

    MFGroundStateComputed = true;

}



/// @brief Write a JSON version of the Fermi-Hubbard System
/// @param folderName emplacement of the folder
Json::Value                    FHSystem::writeFHSystemJson()
{


    Json::Value systemRoot;

    // Lattice 
    systemRoot["lattice"] =  lattice.writeFHLatticeJson();



//Initial condition
    Json::Value jsonIC = Json::Value(Json::objectValue);
    switch (initialCondition)
    {
    case FHInitialConditionType::spinWave:
        jsonIC["type"] = "givenDiagonal";
        jsonIC["subtype"] = "spinWave";
        break;
    case FHInitialConditionType::chargeWave:
        jsonIC["type"] = "givenDiagonal";
        jsonIC["subtype"] = "chargeWave";
        break;
    case FHInitialConditionType::tbGroundState:
        jsonIC["type"] = "tbGroundState";
        break;
    case FHInitialConditionType::groundState:
        jsonIC["type"] = "groundState";
        break;
    case FHInitialConditionType::givenDiagonal:
        jsonIC["type"] = "givenDiagonal"; 
        break;
    case FHInitialConditionType::mfGroundState:
        jsonIC["type"] = "mfGroundState"; 
        break;
    default:
        break;
    } 
    systemRoot["IC"] = jsonIC;
// Energy
    systemRoot["tunnelingInteraction"] = tunnelingInteraction.writeTunnelingInteractionJson();
    systemRoot["coulombInteraction"] = coulombInteraction.writeCoulombInteractionJson();
    /*
    switch (tunnelingType)
    {
    case FHTunnelingType::tHnone:
        systemRoot["EnergyT"] = "None";
        break;   
    case FHTunnelingType::nearestNeighboursConstant:
    {
        systemRoot["EnergyT"] = "nearestNeighbours";
        arma::uvec q1 = find(tMatrix,1);
        if (q1.n_rows !=0)
        {
            double th = std::real( tMatrix(q1(0)) );
            systemRoot["parametersT"] = th;
        }
        else
        {
            systemRoot["parametersT"] = 0;
        }
        break;
    }
    default:
        break;
    }

    switch (UInteractionType)
    {
    case FHOppositeSpinInteractionType::Unone:
        systemRoot["EnergyU"] = "None";
        break;   
    case FHOppositeSpinInteractionType::UonSiteConstant:
    {
        systemRoot["EnergyU"] = "onSite";
        systemRoot["parametersU"] = std::real( UMatrix(0,0));
        break;
    }
    default:
        break;
    }
    // V interaction
    switch (VInteractionType)
    {
    case FHSameSpinInteractionType::Vnone:
        systemRoot["EnergyV"] = "None";
        break;   
    case FHSameSpinInteractionType::VonSiteConstant:
    {
        systemRoot["EnergyV"] = "onSite";
        systemRoot["parametersV"] = std::real( VMatrix(0,0));
        break;
    }
    default:
        break;
    }
*/
        // External potential
    Json::Value ExtPotJ = externalPotential.writePotentialJson();
    systemRoot["externalPotential"] = ExtPotJ;


    // TB Ground State
    if (isTBGroundStateComputed())
    {
        systemRoot["tbGroundState"] = writeTBGroundStateJson();
    }

    // MF Ground State
    if (isMFGroundStateComputed())
    {
        systemRoot["mfGroundState"] = writeMFGroundStateJson();
    }

    return systemRoot;
}


Json::Value                         FHSystem::writeTBGroundStateJson()
{
    if (!isTBGroundStateComputed())
    {
        computeTBGroundState();
    }

    Json::Value tbgsJson(Json::objectValue);
    Json::Value jsonVecE(Json::arrayValue);
    const arma::dvec* vecE = get_tbgs_E_Pointer();
    for (uint i = 0; i < vecE->size(); i++)
    {
        jsonVecE.append(vecE->at(i));
    }  
    tbgsJson["energyVector"] = jsonVecE;
    Json::Value jsonMatEig(Json::arrayValue);
    const arma::cx_dmat* matV = get_tbgs_V_Pointer();
    for (uint i = 0; i < matV->n_rows; i++)
    {
        Json::Value jsonMatEigJ(Json::arrayValue);
        for (uint j = 0; j < matV->n_cols; j++)
        {
            jsonMatEigJ.append(matV->at(i,j).real());
        }
        jsonMatEig.append(jsonMatEigJ);
    }  
    tbgsJson["rotationMatrix"] = jsonMatEig;
    return tbgsJson;
}


Json::Value                         FHSystem::writeMFGroundStateJson()
{

    if (!isMFGroundStateComputed())
    {
        computeMFGroundState();
    }


    Json::Value mfgsJson(Json::objectValue);

    Json::Value jsonVecEu(Json::arrayValue);
    const arma::dvec* vecEu = get_mfgs_Eu_Pointer();
    for (uint i = 0; i < vecEu->size(); i++)
    {
        jsonVecEu.append(vecEu->at(i));
    }  
    mfgsJson["energyVectorUp"] = jsonVecEu;

    Json::Value jsonVecEd(Json::arrayValue);
    const arma::dvec* vecEd = get_mfgs_Ed_Pointer();
    for (uint i = 0; i < vecEd->size(); i++)
    {
        jsonVecEd.append(vecEd->at(i));
    }  
    mfgsJson["energyVectorDown"] = jsonVecEd;

    Json::Value jsonMatEigU(Json::arrayValue);
    const arma::cx_dmat* matV = get_mfgs_Vu_Pointer();
    for (uint i = 0; i < matV->n_rows; i++)
    {
        Json::Value jsonMatEigJ(Json::arrayValue);
        for (uint j = 0; j < matV->n_cols; j++)
        {
            jsonMatEigJ.append(matV->at(i,j).real());
        }
        jsonMatEigU.append(jsonMatEigJ);
    }  
    mfgsJson["rotationMatrixUp"] = jsonMatEigU;

    Json::Value jsonMatEigD(Json::arrayValue);
    matV = get_mfgs_Vd_Pointer();
    for (uint i = 0; i < matV->n_rows; i++)
    {
        Json::Value jsonMatEigJ(Json::arrayValue);
        for (uint j = 0; j < matV->n_cols; j++)
        {
            jsonMatEigJ.append(matV->at(i,j).real());
        }
        jsonMatEigD.append(jsonMatEigJ);
    }  
    mfgsJson["rotationMatrixDown"] = jsonMatEigD;

    mfgsJson["iterationsNeeded"] = iterationsForMFState;
    mfgsJson["residueComputation"] = residueForMFState;

    return mfgsJson;
}

Json::Value                    FHSystem::getOriginalJsonSystem()
{
    return originalJsonSystem;
}




 void FHSystem::writeWannierToolsDatFile(std::string wtDatFile, std::string folderName)
 {
    std::string fileAdressName = wtDatFile; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + wtDatFile;
    }

    uint nSites = getnSites();

    std::ofstream datFile(fileAdressName);
    if (datFile.is_open())
    {

        datFile << "  ! HmnR file from sparse hr file" << std::endl; 
        datFile << "\t\t" << nSites <<  "\t! Num_wann: number of wannier orbitals" << std::endl;       
        datFile << "\t\t" << 9 <<  "\t! NRPTS: number of R points" << std::endl;  
        datFile << "    1    1    1    1    1    1    1    1    1" << std::endl; 

        //std::cout << std::fixed;
        for (uint iSite = 0; iSite < nSites; iSite ++ )
        {
            for (uint jSite = 0; jSite < nSites; jSite ++ )
            {
                datFile << "    0    0    0    " << iSite+1 << "    " << jSite+1 << "    " <<  std::to_string(tunnelingInteraction.getTunnelingInteraction(iSite, jSite).real()) << "    " << std::to_string((double)0) << std::endl; 
            }
        }

        datFile.close();
    }

 }

 void FHSystem::writeWannierToolsInBandFile(std::string wtIBFile, std::string wtDatFile, std::string folderName)
 {
    std::string fileAdressName = wtIBFile; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + wtIBFile;
    }

    //uint nSites = getnSites();

    std::ofstream IBStream(fileAdressName);
    if (IBStream.is_open())
    {
        IBStream << "&TB_FILE" << std::endl; 
        IBStream << " Hrfile = \'" << wtDatFile <<  "\'" << std::endl << "/" << std::endl << std::endl ;   

        IBStream << "!> bulk band structure calculation flag" << std::endl << "&CONTROL" << std::endl << "BulkBand_calc    = T" << std::endl << "/" << std::endl << std::endl ;

        IBStream << " &SYSTEM" << std::endl;
        IBStream << "NumOccupied =  " << 10 << "     \t! NumOccupied" << std::endl;
        IBStream << "SOC = " << 0 << "               \t! without" << std::endl;
        IBStream << "E_FERMI =  0.00     \t! e-fermi" << std::endl;
        IBStream << "/" << std::endl << std::endl;

        IBStream << "&PARAMETERS" << std::endl;
        IBStream << "E_arc = 0.00        \t! obtained NumSelectedEigenVals bands close to E_arc" << std::endl;
        IBStream << "Nk1 = " << 41 << "          \t! number k points" << std::endl;
        IBStream << "/" << std::endl << std::endl;

        WannierToolsSystemInofStream(IBStream);
        
        IBStream.close();
    }
 }

 void FHSystem::writeWannierToolsDosFile(std::string wtDosFile, std::string wtDatFile, std::string folderName)
 {
    std::string fileAdressName = wtDosFile; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + wtDosFile;
    }

    //uint nSites = getnSites();

    std::ofstream DOSStream(fileAdressName);
    if (DOSStream.is_open())
    {
        DOSStream << "&TB_FILE" << std::endl; 
        DOSStream << " Hrfile = \'" << wtDatFile <<  "\'" << std::endl << "/" << std::endl << std::endl ;   

        DOSStream << "!> bulk band structure calculation flag" << std::endl << "&CONTROL" << std::endl << "DOS_calc    = T" << std::endl << "/" << std::endl << std::endl ;

        DOSStream << " &SYSTEM" << std::endl;
        DOSStream << "NumOccupied =  " << 10 << "      \t! NumOccupied" << std::endl;
        DOSStream << "SOC = " << 0 << "                \t! without" << std::endl;
        DOSStream << "E_FERMI =  0.00        \t! e-fermi" << std::endl;
        DOSStream << "/" << std::endl;

        DOSStream << "&PARAMETERS" << std::endl;
        DOSStream << "Eta_arc = "     << std::to_string(0.05)     << "   \t! infinite small value, like brodening" << std::endl;
        DOSStream << "E_arc = "       << std::to_string(0.00)     << "   \t! obtained NumSelectedEigenVals bands close to E_arc" << std::endl;
        DOSStream << "OmegaNum = "    << std::to_string(401)      << "   \t! omega number" << std::endl;
        DOSStream << "OmegaMin = "    << std::to_string(-0.3)     << "   \t! energy interval" << std::endl;
        DOSStream << "OmegaMax = "    << std::to_string(0.3)      << "   \t! energy interval" << std::endl;
        DOSStream << "Nk1 = "         << std::to_string(41)       << "   \t! number k points" << std::endl;
        DOSStream << "Nk2 = "         << std::to_string(41)       << "   \t! number k points" << std::endl;
        DOSStream << "/" << std::endl << std::endl;

        WannierToolsSystemInofStream(DOSStream);
        
        DOSStream.close();
    }
 }

 void FHSystem::writeWannierToolsInBandUnfoldFile(std::string wtIBUFile, std::string wtDatFile, std::string folderName)
 {
        std::string fileAdressName = wtIBUFile; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + wtIBUFile;
    }

   // uint nSites = getnSites();

    std::ofstream IBUStream(fileAdressName);
    if (IBUStream.is_open())
    {
        IBUStream << "&TB_FILE" << std::endl; 
        IBUStream << " Hrfile = \'" << wtDatFile <<  "\'" << std::endl << "/" << std::endl << std::endl ;   

        IBUStream << "!> bulk band structure calculation flag" << std::endl << "&CONTROL" << std::endl;
        IBUStream << "BulkBand_calc    = T" << std::endl;
        IBUStream << "BulkFS_plane_calc    = F" << std::endl;
        IBUStream << "BulkBand_Unfold_line_calc    = T" << std::endl;
        IBUStream << "BulkBand_Unfold_plane_calc    = F" << std::endl;
        IBUStream << "/" << std::endl << std::endl ;

        IBUStream << " &SYSTEM" << std::endl;
        IBUStream << "NumOccupied =  " << 10 << "     \t! NumOccupied" << std::endl;
        IBUStream << "SOC = " << 0 << "               \t! without" << std::endl;
        IBUStream << "E_FERMI =  0.00       \t! e-fermi" << std::endl;
        IBUStream << "/" << std::endl;

        IBUStream << "&PARAMETERS" << std::endl;
        IBUStream << "Eta_arc = "           << std::to_string(0.05)     << "   \t! infinite small value, like brodening" << std::endl;
        IBUStream << "E_arc = "             << std::to_string(0.00)     << "   \t! obtained NumSelectedEigenVals bands close to E_arc" << std::endl;
        IBUStream << "OmegaNum_unfold = "   << std::to_string(600)      << "   \t! omega number" << std::endl;
        IBUStream << "OmegaMin = "          << std::to_string(-13)      << "   \t! energy interval" << std::endl;
        IBUStream << "OmegaMax = "          << std::to_string(7)        << "   \t! energy interval" << std::endl;
        IBUStream << "Nk1 = "               << std::to_string(41)       << "   \t! number k points" << std::endl;
        IBUStream << "Nk2 = "               << std::to_string(41)       << "   \t! number k points" << std::endl;
        IBUStream << "/" << std::endl << std::endl;

        WannierToolsSystemInofStream(IBUStream);
        
        IBUStream.close();
    }
 }

 void FHSystem::WannierToolsSystemInofStream(std::ofstream& wannierFileStream)
 {
        uint nSites = getnSites();
        uint nProjectors_unfold = 2;


        wannierFileStream << "LATTICE" << std::endl << "Angstrom" << std::endl;
        for (uint i = 0; i < 3; i++)
        {
            for (uint j = 0; j < 3; j++)
            {
                 wannierFileStream << "       " << std::to_string(getLatticeMatrix()(i,j));
            }
            wannierFileStream << std::endl;
        }
        wannierFileStream << std::endl;

        wannierFileStream << "ATOM_POSITIONS" << std::endl << "           " << nSites << std::endl;
        wannierFileStream << "Direct" << std::endl;
        for (uint iSite = 0; iSite < nSites; iSite++ )
        {

            arma::drowvec latPos = getLatticeSitePosition(iSite);
            wannierFileStream << getSiteName(iSite);
            for (uint j = 0; j<3; j++)
            {
                wannierFileStream << "      " << std::to_string(latPos(j));
            } 
            wannierFileStream << std::endl;
        }

        wannierFileStream << "PROJECTORS" << std::endl << "           " << nSites << "*1" << std::endl;
        for (uint iSite = 0; iSite < nSites; iSite++ )
        {
            wannierFileStream << getSiteName(iSite) << " pz" << std::endl;;
        }

        wannierFileStream << std::endl << std::endl ;

        wannierFileStream << "LATTICE_UNFOLD" << std::endl << "Angstrom" << std::endl;
        for (uint i = 0; i < 3; i++)
        {
            for (uint j = 0; j < 3; j++)
            {
                 wannierFileStream << "       " << std::to_string(getLatticeMatrix()(i,j));
            }
            wannierFileStream << std::endl;
        }
        wannierFileStream << std::endl;

        wannierFileStream << "ATOM_POSITION_UNFOLD" << std::endl << nProjectors_unfold << std::endl << "Direct" << std::endl;
        for (uint iSite = 0; iSite < nProjectors_unfold; iSite++ )
        {
            arma::drowvec latPos = getLatticeSitePosition(iSite);
            wannierFileStream << getSiteName(iSite);
            for (uint j = 0; j<3; j++)
            {
                wannierFileStream << "      " << std::to_string(latPos(j));
            } 
            wannierFileStream << std::endl;
        }
        wannierFileStream << std::endl;

        wannierFileStream << "PROJECTORS_UNFOLD" << std::endl << nProjectors_unfold << "*1" << "         ! number of projectors" << std::endl;
        for (uint iSite = 0; iSite < nProjectors_unfold; iSite++ )
        {
            wannierFileStream << getSiteName(iSite) << " pz" << std::endl;;
        }
        wannierFileStream << std::endl;

        wannierFileStream << "SELECTED_ATOMS" << std::endl;
        wannierFileStream << 1 << std::endl;
        wannierFileStream << "1-2" << std::endl;
        wannierFileStream << std::endl;


        wannierFileStream << "SURFACE" << std::endl;
        wannierFileStream << "  "<< 1 << " " << 0 << " " << 0 << std::endl;
        wannierFileStream << "  "<< 0 << " " << 1 << " " << 0 << std::endl;
        wannierFileStream << "  "<< 0 << " " << 0 << " " << 1 << std::endl;
        wannierFileStream << std::endl;

        wannierFileStream << "KPATH_BULK            ! k point path" << std::endl;
        wannierFileStream << "3              ! number of k line only for bulk band" << std::endl;
        wannierFileStream << "   M   0.00000  0.50000  0.00000   K   0.33333  0.66667  0.00000" << std::endl;
        wannierFileStream << "   K   0.33333  0.66667  0.00000   G   0.00000  0.00000  0.00000" << std::endl;
        wannierFileStream << "   G   0.00000  0.00000  0.00000   M   0.00000  0.50000  0.00000" << std::endl;
        wannierFileStream << std::endl;

        wannierFileStream << "KPLANE_BULK" << std::endl;
        wannierFileStream << "0.00  0.00  0.00   ! center of the k plane" << std::endl;
        wannierFileStream << "2.00  0.00  0.00   ! The first vector to define 3d k space plane" << std::endl;
        wannierFileStream << "0.00  2.00  0.00   ! The second vector to define 3d k space plane" << std::endl;
        wannierFileStream << std::endl;

        wannierFileStream << "KCUBE_BULK" << std::endl;
        wannierFileStream << "0.00  0.00  0.00   ! Original point for 3D k plane" << std::endl;
        wannierFileStream << "1.00  0.00  0.00   ! The first vector to define 3d k space plane" << std::endl;
        wannierFileStream << "0.00  1.00  0.00   ! The second vector to define 3d k space plane" << std::endl;
        wannierFileStream << "0.00  0.00  1.00   ! The third vector to define 3d k cube" << std::endl;


}














FHSystem::~FHSystem()
{

}






// Other functions





std::complex<double> U_os_nij(double U0, double alpha, double d)
{
    double u = 0;
    if (d==0.)
        u = U0;
    else if (alpha != -1)
    {      
        u = U0/pow(d,alpha);      
    }
    return std::complex<double>(u);
}









// DEPRECIATED FUNCTIONS

/*


FHSystem::FHSystem(std::string fileName, std::string folderName)
{

    enum troolean
    {
        yes,
        no,
        dontknow
    };

    enum WhatToRead
    {
        nothing,
        positions,
        symmetries,
        energyParameter_U_constant,
        energyParameter_U_bySite,
        energyParameter_U_bySiteType,
        energyParameter_U_GivenMatrix,
        energyParameter_tH_NNConstant,
        energyParameter_tH_GivenMatrix,
        energyParameter_J,
        diagonalIC_UpSpin,
        diagonalIC_DownSpin,
        numberOfParticles_UpSpin,
        numberOfParticles_DownSpin
    };

    systemName = fileName;
    systemFolder = folderName;
    std::string fileAdressName = fileName; 
    if (folderName !="")
    {
        fileAdressName = folderName + "/" + fileName;
    }

    UInteractionType = FHOppositeSpinInteractionType::Unone;
    VInteractionType = FHSameSpinInteractionType::Vnone;
    tunnelingType = FHTunnelingType::tHnone;
    tJExCType = tJExchangeCouplingType::tJJnone;
    externalPotentialType   = ExternalPotentialType::ExtPotNone;

    WhatToRead whatToR =  WhatToRead::nothing;  
    int nDim = -1;
    troolean orbitalHaveNames = troolean::dontknow;
    bool positionsReady = false;

    std::ifstream systFile (fileAdressName);
    if (systFile.is_open())
    {
        printf("== Reading %s \n", fileName.c_str());

        std::vector<unsigned int>           orbitalTypes    = std::vector<unsigned int>();
        std::vector<std::string>            orbitalNames    = std::vector<std::string>();
        std::vector<std::vector<double>>    posVecVec       = std::vector<std::vector<double>>();
        std::vector<std::vector<double>>    symVecVec       = std::vector<std::vector<double>>();
        std::vector<unsigned int>           initCond_up     = std::vector<unsigned int>();
        std::vector<unsigned int>           initCond_down   = std::vector<unsigned int>();
        std::vector<unsigned int>           particlesNumber(2, 0);

        std::vector<double> Jvector                         = std::vector<double>();
        std::vector<std::vector<double>> tHMatrix_vec           = std::vector<std::vector<double>>();
        std::vector<std::vector<double>> UMatrix_vec            = std::vector<std::vector<double>>();


        std::string line;
        while (std::getline(systFile, line))
        {
           if (line == "#POSITIONS")                                                //<- ready to read the positions
           {
                whatToR = WhatToRead::positions;
           }
           else if (line == "#POSITIONS_END")                                       //<- the positions stored in "pos" 
           {
                arma::dmat sitePos = armaMatfromVectorsVector(posVecVec);
                setPositions(sitePos);

                if (orbitalHaveNames == troolean::yes)
                {
                    std::vector<std::string>            orbitalNameTypes    = std::vector<std::string>();
                    uint numberOfNames = 0;
                    uint index = 0;
                    for (uint nameInt = 0; nameInt < orbitalNames.size(); nameInt++)
                    {
                        if (isIn( orbitalNames[nameInt], orbitalNameTypes, index))
                        {
                            orbitalTypes[nameInt] = index;
                        }
                        else
                        {
                            orbitalNameTypes.push_back(orbitalNames[nameInt]);
                            orbitalTypes[nameInt] = orbitalNameTypes.size() - 1;
                        }
                    }
                }

                initiateEmptyEnergyMatrices();
                positionsReady = true;
                whatToR = WhatToRead::nothing;
           }
           else if (line == "#SYMMETRIES")
           {
                whatToR = WhatToRead::symmetries;
           }
           else if (line == "#SYMMETRIES_END")
           {
                arma::dmat symmetryArmaMatrix = armaMatfromVectorsVector(symVecVec);
                setSymmetry(symmetryArmaMatrix);
                whatToR = WhatToRead::nothing;
           }
           else if (line == "#INITIALCONDITIONS_GIVENDIAGONAL")
           {
                initialCondition = FHInitialConditionType::givenDiagonal;
                whatToR = WhatToRead::diagonalIC_UpSpin;
           }
           else if (line == "#INITIALCONDITIONS_tbGroundState")
           {
                initialCondition = FHInitialConditionType::tbGroundState;
                whatToR = WhatToRead::numberOfParticles_UpSpin;
           }
           else if (line == "#ENERGYPARAMETER_U_CONSTANT")
           {
                whatToR = WhatToRead::energyParameter_U_constant;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteConstant;
           }
           else if (line == "#ENERGYPARAMETER_U_BYSITE")
           {
                whatToR = WhatToRead::energyParameter_U_bySite;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteGivenDiagonal;
           }
           else if (line == "#ENERGYPARAMETER_U_BYSITETYPE")
           {
                whatToR = WhatToRead::energyParameter_U_bySiteType;
                UInteractionType = FHOppositeSpinInteractionType::UonSiteBySiteType;
           }
           else if (line == "#ENERGYPARAMETER_tH_NEARESTNEIGHBORS")
           {
                whatToR = WhatToRead::energyParameter_tH_NNConstant;
                tunnelingType = FHTunnelingType::nearestNeighboursConstant; 
           }
           else if (line == "#ENERGYPARAMETER_tH_GIVENMATRIX")
           {
                whatToR = WhatToRead::energyParameter_tH_GivenMatrix;
                tunnelingType = FHTunnelingType::tunnelingGivenMatrix;
           }
            else if (line == "#ENERGYPARAMETER_J")
           {
                whatToR = WhatToRead::energyParameter_J;
           }
           else switch (whatToR)
           {
            case WhatToRead::positions :  // Read positions
            {

                if (orbitalHaveNames == troolean::dontknow)
                {
                    std::stringstream st(line);
                    std::string coordt;
                    orbitalHaveNames = troolean::no;
                    std::getline (st, coordt, ' ');
                    try {
                        std::stod(coordt);
                    }catch(const std::invalid_argument& e)  {
                        std::cout << "Orbitals have names !" << std::endl;
                        orbitalHaveNames = troolean::yes;
                    }
                }

                std::stringstream ss(line);
                std::string coord;
                std::vector<double> coordinates = std::vector<double>();

                if (orbitalHaveNames == troolean::yes)
                {
                    std::getline (ss, coord, ' ');
                    orbitalNames.push_back(coord);
                }

                while (std::getline (ss, coord, ' ')) 
                {
                    double d = std::stod(coord);
                    coordinates.push_back(d);
                }

                if (nDim == -1)
                {
                    nDim = coordinates.size();
                    posVecVec.push_back(coordinates);
                }
                else
                {
                    if (coordinates.size() == nDim)
                    {
                        posVecVec.push_back(coordinates);
                    }
                    else
                    {
                        printf("Dimension problem with line \"%s\"", line.c_str());
                        return;
                    }
                }


                break;
            }
            case WhatToRead::symmetries :
            {
                std::stringstream ss(line);
                std::string coord;
                std::vector<double> coordinates = std::vector<double>();
                while (std::getline (ss, coord, ' ')) 
                {
                    double d = std::stod(coord);
                    coordinates.push_back(d);
                }
                if (coordinates.size() == nDim)
                {
                    symVecVec.push_back(coordinates);
                }
                else
                {
                    printf("Dimension problem with line \"%s\"", line.c_str());
                    return;
                }
            }
            case WhatToRead::energyParameter_U_constant :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::getline (ss, U_string, ' ');
                double U0 = std::stod(U_string);
                
                setOppositeSpinInteractionType(FHOppositeSpinInteractionType::UonSiteConstant);
                setUInteraction_OnSite(U0);
                break;
            }
            case WhatToRead::energyParameter_U_bySite :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::vector<double> Uvector = std::vector<double>();
                while (std::getline (ss, U_string, ' ')) 
                {
                     Uvector.push_back(std::stod(U_string));
                }     
                arma::dvec Udiag = arma::dvec(Uvector);
                setOppositeSpinInteractionType(FHOppositeSpinInteractionType::UonSiteGivenDiagonal);
                setUInteraction_GivenDiagonal(Udiag);
                break;
            }
            case WhatToRead::energyParameter_U_bySiteType :
            {

                std::stringstream ss(line);
                std::string site_name;
                std::string U_value;
                std::getline (ss, site_name, ' ');
                std::getline (ss, U_value, ' ');
                std::vector<std::string>::iterator itName = std::find(orbitalNames.begin(), orbitalNames.end(), site_name);
                std::vector<double> Uvector = std::vector<double>(orbitalNames.size());
                Uvector[itName - orbitalNames.begin()] = std::stod(U_value);          
                break;
            }
            case WhatToRead::energyParameter_U_GivenMatrix :
            {
                std::stringstream ss(line);
                std::string U_string;
                std::vector<double> Uvector = std::vector<double>();
                while (std::getline (ss, U_string, ' ')) 
                {
                     Uvector.push_back(std::stod(U_string));
                }
                UMatrix_vec.push_back(Uvector);
                if (UMatrix.size() == posVecVec.size())
                {
                    arma::dmat UM = armaMatfromVectorsVector(UMatrix_vec);
                    setUInteraction_GivenMatrix(UM);
                }    
                break;
            }
            case WhatToRead::energyParameter_tH_NNConstant :
            {
                std::stringstream ss(line);
                std::string t_string;
                std::getline (ss, t_string, ' ');
                double tH1 = std::stod(t_string);
                setTunneling_BetweenNeighbours(tH1);
                break;
            }
            case WhatToRead::energyParameter_tH_GivenMatrix :
            {
                std::stringstream ss(line);
                std::string t_string;
                std::vector<double> tHvector = std::vector<double>();
                while (std::getline (ss, t_string, ' ')) 
                {
                     tHvector.push_back(std::stod(t_string));
                }
                tHMatrix_vec.push_back(tHvector);
                if (tHMatrix_vec.size() == posVecVec.size())
                {
                    arma::dmat tHM = armaMatfromVectorsVector(tHMatrix_vec);
                    setTunneling_GivenMatrix(tHM);
                }    
                break;
            }
            case WhatToRead::energyParameter_J :
            {
                break;
            }
           case WhatToRead::diagonalIC_UpSpin :  // Read  up spin particles
           {
                std::stringstream ss(line);
                std::string siteN;
                
                while (std::getline(ss, siteN, ' '))            //<- We read the line word by word
                {
                    //std::cout << "Reading initial positions of up-spin particles, siteN = " << siteN << std::endl;
                    initCond_up.push_back(std::stoi(siteN));    //<- And fill the inital condition vector with the site numbers
                    particlesNumber[0]++;
                }
                whatToR = WhatToRead::diagonalIC_DownSpin;      //<- At the end, the up-spin particles positions have been read.
                break;
            }
            case WhatToRead::diagonalIC_DownSpin :              //<- Else, if the down-spin particles positions have not been read yet
            {
                std::stringstream ss(line);
                std::string siteN;
                while (std::getline(ss, siteN, ' ')) 
                {
                    // std::cout << "Reading initial positions of down-spin particles, siteN = " << siteN << std::endl;
                    initCond_down.push_back(std::stoi(siteN)); 
                    particlesNumber[1]++;
                }
                whatToR = WhatToRead::nothing; 
                break;
            }
            case WhatToRead::numberOfParticles_UpSpin :
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[0] = std::stoi(nPart);
                whatToR = WhatToRead::numberOfParticles_DownSpin;
            }
            break;
            case WhatToRead::numberOfParticles_DownSpin :
            {
                std::stringstream ss(line);
                std::string nPart;

                std::getline(ss, nPart, ' ');
                particlesNumber[1] = std::stoi(nPart);
                whatToR = WhatToRead::nothing;
            }
            break;


           }
        }

        systFile.close();



    if (initialCondition == FHInitialConditionType::givenDiagonal)
    {
            std::vector<std::vector<unsigned int>> initCond = std::vector<std::vector<unsigned int>>(2);
            initCond[0] = initCond_up;
            initCond[1] = initCond_down;

            setInitialConditionsDiagonal(initCond);
    }
    else if (initialCondition == FHInitialConditionType::tbGroundState)
    {
        setInitialConditions(initialCondition, particlesNumber);
    }


    computeMatrices();

    }
    else
    {
        printf("Can't open %s", fileName.c_str());
    }


    originalJsonSystem = Json::Value();
}




/// @brief Write an ASCII version of the Fermi-Hubbard System
/// @param folderName emplacement of the folder
void                    FHSystem::writeFHSystemBis(std::string folderName)
{
    // Writing of everything on the file :
    std::string fileAdress =  systemName;

    std::cout << "In system, system name : "<< systemName << std::endl;
    if (folderName != "")
    {
        fileAdress = folderName + "/" + systemName;
    }
   std::ofstream streamSystem(fileAdress);

   std::cout << "In system, file adress : "<< fileAdress << std::endl;
   
    if (streamSystem.is_open())
    {
// I - POSITIONS
    // I - 1 POSITION OF SITES
        bool sitesHaveTags = sitesHaveNames();

        streamSystem << "#POSITIONS" << std::endl;
        for (uint itpos = 0; itpos < getnSites(); itpos++)
        {
            if (sitesHaveTags)
            {
                streamSystem << siteTypeNames[siteTypes[itpos]] << "\t" ;
            }
            for(uint itdim = 0 ; itdim < sitePositions.n_cols; itdim++)
            {
                streamSystem << sitePositions(itpos, itdim) << "\t";
            }
            streamSystem << std::endl;          
        }
        streamSystem << "#POSITIONS_END" << std::endl;

    // I - 2 SYMMETRIES
        if (isSymmetric())
        {
            streamSystem << "#SYMMETRIES" << std::endl;
            for (uint itSym = 0; itSym < symmetryVectors.n_rows; itSym++)
            {
                for(uint itdim = 0 ; itdim < symmetryVectors.n_cols; itdim++)
                {
                    streamSystem << sitePositions(itSym, itdim) << "\t";
                }
                streamSystem << std::endl;
            }    
            streamSystem << "#SYMMETRIES_END" << std::endl;       
        }
        
// II - ENERGY PARAMETERS
    // II-1 OPPOSITE SPIN INTERACTION
        switch(UInteractionType)
        {
            case FHOppositeSpinInteractionType::Unone :
            {
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteConstant:
            {
                streamSystem << "#ENERGYPARAMETER_U_CONSTANT" << std::endl;
                streamSystem << UMatrix(0,0).real() << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteGivenDiagonal:
            {
                streamSystem << "#ENERGYPARAMETER_U_BYSITE" << std::endl;
                for(uint iSite = 0; iSite<getnSites(); iSite++)
                {
                    streamSystem << UMatrix(iSite, iSite).real() << "\t";
                }
                streamSystem << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UonSiteBySiteType:
            {
                streamSystem << "#ENERGYPARAMETER_U_BYSITETYPE" << std::endl;
                for(uint iSiteType = 0; iSiteType < siteTypeNames.size(); iSiteType++)
                {
                    arma::uvec indSite = find(siteTypes = iSiteType, 1);
                    if (indSite.size() !=0)
                    {
                       streamSystem << siteTypeNames[iSiteType]<< "\t" << UMatrix(indSite(0)).real() << std::endl;
                    }            
                }
                streamSystem << std::endl;
                break;
            }
            case FHOppositeSpinInteractionType::UGivenMatrix:
            {
                streamSystem << "#ENERGYPARAMETER_U_GIVENMATRIX" << std::endl;
                for(uint itU1 = 0; itU1 < UMatrix.n_rows; itU1++)
                {
                    for(uint itU2 = 0; itU2 < UMatrix.n_cols; itU2++)
                    {       
                        streamSystem << std::real(UMatrix(itU1, itU2)) << "\t";
                    }
                    streamSystem << std::endl;
                }
                break;
            }
            default :
            {
                break;
            }
        }



        // II-2 TUNNELING
        switch(tunnelingType)
        {
            case FHTunnelingType::tHnone:
            {
                break;
            }
            case FHTunnelingType::nearestNeighboursConstant:
            {
                streamSystem << "#ENERGYPARAMETER_tH_NEARESTNEIGHBORS" << std::endl;
                arma::uvec q1 = find(tMatrix,1);
                //std::cout << q1 << std::endl;
                if (q1.n_rows !=0)
                {
                    double th = std::real( tMatrix(q1(0)) );
                    streamSystem << th << std::endl;
                }
                else
                {
                     streamSystem << 1 << std::endl;
                }
                break;
            }
            case FHTunnelingType::tunnelingGivenMatrix:
            {
                streamSystem << "#ENERGYPARAMETER_tH_GIVENMATRIX" << std::endl;
                for(uint itt1 = 0; itt1 < tMatrix.n_rows; itt1++)
                {
                    for(uint itt2 = 0; itt2 < tMatrix.n_cols; itt2++)
                    {       
                        streamSystem << std::real(tMatrix(itt1, itt2)) << "\t";
                    }
                    streamSystem << std::endl;
                }
                break;
            }
            default :
            {
                break;
            }
        }



// III - INITIAL CONDITIONS
        switch (initialCondition)
        {
            case FHInitialConditionType::tbGroundState:
                streamSystem << "#INITIALCONDITIONS_tbGroundState" << std::endl;
                break;
            case FHInitialConditionType::givenDiagonal:
            {
                streamSystem << "#INITIALCONDITIONS_GIVENDIAGONAL" << std::endl;
                for (uint itposUp = 0; itposUp < initConditions_upSpin.size(); itposUp++)
                {
                    if (initConditions_upSpin[itposUp])
                    {
                        streamSystem << itposUp << "\t";
                    }
                }
                streamSystem << std::endl;
                for (uint itposDown = 0; itposDown < initConditions_downSpin.size(); itposDown++)
                {
                    if (initConditions_downSpin[itposDown])
                    {
                        streamSystem << itposDown << "\t";
                    }
                }
                streamSystem << std::endl;       
                break;
            }   

        }


        streamSystem.close();
    }  
    else
    {
        std::cout << "Stream not open." << std::endl;
    }
}


*/

