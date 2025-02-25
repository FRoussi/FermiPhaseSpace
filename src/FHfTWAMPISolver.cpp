#include "../include/FHfTWAMPISolver.h"


FHfTWAMPISolver::FHfTWAMPISolver(FHSystem* fhS) : FHfTWASolver(fhS)
{

    runMPI = false;

}




/// @brief Launch simulation and write the results
void FHfTWAMPISolver::runSimulation( )
{



   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes

   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Get the rank of the process

   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len); // Get the name of the processor


    if (world_rank == 0)
    {
        std::cout << "   -- Start simulation." << std::endl;
    }




    uint nSites = fhSystem->getnSites();
    uint nEnergies = 3;     
    
    std::vector<double> timeLine = std::vector<double>(nTimePoints,0);
    //std::cout << "nSites=" << nSites << " , nTPoints="<< nTimePoints <<  std::endl;
    arma::cx_dcube  rhou_t,     rhod_t;
    arma::cx_dmat   rhoiiu,     rhoiid;
    arma::cx_dmat   rho2iiu,    rho2iid;
    arma::cx_dmat   rhoiurhoid;
    arma::cx_dmat   Enrgs;
    arma::cx_dmat   mkku, mkkd, mukku, mukkd;

    if (saveWholeDMatrix)
    {
        rhou_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t      = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu    = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid    = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
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
    if (saveTBStates)
    {
        mkku     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mkkd     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveMFStates)
    {
        mukku     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mukkd     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        Enrgs    = arma::cx_dmat(nEnergies,nTimePoints, arma::fill::zeros);
    }


    arma::cx_dmat halfDiagMatrix = arma::cx_dmat(0.5 * arma::dmat(nSites, nSites, arma::fill::eye), arma::dmat(nSites, nSites, arma::fill::zeros)); //<- Matrix of diagnoal 1/2


    uint trajPerNode = nTraj/world_size;
    uint modTrajPerNode = nTraj % world_size;

    std::vector<uint> trajNumberSplits = {0};
    for (uint r = 0; r < (world_size - modTrajPerNode); r++)
    {
        trajNumberSplits.push_back(trajNumberSplits[r] + trajPerNode);
    }
    for (uint r = world_size - modTrajPerNode; r < world_size  ; r++)
    {        
        trajNumberSplits.push_back(trajNumberSplits[r] + trajPerNode + 1) ;       
    }


    std::cout << "   -- MPI node " << world_rank << " will compute trajectories from " << trajNumberSplits[world_rank] << " to " << trajNumberSplits[world_rank+1] -1 << std::endl;
// SCATTER important variables :

    const auto startThreadSimTime{std::chrono::steady_clock::now()};

// PARALLEL LOOP ON TRAJECTORIES
#pragma omp parallel shared(rhou0_mean, rhod0_mean)
    {
        // Local (in threads) versions of important variables :
        arma::cx_dcube rhou_t_loc, rhod_t_loc;
        arma::cx_dmat rhoiiu_loc, rhoiid_loc; 
        arma::cx_dmat rho2iiu_loc, rho2iid_loc; 
        arma::cx_dmat rhoiurhoid_loc;
    
        arma::cx_dmat En_loc;   
        arma::cx_dmat mkku_loc, mkkd_loc, mukku_loc, mukkd_loc;
    


    if (saveWholeDMatrix)
    {
        rhou_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
        rhod_t_loc = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiagDMatrix)
    {
        rhoiiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rhoiid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveDiag2DMatrix)
    {
        rho2iiu_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        rho2iid_loc = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveNiuNid)
    {
        rhoiurhoid_loc  = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveTBStates)
    {
        mkku_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mkkd_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveMFStates)
    {
        mukku_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
        mukkd_loc     = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
    }
    if (saveEnergies)
    {
        En_loc          = arma::cx_dmat(nEnergies,nTimePoints, arma::fill::zeros);
    }

    arma::cx_dmat Nu_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat Nd_loc    = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhou      = arma::cx_dmat(nSites, nSites);
    arma::cx_dmat rhod      = arma::cx_dmat(nSites, nSites);


    #pragma omp for 
        for (unsigned int i_traj=trajNumberSplits[world_rank]; i_traj <trajNumberSplits[world_rank+1]; i_traj ++)
        {

            double currentSimulationTime = 0;


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


    /*
            for (unsigned int i = 0; i<nSites; i++)
            {
                for (unsigned int j = i+1; j<nSites; j++)
                {               
                    wnu(i,j) += wnM(i,j) * sqrtNoiseVarianceMatrixU(i,j);
                    wnu(j,i) += std::conj(wnM(i,j)) * sqrtNoiseVarianceMatrixU(j,i);
                    wnd(i,j) += wnM(j,i) * sqrtNoiseVarianceMatrixD(i,j);
                    wnd(j,i) += std::conj(wnM(j,i)) * sqrtNoiseVarianceMatrixD(j,i);
                }
            }*/



            if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::tbGroundState)
            {
                // Initial matrices in momentum space :
                arma::cx_dmat  mU = arma::diagmat(diagMu0) - halfDiag + wnu;
                arma::cx_dmat  mD = arma::diagmat(diagMd0) - halfDiag + wnd;

                // Initial matrices in position space :
                rhou = (*fhSystem->get_tbgs_V_Pointer()) * mU * (*fhSystem->get_tbgs_V_Pointer()).t() ;
                rhod = (*fhSystem->get_tbgs_V_Pointer()) * mD * (*fhSystem->get_tbgs_V_Pointer()).t();

            }
            else if (fhSystem->getInitialCondition() == FHSystem::FHInitialConditionType::mfGroundState)
            {
                // Initial matrices in momentum space :
                arma::cx_dmat  mU = arma::diagmat(diagMu0) - halfDiag + wnu;
                arma::cx_dmat  mD = arma::diagmat(diagMd0) - halfDiag + wnd;

                // Initial matrices in position space :
                rhou = (*fhSystem->get_mfgs_Vu_Pointer()) * mU * ((*fhSystem->get_mfgs_Vu_Pointer()).t()) ;
                rhod = (*fhSystem->get_mfgs_Vd_Pointer()) * mD * ((*fhSystem->get_mfgs_Vd_Pointer()).t()) ;

            }
            else 
            {
               // Initial rho matrices :
               rhou = rhou0_mean + wnu;  
               rhod = rhod0_mean + wnd;

            }


            // We compute the trajectory dynamic and add up the diagonal to the result
            for (unsigned int t = 0;  t<nTimePoints; t++)
            {   
                // We first store the density in the thread-local variables  
                timeLine[t] = t*dt;                                     //<- time line, stored too many times. To Fix ?
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
                    rho2iiu_loc.col(t) += arma::conj(rhou.diag()) % rhou.diag();    //<- Diagonal of the squared norm of up-spin density matrix 
                    rho2iid_loc.col(t) += arma::conj(rhod.diag()) % rhod.diag();    //<- Diagonal of the squared norm of down-spin density matrix
                }
                if (saveNiuNid)
                {
                    rhoiurhoid_loc.col(t)   += rhou.diag() % rhod.diag();
                }
                if (saveTBStates)
                {
                    mkku_loc.col(t)     += arma::diagvec((*fhSystem->get_tbgs_V_Pointer()).t() * rhou *(*fhSystem->get_tbgs_V_Pointer()));
                    mkkd_loc.col(t)     += arma::diagvec((*fhSystem->get_tbgs_V_Pointer()).t() * rhod *(*fhSystem->get_tbgs_V_Pointer()));
                }
                if (saveMFStates)
                {
                    mukku_loc.col(t)     += arma::diagvec((*fhSystem->get_mfgs_Vu_Pointer()).t() * rhou *(*fhSystem->get_mfgs_Vu_Pointer()));
                    mukkd_loc.col(t)     += arma::diagvec((*fhSystem->get_mfgs_Vd_Pointer()).t() * rhod *(*fhSystem->get_mfgs_Vd_Pointer()));
                }
                if (saveEnergies)
                {
                    Nu_loc =  rhou + halfDiagMatrix;
                    Nd_loc =  rhod + halfDiagMatrix;
                    En_loc.at(0,t)=  En_loc.at(0,t) + (arma::dot(Nu_loc.diag(), *U * Nd_loc.diag()) + arma::dot(Nd_loc.diag(), *U * Nu_loc.diag())) / arma::cx_double(2,0);
                    En_loc.at(1,t) = En_loc.at(1,t) + arma::accu((Nu_loc + Nd_loc) % *tM);
                    if (fhSystem->getExternalPotentialType() != FHExternalPotentialType::EPnull)
                    {
                        En_loc.at(2,t) = En_loc.at(2,t) + arma::accu((Nu_loc + Nd_loc) % *extPot);
                    }
                }

                // We do 'nStepsPerTimePoint' simulation steps
                for (unsigned int t_s = 0; t_s < nStepsPerTimePoint; t_s ++)
                {
                    if (integrator==fTWAIntegrator::RungeKutta2) step_RK2(currentSimulationTime, rhou, rhod, dt/nStepsPerTimePoint);                   //<- RK2 integrator chosen
                    else if (integrator==fTWAIntegrator::RungeKutta4_38) step_RK4_38(currentSimulationTime, rhou, rhod, dt/nStepsPerTimePoint);        //<- RK4 integrator chosen
                    currentSimulationTime += dt/nStepsPerTimePoint;
                }
            }


            // Little idea of where we are in the computation :
            unsigned int tt = std::floor(nTraj/10);
            if ((i_traj%tt) == 0)
            {
                std::cout << "     - Trajectory " << i_traj << " computed by thread " << omp_get_thread_num() << " of node " << world_rank <<  std::endl;    
            }

        }

        // We add up all thread-local variables into the global ones :
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
                for (uint ie =0; ie < (nEnergies-1); ie ++)
                {
                    Enrgs.row(ie) +=En_loc.row(ie);
                }
            }
        }

    }


    const auto endThreadSimTime{std::chrono::steady_clock::now()};

    const std::chrono::duration<double> thread_elapsed_seconds{endThreadSimTime - startThreadSimTime};

    std::cout << "   -- Computation of node " << world_rank << " ended, computation time : " << thread_elapsed_seconds.count() << " s." << std::endl;
    

    if ((world_size>1) &&  (world_rank == 0))
    {
        std::cout << "   -- Gathering results in node 0 : " << std::endl;
    }


    // Making array buffers in node 0
    std::complex<double> *nnt_received_buffer;
    std::complex<double> *nt_received_buffer;
    std::complex<double> *et_received_buffer;
    std::complex<double> *t_received_buffer;
    if (world_rank == 0)
    {   
        if (saveWholeDMatrix)
        {
            nnt_received_buffer = new std::complex<double>[nSites * nSites * nTimePoints];
        }
        nt_received_buffer = new std::complex<double>[nSites * nTimePoints];
        et_received_buffer = new std::complex<double>[nEnergies * nTimePoints];
        t_received_buffer = new std::complex<double>[nTimePoints];
    }

    // For each secondary  node : 
    for (uint procIndex = 1; procIndex < world_size; procIndex++)
    {
        if (saveWholeDMatrix)
        {
            if (world_rank == procIndex)         //<- Sending rho_u
            {           
                //std::cout << "     - Node " << procIndex << " sending rho_u to node 0." << std::endl;
                MPI_Send(rhou_t.memptr(), nSites * nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)                //<- Receiving rho_u
            {   
                //std::cout << "     - Node 0 receiving rho_u from node " << procIndex << "." << std::endl;
                MPI_Recv(nnt_received_buffer, nSites * nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dcube rho_received_arma = arma::cx_dcube(nnt_received_buffer, nSites, nSites, nTimePoints); 
                rhou_t +=  rho_received_arma;
            }
        

            if (world_rank == procIndex)         //<- Sending rho_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho_d to node 0." << std::endl;
                MPI_Send(rhod_t.memptr(), nSites * nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)                //<- Receiving rho_d
            {   
                //std::cout << "     - Node 0 receiving rho_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nnt_received_buffer, nSites * nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dcube rho_received_arma = arma::cx_dcube(nnt_received_buffer, nSites, nSites, nTimePoints); 
                rhod_t +=  rho_received_arma;
            }

        }
    

    
        if (saveDiagDMatrix)
        {
            if (world_rank == procIndex)         //<- Sending rhoii_u
            {
                std::cout << "     - Node " << procIndex << " sending data to node 0." << std::endl;
                //std::cout << "     - Node " << procIndex << " sending rho2ii_u to node 0." << std::endl;
                MPI_Send(rhoiiu.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rhoii_u
            {   
                std::cout << "     - Node 0 receiving data from node " << procIndex << "." << std::endl;
                //std::cout << "     - Node 0 receiving rho2ii_u from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat rhoii_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                rhoiiu += rhoii_received_arma;
            }

            if (world_rank == procIndex)         //<- Sending rhoii_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_d to node 0." << std::endl;
                MPI_Send(rhoiid.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rhoii_d
            {   
                //std::cout << "     - Node 0 receiving rho2ii_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat rhoii_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                rhoiid += rhoii_received_arma;          
            }
        }

        if (saveDiag2DMatrix)
        {
     
            if (world_rank == procIndex)         //<- Sending rho2ii_u
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_u to node 0." << std::endl;
                MPI_Send(rho2iiu.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_u
            {   
                //std::cout << "     - Node 0 receiving rho2ii_u from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat rhoii_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                rho2iiu += rhoii_received_arma;
            }

            if (world_rank == procIndex)         //<- Sending rho2ii_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_d to node 0." << std::endl;
                MPI_Send(rho2iid.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_d
            {   
                //std::cout << "     - Node 0 receiving rho2ii_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat rhoii_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                rho2iid += rhoii_received_arma;          
            }
        }

        if (saveNiuNid)
        {
            if (world_rank == procIndex)         //<- Sending rhoirhoid_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_d to node 0." << std::endl;
                MPI_Send(rhoiurhoid.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rhoii_d
            {   
                //std::cout << "     - Node 0 receiving rho2ii_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat rhoiurhoid_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                rhoiurhoid += rhoiurhoid_received_arma;          
            }
        }

        if (saveTBStates)
        {
     
            if (world_rank == procIndex)         //<- Sending rho2ii_u
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_u to node 0." << std::endl;
                MPI_Send(mkku.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_u
            {   
                //std::cout << "     - Node 0 receiving rho2ii_u from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat mkku_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                mkku += mkku_received_arma;
            }

            if (world_rank == procIndex)         //<- Sending rho2ii_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_d to node 0." << std::endl;
                MPI_Send(mkkd.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_d
            {   
                //std::cout << "     - Node 0 receiving rho2ii_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat mkkd_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                mkkd += mkkd_received_arma;          
            }
        }

        if (saveMFStates)
        {
     
            if (world_rank == procIndex)         //<- Sending rho2ii_u
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_u to node 0." << std::endl;
                MPI_Send(mukku.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_u
            {   
                //std::cout << "     - Node 0 receiving rho2ii_u from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat mukku_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                mukku += mukku_received_arma;
            }

            if (world_rank == procIndex)         //<- Sending rho2ii_d
            {
                //std::cout << "     - Node " << procIndex << " sending rho2ii_d to node 0." << std::endl;
                MPI_Send(mukkd.memptr(), nSites * nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving rho2ii_d
            {   
                //std::cout << "     - Node 0 receiving rho2ii_d from node " << procIndex << "." << std::endl;
                MPI_Recv(nt_received_buffer, nSites * nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat mukkd_received_arma = arma::cx_dmat(nt_received_buffer, nSites, nTimePoints); 
                mukkd += mukkd_received_arma;          
            }
        }

        if (saveEnergies)
        {
            if (world_rank == procIndex)         //<- Sending EU
            {
                //std::cout << "     - Node " << procIndex << " sending EU to node 0." << std::endl;
                MPI_Send(Enrgs.memptr(), nEnergies*nTimePoints, MPI_C_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
            if (world_rank == 0)        //<- Receiving EU
            {   
                //std::cout << "     - Node 0 receiving EU from node " << procIndex << "." << std::endl;
                MPI_Recv(et_received_buffer, nEnergies*nTimePoints, MPI_C_COMPLEX, procIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                arma::cx_dmat E_received_arma = arma::cx_dmat(et_received_buffer, nEnergies, nTimePoints); 
                Enrgs += E_received_arma;
            }
       
        }


    }

    if (world_rank == 0)
    {   
        if (saveWholeDMatrix) {delete [] nnt_received_buffer;}
        delete [] et_received_buffer;
        delete [] nt_received_buffer;
        delete [] t_received_buffer;
    }


    


    if (world_rank == 0)
    {

       // std::cout << "Computin results" << std::endl;

        // Computation of different observables :

        arma::cx_double nTraj_arma = arma::cx_double(nTraj,0);
        arma::cx_dcube  Nu, Nd   ;
        arma::cx_dmat   Niiu, Niid ;   
        arma::cx_dmat   N2iiu, N2iid;   
        arma::cx_dmat   Mkku ,  Mkkd ;
        arma::cx_dmat   Mukku,  Mukkd;
        arma::cx_dmat   NiuNid;
        arma::cx_dmat   Energies;


        if (saveWholeDMatrix)
        {
            Nu = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
            Nd = arma::cx_dcube(nSites, nSites, nTimePoints, arma::fill::zeros);
            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                Nu.slice(t) = rhou_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhou_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
                Nd.slice(t) = rhod_t.slice(t) / nTraj_arma + arma::cx_dmat(arma::size(rhod_t.slice(t)), arma::fill::eye) / arma::cx_double(2,0);
            }
        }
        if (saveDiagDMatrix)
        {
            Niiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            Niid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            arma::cx_dcolvec col1Half = arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(2,0);
            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                Niiu.col(t) = rhoiiu.col(t) / nTraj_arma + col1Half;
                Niid.col(t) = rhoiid.col(t) / nTraj_arma + col1Half;
            }
        }
        if (saveDiag2DMatrix)
        {
            N2iiu = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            N2iid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            arma::cx_dcolvec col1Quarter = arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                N2iiu.col(t) = rho2iiu.col(t) / nTraj_arma + (rhoiiu.col(t) + arma::conj(rhoiiu.col(t)))/ arma::cx_double(2*nTraj,0) + col1Quarter;
                N2iid.col(t) = rho2iid.col(t) / nTraj_arma + (rhoiid.col(t) + arma::conj(rhoiid.col(t)))/ arma::cx_double(2*nTraj,0) + col1Quarter;
            }
        }
        if (saveNiuNid)
        {
            NiuNid = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            arma::cx_dcolvec vec_4 = arma::cx_dcolvec(nSites, arma::fill::ones) / arma::cx_double(4,0);
            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                NiuNid.col(t) = rhoiurhoid.col(t) / nTraj_arma + (rhoiiu.col(t) + rhoiid.col(t))/ arma::cx_double(2*nTraj,0) + vec_4;
            }
        }
        if (saveTBStates)
        {

            Mkku = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            Mkkd = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                Mkku.col(t) = mkku.col(t) / nTraj_arma + halfVec;
                Mkkd.col(t) = mkkd.col(t) / nTraj_arma + halfVec;
            } 
        }
        if (saveMFStates)
        {

            Mukku = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);
            Mukkd = arma::cx_dmat(nSites, nTimePoints, arma::fill::zeros);

            for (unsigned int t = 0; t < nTimePoints; t++)
            {
                Mukku.col(t) = mukku.col(t) / nTraj_arma + halfVec;
                Mukkd.col(t) = mukkd.col(t) / nTraj_arma + halfVec;
            } 
        }
        if (saveEnergies)
        {
            Energies = arma::cx_dmat(nEnergies + 1, nTimePoints, arma::fill::zeros);
            for (uint ie = 0 ; ie < nEnergies; ie++)
            {
                Energies.row(ie) = Enrgs.row(ie) / nTraj;
            }
            Energies.row(nEnergies-1) = arma::sum(Energies, 0);
        }


       writeResults(Nu, Nd, Niiu, Niid, Mkku, Mkkd, Mukku, Mukkd, N2iiu, N2iid, NiuNid, Energies, timeLine);
    }


}

// MPI functions
void                        FHfTWAMPISolver::set_run_MPI(bool yn, int argc, char** argv)
{

    runMPI = yn;
    argc_copy = argc;
    argv_copy = argv;

}


FHfTWAMPISolver::~FHfTWAMPISolver()
{
}
