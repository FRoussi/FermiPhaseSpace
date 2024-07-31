#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <omp.h>

#include <mpi.h>


#include "FHSystem.h"
#include "FHfTWASolver.h"



class FHfTWAMPISolver : public FHfTWASolver
{
public:



    FHfTWAMPISolver(FHSystem* fhs);


    void                        runSimulation();

    void                        set_run_MPI(bool yn, int argc, char** argv);

    ~FHfTWAMPISolver();

private:


// MPI parameters
    bool                        runMPI;
    int                         argc_copy;
    char**                      argv_copy;


};