//#include <stdio>
#include <iostream>
#include <string>
#include <armadillo>
#include <omp.h>

#include "FHSystem.h"
#include "FHEDSolver.h"
#include "FHHFSolver.h"
#include "FHfTWASolver.h"
#include "FHGPSRSolver.h"

#define THREAD_NUM 2

int main()
{

   std::cout << "--- Fermi-Hubbard Dynamic with Phase-space methods --- " << std::endl;
   arma::arma_version ver;
   std::cout << "ARMA version: "<< ver.as_string() << std::endl;

   std::cout << "OpenMP active, " << THREAD_NUM << " threads used on " << omp_get_max_threads() << " avalaible"<< std::endl;
   omp_set_num_threads(THREAD_NUM);





   FHSystem twoSites =  FHSystem("syst2Sites.txt", "systemsFolder");
   twoSites.setInteraction_OnSite(1.);
   twoSites.setTunneling_BetweenNeighbours(1.);
   unsigned int n = twoSites.getnSites();

   std::vector<unsigned int> pU = {0};
   std::vector<unsigned int> pD = {1};

   printf("-- ED dynamic \n");
   FHEDSolver edTry = FHEDSolver(&twoSites);
   edTry.initiateSolver();
   edTry.setTimeParameters(1., 0.01);
   edTry.runSimulation();


   printf("-- HF dynamic \n");
   FHHFSolver fhTry = FHHFSolver(&twoSites);  
   fhTry.initiateSolver();
   fhTry.setTimeParameters(1., 0.01, 5);
   fhTry.runSimulation();

   printf("-- fTWA dynamic \n");
   FHfTWASolver fTWATry = FHfTWASolver(&twoSites);  
   fTWATry.setTimeParameters(1., 0.01, 5);
   fTWATry.setfTWAParameters(10);
   fTWATry.initiateSolver(pU, pD);
   fTWATry.runSimulation();

   printf("-- GPSR dynamic \n");
   FHGPSRSolver GPSRTry = FHGPSRSolver(&twoSites);  
   GPSRTry.setTimeParameters(1., 0.01, 1);
   GPSRTry.setGPSRParameters(10);
   GPSRTry.initiateSolver(pU, pD);
   GPSRTry.runSimulation();
   

   std::cout <<"Final entry. Looks good." << " \n";


   return 0;
}
