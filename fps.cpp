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

#define THREAD_NUM 4

int main(int argc, char* argv[])
{

   std::cout << "--- Fermi-Hubbard Dynamic with Phase-space methods --- " << std::endl;


   arma::arma_version ver;
   std::cout << "ARMA version: "<< ver.as_string() << std::endl;

   std::cout << "OpenMP active, "  << THREAD_NUM << " threads used on " << omp_get_max_threads() << " avalaible"<< std::endl;

   std::cout << "Double precision = " << std::numeric_limits<double>::digits10 << std::endl;

   omp_set_num_threads(THREAD_NUM);


   std::string work = "2Sites";
    if (argc > 1)
    {
        work = std::string(argv[1]);
        std::cout << "== Work file : " << work << std::endl;
    }

   enum task
   {
      nothing,
      readSystemFile,
      compute_ED,
      compute_FH,
      compute_fTWA,
      compute_GPSR
   };


   std::string workFolder = "systemsFolder";
   std::string instructionsFileName = work;
   std::stringstream instructionsStream("");
   std::string systemFileName;
   FHSystem fhSystem;
   std::vector<double> nrjParameters;

   std::cout << "-- Reading work file "<< instructionsFileName << " --" << std::endl;
   std::ifstream instructionsFile = std::ifstream(workFolder + "/" + instructionsFileName);
   std::string line, word;
   std::stringstream lineStream;

   while(std::getline(instructionsFile, line))
   {
      instructionsStream << line << "\n";
      std::cout << line.c_str() << "\n";
   }
   instructionsFile.close();
   std::cout << "---------------------------" << std::endl;

   task TODO = task::nothing;
   bool systemLoaded = false;
   while(std::getline(instructionsStream, line))
   {
      if (line == "#SYSTEM_FILE")
        TODO = task::readSystemFile;        
      else if (line == "#ED")
      {
         if (systemLoaded )
            TODO = task::compute_ED;
         else
            std::cout << "Not reading ED instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#HF")
      {
         if (systemLoaded )
            TODO = task::compute_FH;
         else
            std::cout << "Not reading HF instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#fTWA")
      {
         if (systemLoaded)
            TODO = task::compute_fTWA;
         else
            std::cout << "Not reading fTWA instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#GPSR")
      {
         if (systemLoaded )
            TODO = task::compute_GPSR;
         else
            std::cout << "Not reading GPSR instructions because something not loaded." << std::endl;      
      }   

      else if (TODO == task::readSystemFile)
      {
         systemFileName = line;
         fhSystem = FHSystem(systemFileName, workFolder);
         systemLoaded = true;
         std::cout << "System loaded" << std::endl;
         //std::cout << "Initial conditions spin up   : [" << fhSystem.getInitialCondition_upSpin()[0] << fhSystem.getInitialCondition_upSpin()[1] << "]"<< std::endl;
         //std::cout << "Initial conditions spin down : [" << fhSystem.getInitialCondition_downSpin()[0] << fhSystem.getInitialCondition_downSpin()[1] << "]"<< std::endl;
         
      }

      else if (TODO == task::compute_ED)
      {
         lineStream = std::stringstream(line);
         std::vector<double> edParameters;
         while(std::getline(lineStream, word, ' '))
            edParameters.push_back(std::stod(word));
         
         printf("-- ED dynamic \n");
         FHEDSolver edSolver = FHEDSolver(&fhSystem);
         edSolver.setResultFile(work + "_ED", workFolder);
         edSolver.initiateSolver();
         edSolver.setTimeParameters(edParameters[0], edParameters[1], edParameters[2]);
         edSolver.runSimulation();

      }
      else if (TODO == task::compute_FH)
      {
         lineStream = std::stringstream(line);
         std::vector<double> hfParameters;
         while(std::getline(lineStream, word, ' '))
            hfParameters.push_back(std::stod(word));
         
         printf("-- HF dynamic \n");
         FHHFSolver hfSolver = FHHFSolver(&fhSystem);
         hfSolver.setResultFile(work + "_HF", workFolder);
         hfSolver.initiateSolver();
         hfSolver.setTimeParameters(hfParameters[0], hfParameters[1], hfParameters[2]);
         hfSolver.runSimulation();
      }
      else if (TODO == task::compute_fTWA)
      {
         lineStream = std::stringstream(line);
         std::vector<double> ftwaParameters;
         while(std::getline(lineStream, word, ' '))
            ftwaParameters.push_back(std::stod(word));
         
         printf("-- fTWA dynamic \n");
         FHfTWASolver ftwaSolver = FHfTWASolver(&fhSystem);
         ftwaSolver.setResultFile(work + "_fTWA", workFolder);
         ftwaSolver.initiateSolver();
         ftwaSolver.setTimeParameters(ftwaParameters[0], ftwaParameters[1], ftwaParameters[2]);
         ftwaSolver.setfTWAParameters(ftwaParameters[3]);
         ftwaSolver.runSimulation();
      }
      else if (TODO == task::compute_GPSR)
      {
         lineStream = std::stringstream(line);
         std::vector<double> gpsrParameters;
         while(std::getline(lineStream, word, ' '))
            gpsrParameters.push_back(std::stod(word));
         
         printf("-- GPSR dynamic \n");
         FHGPSRSolver gpsrSolver = FHGPSRSolver(&fhSystem);
         gpsrSolver.setResultFile(work + "_GPSR", workFolder);
         gpsrSolver.initiateSolver();
         gpsrSolver.setTimeParameters(gpsrParameters[0], gpsrParameters[1], gpsrParameters[2]);
         gpsrSolver.setGPSRParameters(gpsrParameters[3]);
         gpsrSolver.runSimulation();
      }
   }

   std::cout <<"Final entry. Looks good." << " \n";

   return 0;
}
