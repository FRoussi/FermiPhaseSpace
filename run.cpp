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

   std::cout << "OpenMP active, "  << THREAD_NUM << " threads used on " << omp_get_max_threads() << " avalaible"<< std::endl;

   omp_set_num_threads(THREAD_NUM);


   std::string workFolder = "systemsFolder";
   std::string work = "2Sites";
   std::string instructionsFileName = work + ".txt";
   std::stringstream instructionsStream("");
   std::string systemFileName;
   FHSystem fhSystem;
   std::vector<double> nrjParameters;


   std::ifstream instructionsFile = std::ifstream(workFolder + "/" + instructionsFileName);
   std::string line, word;
   std::stringstream lineStream;
   while(std::getline(instructionsFile, line))
   {
      instructionsStream << line << "\n";
      std::cout << line.c_str() << "\n";
   }
   instructionsFile.close();

   int taskNumber = -1;
   bool systemLoaded = false;
   bool energyParametersLoaded = false;
   while(std::getline(instructionsStream, line))
   {
      if (line == "#SYSTEM_FILE")
         taskNumber = 0;
      else if (line == "#ENERGY_PARAMETERS")
      {
         if (systemLoaded)
            taskNumber = 1;   
         else
            std::cout << "Not reading energy parameters because no system loaded." << std::endl;
      }          
      else if (line == "#ED")
      {
         if (systemLoaded && energyParametersLoaded)
            taskNumber = 2;
         else
            std::cout << "Not reading ED instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#HF")
      {
         if (systemLoaded && energyParametersLoaded)
            taskNumber = 3;
         else
            std::cout << "Not reading HF instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#fTWA")
      {
         if (systemLoaded && energyParametersLoaded)
            taskNumber = 4;
         else
            std::cout << "Not reading fTWA instructions because something not loaded." << std::endl;      
      } 
      else if (line == "#GPSR")
      {
         if (systemLoaded && energyParametersLoaded)
            taskNumber = 5;
         else
            std::cout << "Not reading GPSR instructions because something not loaded." << std::endl;      
      }   

      else if (taskNumber == 0)
      {
         systemFileName = line;
         fhSystem = FHSystem(systemFileName, workFolder);
         systemLoaded = true;
         std::cout << "System loaded" << std::endl;
      }
      else if (taskNumber == 1)
      {
         lineStream = std::stringstream(line);
         while(std::getline(lineStream, word, ' '))
            nrjParameters.push_back(std::stod(word));
         
         fhSystem.setInteraction_OnSite(nrjParameters[0]);
         fhSystem.setTunneling_BetweenNeighbours(nrjParameters[1]);
         energyParametersLoaded = true;
         std::cout << "Energy parameters loaded, U=" << nrjParameters[0] << ", J=" << nrjParameters[1] << std::endl;
      }
      else if (taskNumber == 2)
      {
         lineStream = std::stringstream(line);
         std::vector<double> edParameters;
         while(std::getline(lineStream, word, ' '))
            edParameters.push_back(std::stod(word));
         
         printf("-- ED dynamic \n");
         FHEDSolver edSolver = FHEDSolver(&fhSystem);
         edSolver.initiateSolver();
         edSolver.setTimeParameters(edParameters[0], edParameters[1], edParameters[2]);
         edSolver.runSimulation();
      }
      else if (taskNumber == 3)
      {
         lineStream = std::stringstream(line);
         std::vector<double> hfParameters;
         while(std::getline(lineStream, word, ' '))
            hfParameters.push_back(std::stod(word));
         
         printf("-- HF dynamic \n");
         FHHFSolver hfSolver = FHHFSolver(&fhSystem);
         hfSolver.initiateSolver();
         hfSolver.setTimeParameters(hfParameters[0], hfParameters[1], hfParameters[2]);
         hfSolver.runSimulation();
      }
      else if (taskNumber == 4)
      {
         lineStream = std::stringstream(line);
         std::vector<double> ftwaParameters;
         while(std::getline(lineStream, word, ' '))
            ftwaParameters.push_back(std::stod(word));
         
         printf("-- fTWA dynamic \n");
         FHfTWASolver ftwaSolver = FHfTWASolver(&fhSystem);
         ftwaSolver.initiateSolver();
         ftwaSolver.setTimeParameters(ftwaParameters[0], ftwaParameters[1], ftwaParameters[2]);
         ftwaSolver.setfTWAParameters(ftwaParameters[3]);
         ftwaSolver.runSimulation();
      }
      else if (taskNumber == 5)
      {
         lineStream = std::stringstream(line);
         std::vector<double> gpsrParameters;
         while(std::getline(lineStream, word, ' '))
            gpsrParameters.push_back(std::stod(word));
         
         printf("-- GPSR dynamic \n");
         FHGPSRSolver gpsrSolver = FHGPSRSolver(&fhSystem);
         gpsrSolver.initiateSolver();
         gpsrSolver.setTimeParameters(gpsrParameters[0], gpsrParameters[1], gpsrParameters[2]);
         gpsrSolver.setGPSRParameters(gpsrParameters[3]);
         gpsrSolver.runSimulation();
      }
   }

   std::cout <<"Final entry. Looks good." << " \n";

   return 0;
}
