//#include <stdio>
#include <iostream>
#include <string>
#include <armadillo>
#include <omp.h>

#include "include/FHSystem.h"
#include "include/FHEDSolver.h"
#include "include/FHHFSolver.h"
#include "include/FHfTWASolver.h"
#include "include/FHGPSRSolver.h"
#include "include/json/json.h"

#define THREAD_NUM 4

int main(int argc, char* argv[])
{

   std::cout << "--- Fermi-Hubbard Dynamic with Phase-space methods --- " << std::endl;


   arma::arma_version ver;
   std::cout << "ARMA version: "<< ver.as_string() << std::endl;

   std::cout << "OpenMP active, "  << THREAD_NUM << " threads used on " << omp_get_max_threads() << " avalaible"<< std::endl;

   std::cout << "Double precision = " << std::numeric_limits<double>::digits10 << std::endl;

   omp_set_num_threads(THREAD_NUM);


   std::string work = "2Sites.json";
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

   enum dynamicMethod
   {
      HF,
      fTWA,
      GPSR,
      ED
   };


   std::string workFolder = "systemsFolder";
   std::string instructionsFileName = work;
   std::stringstream instructionsStream("");
   std::string systemFileName;
   FHSystem fhSystem;
   bool systemLoaded;
   Json::Value rootSim;


   std::cout << "-- Reading work file "<< instructionsFileName << " --" << std::endl;  
   
   std::ifstream instructionsFile = std::ifstream(workFolder + "/" + instructionsFileName);

   if (instructionsFile.good())
   {
      std::string line, word;
      std::stringstream lineStream;

      Json::CharReaderBuilder builder;
      builder["collectComments"] = true;
      JSONCPP_STRING errs;
      if(!Json::parseFromStream(builder, instructionsFile, &rootSim, &errs))
      {
         std::cout << errs << std::endl;
      }
   }
   else
   {
      std::cout << instructionsFileName << " not found." << std::endl;
   }



   std::cout << rootSim << std::endl;

   std::string simName = rootSim.get("sim_name", "defaultSimName").asString();
   std::string stringMethod = rootSim.get("method", "HF").asString();

   dynamicMethod chosenMethod = dynamicMethod::HF;
   if (stringMethod == "fTWA")         chosenMethod = dynamicMethod::fTWA; 
   else if (stringMethod == "GPSR")    chosenMethod = dynamicMethod::GPSR;
   else if (stringMethod == "ED")      chosenMethod = dynamicMethod::ED;


   Json::Value systemJSON = rootSim.get("system", "");
   fhSystem = FHSystem(systemJSON, workFolder);
   systemLoaded = true;
   std::cout << "System loaded" << std::endl;

   switch (chosenMethod)
   {
   case dynamicMethod::HF:
   {
         double param_dt = rootSim["parameters"].get("dt", 0.01).asDouble();
         double param_tm = rootSim["parameters"].get("t_max", 2).asDouble();
         int param_sptp = rootSim["parameters"].get("stepsParTimePoints", 1).asInt();
         
         printf("-- HF dynamic \n");
         FHHFSolver hfSolver = FHHFSolver(&fhSystem);
         hfSolver.setResultFile(simName + "_HF_res", workFolder);
         hfSolver.initiateSolver();
         hfSolver.setTimeParameters(param_tm, param_dt, param_sptp);
         hfSolver.runSimulation();

         break;
   }   
   case dynamicMethod::fTWA:
   {
         double param_dt = rootSim["parameters"].get("dt", 0.01).asDouble();
         double param_tm = rootSim["parameters"].get("t_max", 2).asDouble();
         int param_sptp = rootSim["parameters"].get("stepsParTimePoints", 1).asInt();
         int param_nTraj = rootSim["parameters"].get("n_traj", 10).asInt();
         
         printf("-- fTWA dynamic \n");
         FHfTWASolver fTWASolver = FHfTWASolver(&fhSystem);
         fTWASolver.setResultFile(simName + "_fTWA_res", workFolder);
         fTWASolver.initiateSolver();
         fTWASolver.setTimeParameters(param_tm, param_dt, param_sptp);
         fTWASolver.setfTWAParameters(param_nTraj);
         fTWASolver.runSimulation();

         break;
   }
   case dynamicMethod::GPSR:
   {
         double param_dt = rootSim["parameters"].get("dt", 0.01).asDouble();
         double param_tm = rootSim["parameters"].get("t_max", 2).asDouble();
         int param_sptp = rootSim["parameters"].get("stepsParTimePoints", 1).asInt();
         int param_nTraj = rootSim["parameters"].get("n_traj", 10).asInt();

         printf("-- GPSR dynamic \n");
         FHGPSRSolver GPSRSolver = FHGPSRSolver(&fhSystem);
         GPSRSolver.setResultFile(simName + "_GPSR_res", workFolder);
         GPSRSolver.initiateSolver();
         GPSRSolver.setTimeParameters(param_tm, param_dt, param_sptp);
         GPSRSolver.setGPSRParameters(param_nTraj);
         GPSRSolver.runSimulation();

         break;
   }
   case dynamicMethod::ED:
   {
         double param_dt = rootSim["parameters"].get("dt", 0.01).asDouble();
         double param_tm = rootSim["parameters"].get("t_max", 2).asDouble();
         int param_sptp = rootSim["parameters"].get("stepsParTimePoints", 1).asInt();
         
         printf("-- ED dynamic \n");
         FHEDSolver EDSolver = FHEDSolver(&fhSystem);
         EDSolver.setResultFile(simName + "_ED_res", workFolder);
         EDSolver.initiateSolver();
         EDSolver.setTimeParameters(param_tm, param_dt, param_sptp);
         EDSolver.runSimulation();

         break;
   }
   default:
      std::cout << "Simulation method not known." << std::endl;
      break;
   }




   std::cout <<"Final entry. Looks good." << " \n";



   return 0;
}
