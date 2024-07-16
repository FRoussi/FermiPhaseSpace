//#include <stdio>
#include <iostream>
#include <string>
#include <omp.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>

//#include <filesystem>

#include "include/FHSystem.h"
#include "include/FHEDSolver.h"
#include "include/FHHFSolver.h"
#include "include/FHfTWASolver.h"
#include "include/FHGPSRSolver.h"
#include "include/json/json.h"

#define DEF_MAX_THREAD_NUM 4

int main(int argc, char* argv[])
{

   std::cout << "--- Fermi-Hubbard Dynamic with Phase-space methods --- " << std::endl;


   arma::arma_version ver;
   std::cout << "== ARMA version: "<< ver.as_string() << std::endl;

// Threads of openmp :
   char* env_OMP = std::getenv("OMP_NUM_THREADS");

   if (env_OMP==NULL)
   {
      omp_set_num_threads(DEF_MAX_THREAD_NUM);
   }
   else
   {
      std::cout << "   -- OMP_NUM_THREADS :        " << env_OMP << std::endl;
   }

  std::cout << "== OpenMP active, version "  << std::to_string(_OPENMP) << std::endl; 
  std::cout << "   -- Number of processors :        " << omp_get_num_procs() << std::endl;
  std::cout << "   -- Number of threads available : " << omp_get_max_threads() << std::endl;


   char* env_OPENBLAS_OMP = std::getenv("OPENBLAS_NUM_THREADS");

   char* env_ARMA_OMP_1 = std::getenv("ARMA_USE_OPENMP");
   char* env_ARMA_OMP_2 = std::getenv("ARMA_DONT_USE_OPENMP");
   char* env_ARMA_OMP_3 = std::getenv("ARMA_OPENMP_THREADS");


  std::cout << "== arma-openmp variables : " << std::endl;
   if (env_ARMA_OMP_1!=NULL)
   {
        std::cout << "   -- ARMA_USE_OPENMP :        " << env_ARMA_OMP_1 << std::endl;
   }
      if (env_ARMA_OMP_2!=NULL)
   {
        std::cout << "   -- ARMA_DONT_USE_OPENMP :        " << env_ARMA_OMP_2 << std::endl;
   }
      if (env_ARMA_OMP_3!=NULL)
   {
        std::cout << "   -- ARMA_OPENMP_THREADS :        " << env_ARMA_OMP_3 << std::endl;
   }



// Threads used by openblas :
   //openblas_set_num_threads(1);


   std::cout << "== Double precision :   2^ -" << std::numeric_limits<double>::digits10 << std::endl;




   //std::string work = "2Sites.json";
   if (argc > 1)
   {

      std::string jsonFilePath = std::string(argv[1]);
      std::size_t botDirPos = jsonFilePath.find_last_of("/") ;
      std::string workFolder = jsonFilePath.substr(0, botDirPos); // get directory
      std::string work = jsonFilePath.substr(botDirPos + 1, jsonFilePath.length()); // get file

        std::cout << "== Work file : " << work << std::endl;
        std::cout << "== Work folder : " << workFolder << std::endl;



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


      std::string instructionsFileName = work;
      std::stringstream instructionsStream("");
      std::string systemFileName;
      FHSystem fhSystem;
      bool systemLoaded;
      Json::Value rootSim;


      std::cout << "- Reading work file "<< instructionsFileName << std::endl;  

      std::string file = workFolder + "/" + instructionsFileName;

      std::ifstream instructionsFile = std::ifstream(file);



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



      //std::cout << rootSim << std::endl;

      std::string simName = rootSim.get("sim_name", "defaultSimName").asString();
      std::string stringMethod = rootSim.get("method", "HF").asString();

      dynamicMethod chosenMethod = dynamicMethod::HF;
      if (stringMethod == "fTWA")         chosenMethod = dynamicMethod::fTWA; 
      else if (stringMethod == "GPSR")    chosenMethod = dynamicMethod::GPSR;
      else if (stringMethod == "ED")      chosenMethod = dynamicMethod::ED;


      Json::Value systemJSON = rootSim.get("system", "");
      fhSystem = FHSystem(systemJSON, workFolder);
      systemLoaded = true;
      std::cout << "- System loaded" << std::endl;

      const auto startSimTime{std::chrono::steady_clock::now()};

      switch (chosenMethod)
      {
      case dynamicMethod::HF:
      {
            double param_dt = rootSim["parameters"].get("dt", 0.01).asDouble();
            double param_tm = rootSim["parameters"].get("t_max", 2).asDouble();
            int param_sptp = rootSim["parameters"].get("stepsParTimePoints", 1).asInt();

            std::cout << "- Computation of HF dynamic " << std::endl;
            FHHFSolver hfSolver = FHHFSolver(&fhSystem);
            hfSolver.setResultFile(simName + "_HF_res", workFolder);
            hfSolver.setJsonFile(rootSim);
            hfSolver.writeResults("json");
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

            std::cout << "- Computation of fTWA dynamic "<< std::endl;
            FHfTWASolver fTWASolver = FHfTWASolver(&fhSystem);
            fTWASolver.setResultFile(simName + "_fTWA_res", workFolder);
            fTWASolver.setJsonFile(rootSim);
            fTWASolver.writeResults("json");
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

            std::cout <<  "- Computation of GPSR dynamic "<< std::endl;
            FHGPSRSolver GPSRSolver = FHGPSRSolver(&fhSystem);
            GPSRSolver.setResultFile(simName + "_GPSR_res", workFolder);
            GPSRSolver.setJsonFile(rootSim);
            GPSRSolver.writeResults("json");
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

            std::cout <<  "- Computation of ED dynamic "<< std::endl;
            FHEDSolver EDSolver = FHEDSolver(&fhSystem);
            EDSolver.setResultFile(simName + "_ED_res", workFolder);
            EDSolver.setJsonFile(rootSim);
            EDSolver.writeResults("json");
            EDSolver.initiateSolver();
            EDSolver.setTimeParameters(param_tm, param_dt, param_sptp);
            EDSolver.runSimulation();

            break;
      }
      default:
         std::cout << "Simulation method not known." << std::endl;
         break;
      }

      const auto endSimTime{std::chrono::steady_clock::now()};

      const std::chrono::duration<double> elapsed_seconds{endSimTime - startSimTime};



      std::cout << "- Computation time : " << elapsed_seconds.count() << " s." << std::endl;


   }
   else
   {
         std::cout << "- Need a json file. "  << std::endl;
   }


   std::cout <<"- Program closing." << std::endl;
   return 0;
}
