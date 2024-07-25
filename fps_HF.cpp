//#include <stdio>
#include <iostream>
#include <string>
//#include <H5Cpp.h>
#include <omp.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>


//#include <filesystem>

#include "include/FHSystem.h"
#include "include/FHHFSolver.h"
#include "include/json/json.h"

#define DEF_MAX_THREAD_NUM 4

int main(int argc, char* argv[])
{



   std::cout << "--- Fermi-Hubbard Dynamic with Hartree-Fock method --- " << std::endl;


   arma::arma_version ver;
   std::cout << "== ARMA version: "<< ver.as_string() << std::endl;








// Environment variables :
   char* env_OMP = std::getenv("OMP_NUM_THREADS");
   char* env_OPENBLAS_OMP = std::getenv("OPENBLAS_NUM_THREADS");



  //std::cout << "== OpenMP active, version "  << std::to_string(_OPENMP) << std::endl; 
  std::cout << "== OpenMP active " << std::endl; 
  std::cout << "   -- Number of processors :        " << omp_get_num_procs() << std::endl;
  std::cout << "   -- Number of threads available : " << omp_get_max_threads() << std::endl;






// Threads used by openblas :
   //openblas_set_num_threads(1);


   std::cout << "== Double precision :   2^ -" << std::numeric_limits<double>::digits10 << std::endl;




   //std::string work = "2Sites.json";
   if (argc > 1)
   {

      std::string jsonFilePath = std::string(argv[1]);
      std::size_t botDirPos = jsonFilePath.find_last_of("/") ;
      std::size_t fileNameSize = jsonFilePath.size() ;
      std::string workFolder = "";
      if (botDirPos < fileNameSize)
      {
         workFolder = jsonFilePath.substr(0, botDirPos+1); // get directory
      }
      std::string work = jsonFilePath.substr(botDirPos + 1, jsonFilePath.length()); // get file

      std::cout << "== Work file : " << work << std::endl;
      std::cout << "== Work folder : " << workFolder << std::endl;




      std::string instructionsFileName = work;
      std::stringstream instructionsStream("");
      std::string systemFileName;
      FHSystem fhSystem;
      bool systemLoaded;
      Json::Value rootSim;


      std::cout << "- Reading work file "<< instructionsFileName << std::endl;  

      std::string file = workFolder + instructionsFileName;

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




      Json::Value systemJSON = rootSim.get("system", "");
      fhSystem = FHSystem(systemJSON, workFolder);
      systemLoaded = true;
      std::cout << "- System loaded" << std::endl;

      const auto startSimTime{std::chrono::steady_clock::now()};



      std::cout << "- Computation of HF dynamic " << std::endl;
      FHHFSolver hfSolver = FHHFSolver(&fhSystem);
      hfSolver.setResultFile(simName + "_HF_res", workFolder);
      hfSolver.setJsonFile(rootSim);
      hfSolver.setResultsWriter("json");
      hfSolver.initiateSolver();
      hfSolver.setMethodParameters(rootSim["parameters"]);
      hfSolver.runSimulation();

         

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
