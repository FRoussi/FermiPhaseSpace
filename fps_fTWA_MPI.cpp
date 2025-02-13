//#include <stdio>
#include <iostream>
#include <string>
//#include <H5Cpp.h>
#include <omp.h>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>

#include <mpi.h>

//#include <filesystem>

#include "include/FHSystem.h"
#include "include/FHfTWAMPISolver.h"
#include "include/json/json.h"

#define DEF_MAX_THREAD_NUM 4

int main(int argc, char* argv[])
{


    MPI_Init(&argc, &argv);


   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes

   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Get the rank of the process

   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len); // Get the name of the processor


   if (world_rank == 0) //<- Main node will write nice things
   {

 

      std::cout << "--- Fermi-Hubbard Dynamic with fTWA and MPI --- " << std::endl;

      std::cout << "== FPS version:  " << FPS_VERSION_MAJOR << "." << FPS_VERSION_MINOR << "." << FPS_VERSION_PATCH << std::endl;

      arma::arma_version ver;
      std::cout << "== ARMA version: "<< ver.as_string() << std::endl;

   // Environment variables :
      char* env_OMP = std::getenv("OMP_NUM_THREADS");
      char* env_OPENBLAS_OMP = std::getenv("OPENBLAS_NUM_THREADS");

      std::cout << "== Environment variables : " << std::endl;
      if (env_OMP==NULL)
      {
         omp_set_num_threads(DEF_MAX_THREAD_NUM);
      }
      else
      {
         std::cout << "   -- OMP_NUM_THREADS      :        " << env_OMP << std::endl;
      }

      if (env_OPENBLAS_OMP==NULL)
      {
         std::cout << "   -- OPENBLAS_NUM_THREADS not defined, this could lead to slower parallelisation." << std::endl;
      }
      else
      {
         std::cout << "   -- OPENBLAS_NUM_THREADS :        " << env_OPENBLAS_OMP << std::endl;
      }

     //std::cout << "== OpenMPI active , version " << std::to_string(OMPI_MAJOR_VERSION) << "." << std::to_string(OMPI_MINOR_VERSION) << "." << std::to_string(OMPI_RELEASE_VERSION) << std::endl; 
      std::cout << "== MPI active"  << std::endl;
      std::cout << "   -- Number of nodes :        " << world_size << std::endl;

  std::cout << "== OpenMP active, version "  << std::to_string(_OPENMP) << std::endl; 
      //std::cout << "== OpenMP active " << std::endl; 
      std::cout << "   -- Number of processors :        " << omp_get_num_procs() << std::endl;
      std::cout << "   -- Number of threads available : " << omp_get_max_threads() << std::endl;


      std::cout << "== Double precision :   2^ -" << std::numeric_limits<double>::digits10 << std::endl;

   }

   //std::string work = "2Sites.json";
   if (argc > 1)
   {

      std::string jsonFilePath = std::string(argv[1]);
      std::size_t botDirPos = jsonFilePath.find_last_of("/");
      std::size_t fileNameSize = jsonFilePath.size() ;
      std::string workFolder = "";
      if (botDirPos < fileNameSize)
      {
         workFolder = jsonFilePath.substr(0, botDirPos+1); // get directory
      }
      std::string work = jsonFilePath.substr(botDirPos + 1, jsonFilePath.length()); // get file

      if (world_rank == 0)
      {
         std::cout << "== Loading data : " << std::endl;
         std::cout << "   -- Work file : " << work << std::endl;
         std::cout << "   -- Work folder : " << workFolder << std::endl;

      }


      std::string instructionsFileName = work;
      std::stringstream instructionsStream("");
      std::string systemFileName;
      FHSystem fhSystem;
      bool systemLoaded;
      Json::Value rootSim;

      if (world_rank == 0)
      {
         std::cout << "   -- Reading work file " << instructionsFileName << std::endl;  
      }

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
         std::cout << instructionsFileName << " not found." <<  " node " << world_rank << std::endl;
      }



      //std::cout << rootSim << std::endl;

      std::string simName = rootSim.get("sim_name", "defaultSimName").asString();
      std::string stringMethod = rootSim.get("method", "HF").asString();


      Json::Value systemJSON = rootSim.get("system", "");
      if (world_rank==0)
      {
         fhSystem = FHSystem(systemJSON, workFolder);
      }
      else
      {
         fhSystem = FHSystem(systemJSON, workFolder, false);
      }

      systemLoaded = true;
      if (world_rank==0)
      {
         std::cout << "   -- System loaded" << std::endl;
      }
  

      const auto startSimTime{std::chrono::steady_clock::now()};


      if (world_rank==0)
      {
            std::cout << "== Computation of fTWA dynamic "<< std::endl;
      }

      FHfTWAMPISolver fTWASolver = FHfTWAMPISolver(&fhSystem);
      fTWASolver.setResultFile(simName + "_fTWA_res", workFolder);
      fTWASolver.setJsonFile(rootSim);
      fTWASolver.setResultsWriter("json");
      fTWASolver.setResultsToSave(rootSim["results"]);
      fTWASolver.setMethodParameters(rootSim["parameters"]);
      fTWASolver.initiateSolver();
      fTWASolver.set_run_MPI(true, argc, argv);
      fTWASolver.runSimulation();



      const auto endSimTime{std::chrono::steady_clock::now()};

      const std::chrono::duration<double> elapsed_seconds{endSimTime - startSimTime};

      if (world_rank ==0)
      {
         std::cout << "== Computation time : " << elapsed_seconds.count() << " s." << std::endl;
      }



   }
   else
   {
      if (world_rank == 0)
      {
         std::cout << "Need a json file. "  << std::endl;
      }
   }

   if (world_rank == 0)
   {
      std::cout <<"== Program closing." << std::endl;
   }

   MPI_Finalize();
   
   return 0;
}
