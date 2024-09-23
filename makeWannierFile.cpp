// EXAMPLE : ./makeWannierFile <jsonSystemFile> 
//           ./makeWannierFile <jsonSystemFile> dat

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>


#include "include/FHSystem.h"

#include "include/json/json.h"


using namespace std;











// MAIN

int main(int argc, char* argv[])
{


   std::cout << "--- Writing entry files for Wannier Tools --- " << std::endl;


   // argc[1] : 
   // jsonFilePath = workFolder/jsonFile = workFolder/fileName.json;
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
        std::string jsonFile = jsonFilePath.substr(botDirPos + 1, jsonFilePath.length()); // get file
        std::string fileName = jsonFile.substr(0, jsonFile.length() - 5);

        std::cout << "== Loading data : " << std::endl;
        std::cout << "   -- Work file : " << jsonFile << std::endl;
        std::cout << "   -- Work folder : " << workFolder << std::endl;
        std::cout << "   -- File name : " << fileName << std::endl;



        std::string instructionsFileName = jsonFile;
        std::stringstream instructionsStream("");
        std::string systemFileName;
        FHSystem fhSystem;
        bool systemLoaded;
        Json::Value rootSim;


        std::cout << "   -- Reading work file "<< instructionsFileName << std::endl;  

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

        std::string wannierDatFileName = fileName + ".dat";
        cout << "Write Wannier file " << wannierDatFileName <<  endl; 
        fhSystem.writeWannierToolsDatFile(wannierDatFileName, "");
        cout << "afini" << endl; 
   }
   else
   {
         std::cout << "- Need a json file. "  << std::endl;
   }

    return 0;
}