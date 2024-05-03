#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

class ResultWriter
{
public:

    ResultWriter(const std::string & fileName, const std::string & folderName = ""); 

    void newResults();

    void writeLine(std::string lineName, std::vector<double> vec);
    void writeLine(std::string lineName, std::vector<unsigned int> vec);
    void writeLine(std::string lineName, arma::drowvec vec);
    void writeLine(std::string lineName, arma::cx_drowvec vec);
    void writeLineInReal(std::string lineName, arma::cx_drowvec vec);


    ~ResultWriter();

   

private:

    std::string resultFolder;
    std::string resultFileName; 
    std::string fileAdressName;

};