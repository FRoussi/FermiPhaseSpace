#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "json/json.h"

class ResultWriterJson
{
public:

    ResultWriterJson(); 


    Json::Value writeLine(std::vector<double> vec);
    Json::Value writeLine(std::vector<unsigned int> vec);
    Json::Value writeLine(arma::drowvec vec);
    Json::Value writeLine(arma::cx_drowvec vec);
    Json::Value writeLineInReal(arma::cx_drowvec vec);


    ~ResultWriterJson();

   

private:


};