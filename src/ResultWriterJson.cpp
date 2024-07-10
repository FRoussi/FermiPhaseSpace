#include "../include/ResultWriterJson.h"



ResultWriterJson::ResultWriterJson()
{



}




Json::Value ResultWriterJson::writeLine(std::vector<double> vec)
{
    Json::Value jsonVec(Json::arrayValue);
    for (uint i = 0; i < vec.size(); i++)
    {
        jsonVec.append(vec[i]);
    }  

    return jsonVec;
}

Json::Value ResultWriterJson::writeLine(std::vector<unsigned int> vec)
{
    Json::Value jsonVec(Json::arrayValue);
    for (uint i = 0; i < vec.size(); i++)
    {
        jsonVec.append(vec[i]);
    }  

    return jsonVec;
}

Json::Value ResultWriterJson::writeLine(arma::drowvec vec)
{
    Json::Value jsonVec(Json::arrayValue);
    for (uint i = 0; i < vec.size(); i++)
    {
        jsonVec.append(vec[i]);
    }  

    return jsonVec;
}

Json::Value ResultWriterJson::writeLine(arma::cx_drowvec vec)
{
    Json::Value jsonVecR(Json::arrayValue);
    Json::Value jsonVecI(Json::arrayValue);
    for (uint i = 0; i < vec.size(); i++)
    {
        jsonVecR.append(vec[i].real());
        jsonVecI.append(vec[i].imag());
    }  

    Json::Value jsonVecC(Json::objectValue);
    jsonVecC["r"] = jsonVecR;
    jsonVecC["i"] = jsonVecI;

    return jsonVecC;
}

Json::Value ResultWriterJson::writeLineInReal(arma::cx_drowvec vec)
{
    Json::Value jsonVec(Json::arrayValue);
    for (uint i = 0; i < vec.size(); i++)
    {
        jsonVec.append(vec[i].real());
    }  

    return jsonVec;
}




ResultWriterJson::~ResultWriterJson(){}