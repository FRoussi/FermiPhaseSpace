#include "ResultWriter.h"



ResultWriter::ResultWriter(const std::string& fileName, const std::string& folderName)
{

    resultFileName = fileName + ".txt";
    resultFolder = folderName;
    fileAdressName = resultFileName;
    if (resultFolder != "")
    {
        fileAdressName = resultFolder + "/" + fileAdressName;
    }

}

void ResultWriter::newResults()
{
    std::ofstream resultFile(fileAdressName);
    resultFile.close();
}


void ResultWriter::writeLine(std::string lineName, std::vector<double> vec)
{
    std::ofstream resultFile(fileAdressName, std::ios::ate| std::ios::app);
    if (resultFile.is_open())
    {
        resultFile << lineName << " ";
        for (std::vector<double>::iterator itv = vec.begin(); itv !=vec.end(); itv++)
        {
            resultFile << std::to_string(*itv) << " ";
        }
        resultFile << std::endl;
    }
    else
    {
        printf("Stream not open.");
    }
    resultFile.close();
}

void ResultWriter::writeLine(std::string lineName, std::vector<unsigned int> vec)
{
    std::ofstream resultFile(fileAdressName, std::ios::ate | std::ios::app);
    if (resultFile.is_open())
    {
        resultFile << lineName << " ";
        for (std::vector<unsigned int>::iterator itv = vec.begin(); itv !=vec.end(); itv++)
        {
            resultFile << std::to_string(*itv) << " ";
        }
        resultFile << std::endl;
    }
    else
    {
        printf("Stream not open.");
    }
    resultFile.close();
}

void ResultWriter::writeLine(std::string lineName, arma::drowvec vec)
{
    std::ofstream resultFile(fileAdressName, std::ios::ate| std::ios::app);
    if (resultFile.is_open())
    {
        resultFile << lineName << " ";
        for (arma::drowvec::iterator itv = vec.begin(); itv != vec.end(); itv++)
        {
           resultFile << std::to_string(*itv) << " ";
        }
        resultFile << std::endl;
    }
    else
    {
        printf("Stream not open.");
    }
    resultFile.close();
}

void ResultWriter::writeLine(std::string lineName, arma::cx_drowvec vec)
{
    std::ofstream resultFile(fileAdressName, std::ios::ate| std::ios::app);
    if (resultFile.is_open())
    {
        resultFile << lineName << " ";
        for (arma::cx_drowvec::iterator itv = vec.begin(); itv != vec.end(); itv++)
        {
           resultFile << std::to_string(std::real(*itv)) << "+i" << std::to_string(std::imag(*itv)) << " ";
        }
        resultFile << std::endl;
    }
    else
    {
        printf("Stream not open.");
    }
    resultFile.close(); 
}

void ResultWriter::writeLineInReal(std::string lineName, arma::cx_drowvec vec)
{
    std::ofstream resultFile(fileAdressName, std::ios::ate| std::ios::app);
    if (resultFile.is_open())
    {
        resultFile << lineName << " ";
        for (arma::cx_drowvec::iterator itv = vec.begin(); itv != vec.end(); itv++)
        {
           resultFile << std::to_string(std::real(*itv)) << " ";
        }
        resultFile << std::endl;
    }
    else
    {
        printf("Stream not open.");
    }
    resultFile.close();
}




ResultWriter::~ResultWriter(){}