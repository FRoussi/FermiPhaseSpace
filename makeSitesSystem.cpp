//#include <stdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;

using namespace std;

/// @brief The suite of integer doubles S defined by  S(t+1) = getNextIndex(S(t))
/// and starting by (0,0) paves the positive plan in an ordely fashion. The
/// aim is for [S(0) .. S(t)] to stay "as square as possible" for every index t.  
void getNextIndex(unsigned int& i, unsigned int& j)
{
    if (i==j)
    {
        i=i+1;
        j = 0;
    }
    else if (i==(j+1))
    {
        j = i;
        i = 0;
    }
    else if (i>j)
    {
        j = j+1;
    }
    else
    {
        i=i+1;
    }

}

vector<vector<double>> makeGrid(unsigned int nx, unsigned int ny, unsigned int nz, double a)
{
    
    unsigned int nPoints = nx*ny*nz;
    vector<vector<double>> pos = vector<vector<double>>(nPoints);

    unsigned int counter = 0;
    for (unsigned int z = 0; z < nz; z++)
    {
        for (unsigned int y = 0; y < ny; y++)
        {
            for (unsigned int x = 0; x < nx; x++)
            {
                vector<double> newPos = {a*x, a*y, a*z};
                pos[counter] = newPos;
                counter +=1;
            }
        }
    }

    return pos;
}





vector<vector<double>> makeHoneyCombs(unsigned int nCombs, double a)
{

    double sqrt3 = std::sqrt(3);
unsigned int nOfPoints = 0;
    vector<vector<double>> pos = vector<vector<double>>();

vector<vector<double>> pos0 = { {-sqrt3/2 *a, -0.5 *a},
                                {-sqrt3/2 *a, 0.5 *a},
                                {0   *a,      1   *a},
                                {0   *a,      -1  *a},
                                {sqrt3/2 *a, -0.5 *a},
                                {sqrt3/2 *a, 0.5  *a}};
vector<vector<double>> posNew = pos0;


double iVec[2] = {sqrt3 *a, 0 *a};
double jVec[2] = {sqrt3/2 *a, 1.5 *a};

unsigned int basisIndex[2] = {0,0};


for (unsigned int bv = 0; bv < nCombs; bv++)
{
    for (unsigned int p = 0; p<6; p++)
    {
        posNew[p][0] = pos0[p][0] + iVec[0] * basisIndex[0] + jVec[0] * basisIndex[1];
        posNew[p][1] = pos0[p][1] + iVec[1] * basisIndex[0] + jVec[1] * basisIndex[1];;
    }
    pos.insert(pos.end(), posNew.begin(), posNew.end());

    getNextIndex(basisIndex[0], basisIndex[1]);

    nOfPoints +=6;

}

// We remove the doublons :
double samePointThreshold = 0.001 * a;


    for (unsigned int i1 = 0; i1< pos.size(); i1++)         //<- we look for all points equal to positions(i1)
    {
        unsigned int nDoublonsFound = 0;                    //<- There can't be more than three times the same point
         
        unsigned int i2 = i1+1;
        while ((i2<pos.size()) && (nDoublonsFound < 2))    //<- as long as there are points and more doublons to find
        {
            double d = 0;
            for (unsigned int iDim = 0; iDim<pos[i1].size(); iDim++)
            {
                d += std::abs(pos[i1][iDim] - pos[i2][iDim]);
            }
            if (abs(d) < samePointThreshold)                //<- if i1 and i2 are the same point :
            {
                pos.erase(pos.begin() + i2);                //<- We remove the point i2
                nDoublonsFound ++;                          //<- We say we found one more doublon
            }
            else
            {
                i2 ++;                                      //<- else, we look at the next point
            }
                
        }

    }
    
    pos.shrink_to_fit();
    return pos;
}

int main(int argc, char* argv[])
{



   cout << "SYSTEEEEEM OFADOOOOOWN" << endl;


    string fileName = "newSystem.txt";
    string fileFolder = "systemsFolder";


   vector<vector<double>> pos = makeGrid(2,1,1,1);

   if (argc > 1)
   {
    string fName = argv[1];
    fileName = fName;
   } 
   if (argc >2)
   {
        string geometry = argv[2];
        if (geometry == "grid")
        {
            if (argc == 3)  pos = makeGrid(2,1,1,1);
            if (argc == 4)  pos = makeGrid(stoi(argv[3]),1,1,1);
            if (argc == 5)  pos = makeGrid(stoi(argv[3]),stoi(argv[4]),1,1);
            if (argc == 6)  pos = makeGrid(stoi(argv[3]),stoi(argv[4]),stoi(argv[5]),1);
            if (argc == 7)  pos = makeGrid(stoi(argv[3]),stoi(argv[4]),stoi(argv[5]),stod(argv[6]));
        }
         if (geometry == "graphene")
        {
            if (argc == 3)  pos = makeHoneyCombs(1,1);
            if (argc == 4)  pos = makeHoneyCombs(stoi(argv[3]),1);
            if (argc == 5)  pos = makeHoneyCombs(stoi(argv[3]),stoi(argv[4]));
        }
   }


string fileAdress = fileFolder + "/" + fileName;


   ofstream streamSystem(fileAdress);
   
    if (streamSystem.is_open())
    {
       streamSystem << "#POSITIONS" << endl;
        for (vector<vector<double>>::iterator itpos = pos.begin(); itpos !=pos.end(); itpos++)
        {
           for(vector<double>::iterator itpoint = itpos->begin() ; itpoint != itpos->end(); itpoint++)
           {
                streamSystem << (*itpoint) << "\t";
           }
            streamSystem << endl;
        }


        streamSystem << "#INITIALCONDITIONS" << endl;
        for (unsigned int uppos = 0; uppos < pos.size(); uppos+=2)
        {
            streamSystem << uppos << "\t";
        }
        streamSystem << endl;
        for (unsigned int downpos = 1; downpos < pos.size(); downpos+=2)
        {
            streamSystem << downpos << "\t";
        }
        streamSystem << endl;
        

    }
    else
    {
        printf("Stream not open.");
    }






    return 0;
}