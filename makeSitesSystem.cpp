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

vector<vector<double>> pos0 = { {-sqrt3/2, -0.5},
                                {-sqrt3/2, 0.5},
                                {0  ,      1  },
                                {0  ,      -1 },
                                {sqrt3/2, -0.5},
                                {sqrt3/2, 0.5 }};
vector<vector<double>> posNew = pos0;


double iVec[2] = {sqrt3, 0};
double jVec[2] = {sqrt3/2, 1.5};

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


    pos.shrink_to_fit();
    return pos;
}

int main(int argc, char* argv[])
{



   cout << "SYSTEEEEEM OFADOOOOOWN" << endl;


    string fileName = "newSystem.txt";
    string fileFolder = "systemsFolder";
    string fileAdress = fileFolder + "/" + fileName;

   vector<vector<double>> pos = makeGrid(3,1,1,1);

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