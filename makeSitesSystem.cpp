//#include <stdio.h>
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

bool isIn(string word, vector<string>& wordList)
{
    bool isInB = false;
    vector<string>::iterator itWord = wordList.begin();
    while ( (!isInB) && (itWord != wordList.end()) )
    {
        if (*itWord == word)  isInB = true;
        itWord ++;
    }

    return isInB;
}


int main(int argc, char* argv[])
{



   cout << "Writing a new system." << endl;


    string fileName = "newSystem";
    string fileFolder = "systemsFolder";

// Default values
    vector<vector<double>> pos = makeGrid(2,1,1,1);
    string geometry = "";
    string IC = "ic_sw";

// Possible values
    vector<string> possibleGeometry = {"grid", "graphene"};
    vector<string> possibleIC = {"ic_sw", "ic_cw"};

    bool readGeomParam = false;                     //<- are we reading the geometry parameters
    vector<double> geomParam = vector<double>();    //<- list of geometry parameters stored in double
    double argDouble = 0;                           //<- the current entry in double
    bool isANumber = true;                          //<- is the current entry a number

    for (int iArg = 1; iArg < argc ; iArg++)        //<- for each entry
    {

        isANumber = true;
        try                                         //<- We check if it's a number with a try/catch block
        {
            argDouble = stod(argv[iArg]);
        }  catch(const std::invalid_argument& ia) 
        {
            isANumber = false;
        }

        if (isANumber && readGeomParam)             //<- if it's a number and we are readinf the geometruy parameters
        {
            geomParam.push_back(argDouble);         //<- we store the number in the goemetry parameter vector
        }
        else if (isIn(string(argv[iArg]), possibleGeometry)) //((argv[iArg] == string("grid")) || (argv[iArg] == string("graphene")))  //<- else, if the word is among the geometry types
        {
            geometry = argv[iArg];                  //<- we note the geometry
            readGeomParam = true;                   //<- We stand ready to read the geometry parameters
        } 
        else if (isIn(string(argv[iArg]), possibleIC)) //((argv[iArg]==string("ic_sw")) || (argv[iArg]==string("ic_cw")))    //<- else, if the word is among the Initial conditions types
        {
            IC = argv[iArg] ;                       //<- we note the initial condition
            readGeomParam = false;                  //<- We stop reading geometry parameters
        }
        else 
        {
             readGeomParam = false;                 //<- Else we don't stand ready for anything
             if (iArg==1)                           //<- and if it's the first argumentm that's the name of the file
             {
                fileName = argv[1];
             }
        }

    }

// Computation of the geometry
for(unsigned int iGeom = geomParam.size(); iGeom<4; iGeom++)
{
    geomParam.push_back(1.);
}

    if (geometry == "grid")
    {
        pos = makeGrid( static_cast<int>(geomParam[0]), static_cast<int>(geomParam[1]), static_cast<int>(geomParam[2]), geomParam[3]);
    }
    if (geometry == "graphene")
    {
        pos = makeHoneyCombs(static_cast<int>(geomParam[0]),geomParam[1]);
    }

// Computation of the Initial condition
vector<unsigned int> IC_UP = vector<unsigned int>();
vector<unsigned int> IC_DOWN = vector<unsigned int>();
if (IC == "ic_sw")
{
    for (unsigned int uppos = 0; uppos < pos.size(); uppos+=2)
    {
            IC_UP.push_back(uppos);
    }
    for (unsigned int downpos = 1; downpos < pos.size(); downpos+=2)
    {
            IC_DOWN.push_back(downpos);
    }
}
else if (IC == "ic_cw")
{
    for (unsigned int uppos = 0; uppos < pos.size(); uppos+=2)
    {
            IC_UP.push_back(uppos);
    }
    for (unsigned int downpos = 0; downpos < pos.size(); downpos+=2)
    {
            IC_DOWN.push_back(downpos);
    }
}




// Writing of everything on the file :
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


        streamSystem << "#INITIALCONDITIONS_GIVENDIAGONAL" << endl;
        for (vector<unsigned int>::iterator itpos = IC_UP.begin(); itpos !=IC_UP.end(); itpos++)
        {
            streamSystem << (*itpos) << "\t";
        }
        streamSystem << endl;
        for (vector<unsigned int>::iterator itpos = IC_DOWN.begin(); itpos !=IC_DOWN.end(); itpos++)
        {
            streamSystem << (*itpos) << "\t";
        }
        streamSystem << endl;
        

        streamSystem.close();
    }
    else
    {
        cout << "Stream not open." << endl;
    }


    cout << fileName << " system created with " << pos.size() << " sites and " << IC_UP.size() + IC_DOWN.size() << " particles." << endl;  



    return 0;
}