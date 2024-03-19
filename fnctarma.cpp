
#include "fnctarma.h"



/// @brief Convert a vector of vector of double into the corresponding arma::dmat. The vector'vector and the matrix are seen in row-major order (reading/writing order).
/// @param vecVec A vector of vector of double 
/// @return The arma::dmat version of the vector
arma::dmat armaMatfromVectorsVector(std::vector<std::vector<double>>& vecVec)
{
    arma::dmat mattrix = arma::dmat();
    if (vecVec.size() > 0)
    {
        if (vecVec[0].size() > 0)
        {
            mattrix = arma::dmat(vecVec.size(), vecVec[0].size(), arma::fill::zeros);
            for (uint iRow = 0; iRow < vecVec.size(); iRow ++)
            {
                for (uint iCol = 0; iCol < std::min(vecVec[iRow].size(), vecVec[0].size()); iCol ++)
                {
                    mattrix(iRow, iCol) = vecVec[iRow][iCol];
                }
            }
        }
    }

    return mattrix;
}