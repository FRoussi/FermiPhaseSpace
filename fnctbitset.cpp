#include "fnctbitset.h"



/// @brief Convert a vector of boolean into a bitset of boolean 
/// @param vec vector of boolean to convert
/// @return the converted bitset 
bitset<MAX_BITS> bitsetFromVector(const vector<bool>& vec)
{
    bitset<MAX_BITS> bs = bitset<MAX_BITS>();
    uint maxV = min((int)vec.size(), MAX_BITS);
    for (uint it = 0; it<maxV; it++)
    {
        if (vec[it])  bs.set(it);
    }

    return bs;
}


/// @brief Swap two bitset 
/// @param bs1 bitset
/// @param bs2 bitset
/// @return void 
void bitsetswap(bitset<MAX_BITS>& bs1, bitset<MAX_BITS>& bs2)
{
   bs1 = bs1 ^ bs2;
   bs2 = bs2 ^ bs1;
   bs1 = bs1 ^ bs2;
}


/// @brief Sorting of a vector of bistset with the tri-fusion algorithm 
/// @param vec vector of particle positions, in bitset
/// @return void 
void fusionSortVector(vector<bitset<MAX_BITS>>& vec)
{
    uint vecSize = vec.size();

    if (vecSize == 1);

    else if (vecSize == 2)
    {
        if (vec[0].to_ulong() > vec[1].to_ulong())
        {
            bitsetswap(vec[0], vec[1]);
        }
        
    }

    else
    {

        uint mid = (uint)((vecSize / 2));

        vector<bitset<MAX_BITS>> vec1(vec.begin(), vec.begin() + mid);
        fusionSortVector(vec1);

        vector<bitset<MAX_BITS>> vec2(vec.begin() + mid + 1, vec.end());
        fusionSortVector(vec2);

        uint countor1 = 0;
        uint countor2 = 0;

        while (countor1 + countor2 < vecSize)
        {
            if (countor1 == vec1.size())
            {
                while (countor2 < vec2.size())
                {
                    vec[countor1 + countor2] = vec2[countor2];
                    countor2 ++;
                }
            }
            else if (countor2 == vec2.size())
            {
                while (countor1 < vec1.size())
                {
                    vec[countor1 + countor2] = vec1[countor1];
                    countor1 ++;
                }
            }
            else if (vec1[countor1].to_ulong() < vec2[countor2].to_ulong())
            {
                vec[countor1 + countor2] = vec1[countor1];
                countor1++;
            }
            else
            {
                vec[countor1 + countor2] = vec2[countor2];
                countor2++;
            }
        }

    }
}



/// @brief Convert a position vector into an empty/filled vector.
/// @param pos vector of particle positions
/// @return Boolean vector of empty/filled sites 
std::bitset<MAX_BITS> positionToBinary(std::vector<uint>& pos)
{
    std::bitset<MAX_BITS> bs = std::bitset<MAX_BITS>();
    for (std::vector<uint>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        bs.set(*it);
    }
    return bs;
}


/// @brief Convert an empty/filled vector into a position vector.
/// @param binaryPos a bitset of particle presence
/// @return Position vector 
std::vector<uint> binaryToPosition(std::bitset<MAX_BITS>& binaryPos)
{
    std::vector<uint> vec = std::vector<uint>();
    for (uint it = 0; it < binaryPos.size(); it++)
    {
        if (binaryPos[it])
        {
            vec.push_back(it);
        }
    }
    return vec;
}





/// @brief Count the similarities between 2 boolean half vectors and store the last difference.
/// @param phiH1 wave-function half basis bitset
/// @param phiH2 wave-function half basis bitset
/// @param lastIndexWrong1 If there is a difference, update it to be the last difference with a particle of phiH1
/// @param lastIndexWrong2 If there is a difference, update it to be the last difference with a particle of phiH2
/// @return Number of similarities between the 2 wave-function basis vector.
uint countSimilarities(bitset<MAX_BITS>& phiH1, bitset<MAX_BITS>& phiH2, uint nSites, uint& lastIndexWrong1, uint& lastIndexWrong2)
{
    uint simCount = 0;
    uint lastIdW1 = 0;
    uint lastIdW2 = 0;

    for (uint l = 0; l < nSites; l++)
    {
        if (phiH1[l] == phiH2[l])
            simCount += 1;
        else
        {
            if (phiH1[l] == true)
                    lastIdW1 = l;
            else if (phiH2[l] == true)
                    lastIdW2 = l;
        }
    }

    lastIndexWrong1 = lastIdW1;
    lastIndexWrong2 = lastIdW2;
    return simCount;
}


/// @brief Count the doublons in 2 bitsets.
/// @param phiH1 wave-function half basis bitset
/// @param phiH2 wave-function half basis bitset
/// @return Number of doublons
uint countDoublons(bitset<MAX_BITS>& phiH1, bitset<MAX_BITS>& phiH2)
{
    bitset<MAX_BITS> andBits = phiH1 & phiH2;
    return andBits.count();
}


/// @brief Count the number of particles in a wave-function basis vector between 2 sites
/// @param posP empty/filled vector
/// @param from index of 1st site
/// @param to index of 2nd site
/// @return Number of particles between 'from' and 'to'
uint getParticleNumberBetween(std::bitset<MAX_BITS>& posP, uint from, uint to)
{
    uint nPart = 0;
    if ((from < MAX_BITS) && (to < MAX_BITS))
    {
        uint f = std::min(from, to);
        uint t = std::max(from, to);
        for (uint p = f + 1; p < t; p++)
        {
            if (posP[p])
                nPart++;
        }
    }

    return nPart;
}