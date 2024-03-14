#pragma once

#include <bitset>
#include <vector>


#define MAX_BITS 32



using namespace std;

bitset<MAX_BITS> bitsetFromVector(const vector<bool>& vec);

void bitsetswap(bitset<MAX_BITS>& bs1, bitset<MAX_BITS>& bs2);

void fusionSortVector(vector<bitset<MAX_BITS>>& vec);

bitset<MAX_BITS> positionToBinary(vector<uint>&);

vector<uint> binaryToPosition(bitset<MAX_BITS>&);

uint getParticleNumberBetween(std::bitset<MAX_BITS>& posP, uint from, uint to);

uint countDoublons(bitset<MAX_BITS>& phiH1, bitset<MAX_BITS>& phiH2);

uint countSimilarities(bitset<MAX_BITS>& phiH1, bitset<MAX_BITS>& phiH2, uint nSites, uint& lastIndexWrong1, uint& lastIndexWrong2);

bitset<MAX_BITS> expandBitSetUnderConstraint(bitset<MAX_BITS> reducedBS, bitset<MAX_BITS> constraintBS, uint sizeReducedBS = MAX_BITS);


