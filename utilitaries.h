#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>


#define M_PI 3.14159265358979323846


/// @brief Check if a vector contains an element.
/// @param w the element to look for
/// @param v the vector of possibilities
/// @return true if w is in v
template <typename T>       bool isIn(T w, std::vector<T>& v)
{
    std::find(v.begin(), v.end(), w) != v.end();
    /*
        bool isInB = false;
        vector<T>::iterator itWord = wordList.begin();
        while ( (!isInB) && (itWord != wordList.end()) )
        {
            if (*itWord == word)  isInB = true;
            itWord ++;
        }

        return isInB;*/
}