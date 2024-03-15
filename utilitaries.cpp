#include "utilitaries.h"


/// @brief Check if a vector contains an element.
/// @param w the element to look for
/// @param v the vector of possibilities
/// @return true if w is in v
template <typename T>       bool isIn(T w, std::vector<T>& v)
{

        bool isInB = false;
        unsigned int index = 0;
        while ( (!isInB) && (index < v.size()) )
        {
            if (v[index] == word)  isInB = true;
            index ++;
        }

        return isInB;
}

/// @brief Check if a vector contains an element.
/// @param w the element to look for
/// @param v the vector of possibilities
/// @param firstIndex index of the first found element
/// @return true if w is in v
template <typename T>       bool isIn(T w, std::vector<T>& v, uint& firstIndex)
{
    bool isInB = false;
    unsigned int index = 0;
    while ( (!isInB) && (index < v.size()) )
    {
        if (v[index] == word)  isInB = true;
        index ++;
    }
    firstIndex = index-1;
    return isInB;
}