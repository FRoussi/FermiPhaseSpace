
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

        bool isInB = false;
        unsigned int index = 0;
        while ( (!isInB) && (index < v.size()) )
        {
            if (v[index] == w)  isInB = true;
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
        if (v[index] == w)  isInB = true;
        index ++;
    }
    firstIndex = index-1;
    return isInB;
}