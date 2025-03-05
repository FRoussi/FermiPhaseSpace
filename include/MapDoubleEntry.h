#pragma once

#include <vector>
#include <armadillo>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>



template <typename T> class MapDoubleEntry
{

public :
    MapDoubleEntry();

    void append(std::string key, T value);

    std::string get(T value);
    T get(std::string key);

    ~MapDoubleEntry<T>();

private :
    std::map<T, std::string>    map;
    std::map<std::string, T>    mapReverse;

};