#include "../include/MapDoubleEntry.h"

template <typename T> 
MapDoubleEntry<T>::MapDoubleEntry()
{
    map = std::map<T, std::string>();
    mapReverse = std::map<std::string, T>();
}

template <typename T> 
void MapDoubleEntry<T>::append(std::string key, T value)
{
    map[key] = value;
    mapReverse[value] = key;
}

template <typename T> 
std::string MapDoubleEntry<T>::get(T value)
{
    return map.at(value);
}

template <typename T> 
T MapDoubleEntry<T>::get(std::string key)
{
    return map.at(key);
}

template <typename T> 
MapDoubleEntry<T>::~MapDoubleEntry<T>()
{
    delete map;
    delete mapReverse;
}