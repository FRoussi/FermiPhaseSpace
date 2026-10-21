#pragma once

#include <string>
#include <map>
#include <stdexcept>


// Bidirectional map : behaves like a std::map<Key, Value> whose keys can also
// be accessed from the values. The two internal maps are kept in sync by
// insert, erase and clear; every lookup is read-only.
// - operator[] returns a default-constructed object if the entry is unknown,
//   use contains / containsValue to test for existence.
// - at throws std::out_of_range if the entry is unknown, like std::map::at.
// - if the same value is inserted for two different keys, the value maps back
//   to the last key inserted.
template <typename Key, typename Value = std::string>
class MapDoubleEntry
{

public :

    typedef typename std::map<Key, Value>::const_iterator const_iterator;


    // modifiers : the only way to add or remove entries, both maps are updated
    void                            insert(Key key, Value value);
    std::size_t                     erase(const Key& key);
    void                            clear();

    // capacity
    std::size_t                     size() const;
    bool                            empty() const;

    // lookup, key -> value
    const Value&                    getValue(const Key& key) const;
    const Value&                    operator[](const Key& key) const;
    const Value&                    at(const Key& key) const;
    bool                            contains(const Key& key) const;
    const_iterator                  find(const Key& key) const;

    // lookup, value -> key
    const Key&                      getKey(const Value& value) const;
    //const Key&                      operator[](const Value& value) const;
    const Key&                      atKey(const Value& value) const;
    bool                            containsValue(const Value& value) const;

    // iteration over the key -> value pairs, in std::map order
    const_iterator                  begin() const;
    const_iterator                  end() const;


private :

    std::map<Key, Value>            map;
    std::map<Value, Key>            mapReverse;

};



template <typename Key, typename Value>
void                        MapDoubleEntry<Key, Value>::insert(Key key, Value value)
{
    typename std::map<Key, Value>::iterator it = map.find(key);
    if (it != map.end())
    {
        mapReverse.erase(it->second);
    }

    map[key] = value;
    mapReverse[value] = key;
}

template <typename Key, typename Value>
std::size_t                 MapDoubleEntry<Key, Value>::erase(const Key& key)
{
    typename std::map<Key, Value>::iterator it = map.find(key);
    if (it == map.end())
    {
        return 0;
    }

    mapReverse.erase(it->second);
    map.erase(it);
    return 1;
}

template <typename Key, typename Value>
void                        MapDoubleEntry<Key, Value>::clear()
{
    map.clear();
    mapReverse.clear();
}

template <typename Key, typename Value>
std::size_t                 MapDoubleEntry<Key, Value>::size() const
{
    return map.size();
}

template <typename Key, typename Value>
bool                        MapDoubleEntry<Key, Value>::empty() const
{
    return map.empty();
}

template <typename Key, typename Value>
const Value&                MapDoubleEntry<Key, Value>::getValue(const Key& key) const
{
    static const Value defaultValue = Value();

    typename std::map<Key, Value>::const_iterator it = map.find(key);
    if (it != map.end())
    {
        return it->second;
    }
    return defaultValue;
}

template <typename Key, typename Value>
const Value&                MapDoubleEntry<Key, Value>::operator[](const Key& key) const
{
    return getValue(key);
}

template <typename Key, typename Value>
const Value&                MapDoubleEntry<Key, Value>::at(const Key& key) const
{
    return map.at(key);
}

template <typename Key, typename Value>
bool                        MapDoubleEntry<Key, Value>::contains(const Key& key) const
{
    return (map.find(key) != map.end());
}

template <typename Key, typename Value>
typename MapDoubleEntry<Key, Value>::const_iterator   MapDoubleEntry<Key, Value>::find(const Key& key) const
{
    return map.find(key);
}

template <typename Key, typename Value>
const Key&                  MapDoubleEntry<Key, Value>::getKey(const Value& value) const
{
    static const Key defaultKey = Key();

    typename std::map<Value, Key>::const_iterator it = mapReverse.find(value);
    if (it != mapReverse.end())
    {
        return it->second;
    }
    return defaultKey;
}

template <typename Key, typename Value>
const Key&                  MapDoubleEntry<Key, Value>::atKey(const Value& value) const
{
    return mapReverse.at(value);
}

template <typename Key, typename Value>
bool                        MapDoubleEntry<Key, Value>::containsValue(const Value& value) const
{
    return (mapReverse.find(value) != mapReverse.end());
}

template <typename Key, typename Value>
typename MapDoubleEntry<Key, Value>::const_iterator   MapDoubleEntry<Key, Value>::begin() const
{
    return map.begin();
}

template <typename Key, typename Value>
typename MapDoubleEntry<Key, Value>::const_iterator   MapDoubleEntry<Key, Value>::end() const
{
    return map.end();
}

