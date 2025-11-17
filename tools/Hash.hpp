/**
 * @file Hash.hpp
 * @author bwu
 * @brief hash related functions
 * @version 0.1
 * @date 2024-11-15
 */
#pragma once
#include <boost/functional/hash.hpp>
#include <iostream>
#include <numeric>

namespace generic::hash {

template <typename T>
inline size_t Hash(const T & v)
{
    return std::hash<T>()(v);
}

template <typename ... Args>
inline size_t HashCombine(size_t seed, Args &&... args)
{
    (boost::hash_combine(seed, args), ...); return seed;
}

template <typename Container, typename Hasher = std::hash<typename Container::value_type>>
inline size_t OrderedHash(const Container & c, const Hasher & hasher = Hasher())
{
    return std::accumulate(c.begin(), c.end(), size_t{0}, [&](size_t seed, const auto & v) { return HashCombine(seed, hasher(v)); });
}

template <typename Container, typename Hasher = std::hash<typename Container::value_type>>
inline size_t UnorderedHash(const Container & c, const Hasher & hasher = Hasher())
{
    return std::accumulate(c.begin(), c.end(), size_t{0}, [&](size_t seed, const auto & v) { return seed ^ hasher(v); });
}



} // namespace generic::hash
