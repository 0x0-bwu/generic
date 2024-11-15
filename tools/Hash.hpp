/**
 * @file Hash.hpp
 * @author bwu
 * @brief hash related functions
 * @version 0.1
 * @date 2024-11-15
 */

#include <boost/functional/hash.hpp>
#include <iostream>
namespace generic::hash {

template <typename Key>
inline size_t Hash(const Key & key)
{
    return std::hash<Key>{}();
}

template <typename ... Args>
inline size_t HashCombine(Args &&... args)
{
    size_t seed{0}; (boost::hash_combine(seed, args), ...); return seed;
}

} // namespace generic::hash