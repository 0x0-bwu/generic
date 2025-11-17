/**
 * @file LinearMap.hpp
 * @author bwu
 * @brief linear map utils implement based on https://github.com/verilog-to-routing/tatum  
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <vector>

namespace generic::utils {

template<typename Key, typename Value>
class LinearMap
{
public:
    using Iterator = typename std::vector<Value>::iterator;
    using ConstIterator = typename std::vector<Value>::const_iterator;
    using ConstReverseIterator = typename std::vector<Value>::const_reverse_iterator;
    
    using Reference = typename std::vector<Value>::reference;
    using ConstReference = typename std::vector<Value>::const_reference;
public:
    LinearMap() = default;
    virtual ~LinearMap() = default;
/**
 * @brief Brief description of LinearMap.
 * @param m_data(n
 * @return explicit
 */
    explicit LinearMap(size_t n) : m_data(n) {}
/**
 * @brief Brief description of LinearMap.
 * @param n
 * @param m_data(n
 * @param value
 * @return explicit
 */
    explicit LinearMap(size_t n, Value value) : m_data(n, value) {}
/**
 * @brief Brief description of LinearMap.
 * @param m_data(values
 * @return explicit
 */
    explicit LinearMap(std::vector<Value> && values) : m_data(values) {}

/**
 * @brief Brief description of Begin.
 * @return Iterator
 */
    Iterator Begin() { return m_data.begin(); }

/**
 * @brief Brief description of End.
 * @return Iterator
 */
    Iterator End() { return m_data.end(); }

    ConstIterator Begin() const { return m_data.begin(); }

    ConstIterator End() const { return m_data.end(); }

    ConstReverseIterator ReverseBegin() const { return m_data.rbegin(); }

    ConstReverseIterator ReverseEnd() const { return m_data.rend(); }

    Reference operator[] (const Key key)
    { 
        GENERIC_ASSERT_MSG(size_t(key) < m_data.size(), "index out of range");
        return this->operator[] (size_t(key));
    }   

    ConstReference operator[] (const Key key) const
    { 
        GENERIC_ASSERT_MSG(size_t(key) < m_data.size(), "index out of range");
        return this->operator[] (size_t(key));
    }

    Reference operator[] (size_t i)
    { 
        GENERIC_ASSERT_MSG(i < m_data.size(), "index out of range");
        return m_data[i]; 
    }

    ConstReference operator[] (size_t i) const
    { 
        GENERIC_ASSERT_MSG(i < m_data.size(), "index out of range");
        return m_data[i]; 
    }

    Iterator Find(const Key key)
    {
        if(size_t(key) < m_data.size())
            return std::advance(m_data.begin(), size_t(key));
        return m_data.end();
    }

    ConstIterator Find(const Key key) const
    {
        if(size_t(key) < m_data.size())
            return std::advance(m_data.begin(), size_t(key));
        return m_data.cend();
    }

    size_t Size() const { return m_data.size(); }

    bool Empty() const { return m_data.empty(); }

    bool Contain(const Key key) const { return size_t(key) < m_data.size(); }

    template<typename... Args>
/**
 * @brief Brief description of Append.
 * @param args
 * @return Key
 */
    Key Append(Args &&... args) { m_data.emplace_back(std::forward<Args>(args)...); return Key(Size() - 1); }

    template<typename... Args>
/**
 * @brief Brief description of Resize.
 * @param args
 * @return void
 */
    void Resize(Args &&... args) { m_data.resize(std::forward<Args>(args)...); }

/**
 * @brief Brief description of Clear.
 * @return void
 */
    void Clear() { m_data.clear(); }

    size_t Capacity() const { return m_data.capacity(); }

/**
 * @brief Brief description of Reserve.
 * @param size
 * @return void
 */
    void Reserve(size_t size) { return m_data.reserve(size); }

/**
 * @brief Brief description of Shrink.
 * @return void
 */
    void Shrink() { m_data.shrink_to_fit(); }

    void Insert(const Key key, const Value & value)
    {
        if(size_t(key) >= m_data.size())
            m_data.resize(size_t(key) + 1);
        m_data[size_t(key)] = value;
    }

    void Insert(const Key key, Value && value)
    {
        if(size_t(key) >= m_data.size())
            m_data.resize(size_t(key) + 1);
        m_data[size_t(key)] = value;
    }

    void Swap(LinearMap<Key,Value> & other)
    {
        std::swap(m_data, other.m_data);
    }

    // std-style interface
    size_t size() const { return Size(); }
/**
 * @brief Brief description of swap.
 * @param LinearMap<Key
 * @param other
 * @return void
 */
    void swap(LinearMap<Key,Value> & other) { return Swap(other); }
/**
 * @brief Brief description of clear.
 * @return void
 */
    void clear() { return Clear(); }
/**
 * @brief Brief description of reserve.
 * @param size
 * @return void
 */
    void reserve(size_t size) { return Reserve(size); }

    template<typename... Args>
/**
 * @brief Brief description of emplace_back.
 * @param args
 * @return Key
 */
    Key emplace_back(Args &&... args) { return Append(std::forward<Args>(args)...); }
/**
 * @brief Brief description of begin.
 * @return Iterator
 */
    Iterator begin() { return Begin(); }
/**
 * @brief Brief description of end.
 * @return Iterator
 */
    Iterator end() { return End(); }
    ConstIterator begin() const { return Begin(); }
    ConstIterator end() const { return End(); }
    ConstIterator cbegin() const { return Begin(); }
    ConstIterator cend() const { return End(); }

#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & boost::serialization::make_nvp("data", m_data); 
    }
#endif//GENERIC_BOOST_SERIALIZATION_SUPPORT

protected:
    std::vector<Value> m_data;
};

} // namespace generic::utils
