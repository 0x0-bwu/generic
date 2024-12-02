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
    explicit LinearMap(size_t n) : m_data(n) {}
    explicit LinearMap(size_t n, Value value) : m_data(n, value) {}
    explicit LinearMap(std::vector<Value> && values) : m_data(values) {}

    Iterator Begin() { return m_data.begin(); }

    Iterator End() { return m_data.end(); }

    ConstIterator Begin() const { return m_data.begin(); }

    ConstIterator End() const { return m_data.end(); }

    ConstReverseIterator ReverseBegin() const { return m_data.rbegin(); }

    ConstReverseIterator ReverseEnd() const { return m_data.rend(); }

    Reference operator[] (const Key key)
    { 
        GENERIC_ASSERT_MSG(size_t(key) < m_data.size(), "index out of range");
        return m_data[size_t(key)]; 
    }   

    ConstReference operator[] (const Key key) const
    { 
        GENERIC_ASSERT_MSG(size_t(key) < m_data.size(), "index out of range");
        return m_data[size_t(key)]; 
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
    void Append(Args &&... args) { m_data.emplace_back(std::forward<Args>(args)...); }

    template<typename... Args>
    void Resize(Args &&... args) { m_data.resize(std::forward<Args>(args)...); }

    void Clear() { m_data.clear(); }

    size_t Capacity() const { return m_data.capacity(); }

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
    void swap(LinearMap<Key,Value> & other) { return Swap(other); }
    void clear() { return Clear(); }
    template<typename... Args>
    void emplace_back(Args &&... args) { Append(std::forward<Args>(args)...); }
    Iterator begin() { return Begin(); }
    Iterator end() { return End(); }

protected:
    std::vector<Value> m_data;
};

} // namespace generic::utils