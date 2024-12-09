/**
 * @file Index.hpp
 * @author bwu
 * @brief index utils implement based on https://github.com/verilog-to-routing/tatum
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <limits>

namespace generic::utils {

template<typename Tag, typename T = size_t>
class Index
{
	static_assert(std::is_integral_v<T> and std::is_unsigned_v<T>/*shoule be unsigned integral type*/);
public:
    using SizeType = T;
    static constexpr T INVALID_ID = std::numeric_limits<T>::max();
    
    constexpr Index() = default;
    virtual ~Index() = default;
    explicit Index(T id) noexcept : m_id(id) {}
    
    static constexpr Index Invalid() { return Index(); }
    
    virtual operator bool() const { return m_id != INVALID_ID; }

    explicit operator size_t() const { return static_cast<size_t>(m_id); }

    void makeInvalid() { m_id = INVALID_ID; }

    friend std::hash<Index<Tag,T>>;
    bool operator== (const Index<Tag,T> other) const;
    bool operator!= (const Index<Tag,T> other) const;
    bool operator < (const Index<Tag,T> other) const;
    bool operator > (const Index<Tag,T> other) const;
protected:
    T m_id{INVALID_ID};
};

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator== (const Index<Tag,T> other) const
{
    return m_id == other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator!= (const Index<Tag,T> other) const
{
    return m_id != other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator < (const Index<Tag,T> other) const
{
    return m_id < other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator > (const Index<Tag,T> other) const
{
    return m_id > other.m_id;
}

} //  namespace generic::utils

namespace std {

template<typename Tag, typename T>
struct hash<generic::utils::Index<Tag,T>>
{
    std::size_t operator() (const generic::utils::Index<Tag,T> index) const noexcept
    {
        return std::hash<T>()(index.m_id);
    }
};

template <typename Tag, typename T>
struct less<generic::utils::Index<Tag, T>>
{
    bool operator() (const generic::utils::Index<Tag,T>  lhs, generic::utils::Index<Tag,T> rhs) const noexcept
    {
        return lhs < rhs;
    }
};

template <typename Tag, typename T>
struct equal_to<generic::utils::Index<Tag, T>>
{
    bool operator() (const generic::utils::Index<Tag,T> lhs, generic::utils::Index<Tag,T> rhs) const noexcept
    {
        return lhs == rhs;
    }
};

} // namespace std