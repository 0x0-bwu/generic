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
/**
 * @brief Brief description of Index.
 * @param m_id(id
 * @return explicit
 */
    explicit Index(T id) noexcept : m_id(id) {}
        
    virtual operator bool() const { return m_id != INVALID_ID; }

    explicit operator size_t() const { return static_cast<size_t>(m_id); }

/**
 * @brief Brief description of makeInvalid.
 * @return void
 */
    void makeInvalid() { m_id = INVALID_ID; }

    friend std::hash<Index<Tag,T>>;
    bool operator== (const Index<Tag,T> & other) const;
    bool operator!= (const Index<Tag,T> & other) const;
    bool operator < (const Index<Tag,T> & other) const;
    bool operator > (const Index<Tag,T> & other) const;

#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int);
#endif//GENERIC_BOOST_SERIALIZATION_SUPPORT
protected:
    T m_id{INVALID_ID};
};

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator== (const Index<Tag,T> & other) const
{
    return m_id == other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator!= (const Index<Tag,T> & other) const
{
    return m_id != other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator < (const Index<Tag,T> & other) const
{
    return m_id < other.m_id;
}

template<typename Tag, typename T>
inline bool Index<Tag,T>::operator > (const Index<Tag,T> & other) const
{
    return m_id > other.m_id;
}

#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT
template<typename Tag, typename T>
template <typename Archive>
inline void Index<Tag,T>::serialize(Archive & ar, const unsigned int)
{
    ar & boost::serialization::make_nvp("id", m_id); 
}
#endif//GENERIC_BOOST_SERIALIZATION_SUPPORT

} //  namespace generic::utils

namespace std {

template<typename Tag, typename T>
struct hash<generic::utils::Index<Tag,T>>
{
    std::size_t operator() (const generic::utils::Index<Tag,T> & index) const noexcept
    {
        return std::hash<T>()(index.m_id);
    }
};

template <typename Tag, typename T>
struct less<generic::utils::Index<Tag, T>>
{
    bool operator() (const generic::utils::Index<Tag,T> & lhs, const generic::utils::Index<Tag,T> & rhs) const noexcept
    {
        return lhs < rhs;
    }
};

template <typename Tag, typename T>
struct equal_to<generic::utils::Index<Tag, T>>
{
    bool operator() (const generic::utils::Index<Tag,T> & lhs, const generic::utils::Index<Tag,T> & rhs) const noexcept
    {
        return lhs == rhs;
    }
};

} // namespace std
