#pragma once

#include "generic/common/Exception.hpp"
#include <concepts>
#include <numeric>
#include <vector>
#include <array>

namespace generic::math {

template <typename Scalar, std::size_t DIM>
class LookupTable
{
#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & boost::serialization::make_nvp("indices", m_indices);
        ar & boost::serialization::make_nvp("values" , m_values );
    }
#endif
public:
    using Values = std::vector<Scalar>;
    using Indices = std::array<Values, DIM>;
    LookupTable() = default;
    LookupTable(Indices indices, Values values)
     : m_indices(std::move(indices)), m_values(std::move(values))
    {
        GENERIC_ASSERT_MSG(isValid(), "value size mismatch with indices");
    }
    
    constexpr size_t Dim() const { return DIM; }

    constexpr Values & operator[] (size_t index)
    {
        return m_indices[index];
    }

    constexpr const Values & operator[] (size_t index) const
    {
        return m_indices.at(index);
    }
    
    template <typename... Args>
    Scalar & operator() (Args... indices)
    {
        static_assert(sizeof...(indices) == DIM, "number of indices must match the number of dimensions");
        return m_values[CalcIndex({static_cast<std::size_t>(indices)...})];
    }

    template <typename... Args>
    const Scalar & operator() (Args... indices) const
    {
        static_assert(sizeof...(indices) == DIM, "number of indices must match the number of dimensions");
        return m_values.at(CalcIndex({static_cast<std::size_t>(indices)...}));
    }

    bool isValid() const
    {
        return CalcSize(m_indices) == m_values.size();
    }

private:
    Indices m_indices;
    Values m_values;

    static std::size_t CalcSize(const Indices & indices)
    {
        return std::accumulate(indices.begin(), indices.end(), 1, [](std::size_t acc, const Values & v) { return acc * v.size(); });
    }

    std::size_t CalcIndex(const std::array<std::size_t, DIM> & indices) const
    {
        std::size_t index = 0;
        std::size_t multiplier = m_values.size();
        for (size_t i = 0; i < DIM; ++i) {
            if (indices[i] >= m_indices[i].size()) {
                GENERIC_ASSERT_MSG(false, "index out of range");
            }
            multiplier /= m_indices[i].size();
            index += indices[i] * multiplier;
        }
        return index;
    }
};

} // namespace generic::math