#pragma once

#include "generic/common/Exception.hpp"
#include "Interpolation.hpp"
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

    constexpr Values & operator* () { return m_values; }

    constexpr const Values & operator* () const { return m_values; }
    
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
    
    Scalar Lookup(Scalar x) const
    {
        static_assert(DIM == 1, "only 1D lookup is supported");
        if (not isValid() or m_values.empty()) {
            GENERIC_ASSERT_MSG(false, "invalid lookup table");
            return 0;
        }
        if (1 == m_values.size()) return m_values.front();
        Scalar x1, x2;
        auto i = InterpIndex(0, x, x1, x2);
        auto q1 = this->operator()(i), q2 = this->operator()(i + 1);
        return LinearInterpolation(q1, q2, x1, x2, x);
    }

    Scalar Lookup(Scalar x, Scalar y) const
    {
        static_assert(DIM == 2, "only 2D lookup is supported");
        if (not isValid() or m_values.empty()) {
            GENERIC_ASSERT_MSG(false, "invalid lookup table");
            return 0;
        }
        if (1 == m_values.size()) return m_values.front();
        Scalar x1, x2, y1, y2;
        auto i1 = InterpIndex(0, x, x1, x2);
        auto i2 = InterpIndex(1, y, y1, y2);
        auto q11 = this->operator()(i1, i2);
        auto q21 = this->operator()(i1 + 1, i2);
        auto q12 = this->operator()(i1, i2 + 1);
        auto q22 = this->operator()(i1 + 1, i2 + 1);
        return BilinearInterpolation(q11, q12, q21, q22, x1, x2, y1, y2, x, y);
    }
    
private:
    size_t InterpIndex(size_t dim, Scalar x, Scalar & x1, Scalar & x2) const
    {
        if (std::isnan(x)) { GENERIC_ASSERT(false); return 0; }
        if (x < m_indices[dim].front()) {
            x1 = m_indices[dim][0];
            x2 = m_indices[dim][1];
            return 0;
        }
        else if (x > m_indices[dim].back()) {
            x1 = m_indices[dim][m_indices[dim].size() - 2];
            x2 = m_indices[dim][m_indices[dim].size() - 1];
            return m_indices[dim].size() - 2;
        }
        else {
            size_t l{0}, r{m_indices[dim].size() - 1}, m;
            while (r > l + 1) {
                m = (l + r) / 2;
                if (m_indices[dim][m] > x) r = m;
                else l = m;
            }
            x1 = m_indices[dim][l];
            x2 = m_indices[dim][r];
            return l;
        }
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