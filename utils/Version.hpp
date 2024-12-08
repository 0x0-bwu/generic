/**
 * @file Version.hpp
 * @author bwu
 * @brief version utility
 * @version 0.1
 * @date 2024-12-08
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <cstdint>
namespace generic::utils {

struct Version // format: [Major(99).Minor(99).Patch(999)]
{
    Version() = default;
    explicit constexpr Version(std::uint32_t version)
     : m_major(version / 100000 % 100), m_minor(version % 100000 / 1000), m_patch(version % 1000)
    {
        GENERIC_ASSERT_MSG(version < 10000000, "max version: 99.99.999");
    }

    constexpr Version(int major, int minor, int patch)
     : m_major(major), m_minor(minor), m_patch(patch)
    {
        GENERIC_ASSERT_MSG(major < 100 , "max major version: 99");
        GENERIC_ASSERT_MSG(minor < 100 , "max minor version: 99");
        GENERIC_ASSERT_MSG(patch < 1000, "max patch version: 999");
    }

    explicit operator std::uint32_t() const { return toInt(); }

    bool operator == (const Version & other) const { return toInt() == other.toInt(); }
    bool operator <  (const Version & other) const { return toInt() <  other.toInt(); }
    bool operator >  (const Version & other) const { return toInt() >  other.toInt(); } 

    std::uint32_t toInt() const
    {
        return m_major % 100 * 100000 + m_minor % 100 * 1000 + m_patch % 1000;
    }

    std::string toString() const
    {
        return std::to_string(m_major % 100) + '.' + std::to_string(m_minor % 100) + '.' + std::to_string(m_patch % 1000);
    }
private:
    std::uint8_t m_major{0};
    std::uint8_t m_minor{0};
    std::uint8_t m_patch{0};
};

} // namespace generic::utils

namespace std {

template <>
struct hash<generic::utils::Version>
{
    std::size_t operator() (const generic::utils::Version & version) const noexcept
    {
        return std::hash<std::uint32_t>()(version.toInt());
    }
};

} // namespace std
