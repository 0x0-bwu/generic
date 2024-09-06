/**
 * @file System.hpp
 * @author bwu
 * @brief System related functions
 * @version 0.1
 * @date 2024-08-30
 */
#pragma once

#include "Version.hpp"
#if GENERIC_CURRENT_CXX_VERSION >= 20
    #include <bit>
#else
    #include <boost/endian/arithmetic.hpp>
#endif

namespace generic::common {

inline static constexpr bool isLittleEndian = [] () constexpr {
#if GENERIC_CURRENT_CXX_VERSION >= 20
    return std::endian::native == std::endian::big;
#else
    return boost::endian::order::native == boost::endian::order::little;
#endif
}();

} // namespace generic::common
