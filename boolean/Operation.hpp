/**
 * @file Operation.hpp
 * @author bwu
 * @brief boolean related operations
 * @version 0.1
 * @date 2024-09-26
 */
#pragma once

namespace generic::boolean {

template <typename... Args> inline bool All(Args... args) { return (... && args); }

template <typename... Args> inline bool Any(Args... args) { return (... || args); }

} // namespace generic::boolean

