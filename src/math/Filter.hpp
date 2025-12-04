/**
 * @file Filter.hpp
 * @author bwu
 * @brief filter related algorithms
 * @version 0.1
 * @date 2024-05-13
 */
#pragma once
#include "MathUtility.hpp"
#include <vector>

namespace generic::math {

template <typename T, bool UseUpdatedValue = true, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
std::vector<T> SimpleMovingAverage(const T * data, size_t size, size_t length)
{
    std::vector<T> results(size);
    for (int i = 0; i < int(size); ++i) {
        T ave = 0;
        for (int j = 0; j < int(length); ++j) {
            auto iL = std::max<int>(0, i - j - 1);
            auto iR = std::min<int>(size - 1, i + j + 1);
            if constexpr (UseUpdatedValue) {
                ave += (iL < i ? results[iL] : data[iL]) + data[iR];
            }
            else ave += data[iL] + data[iR];
        }
        results[i] = (ave + data[i]) / (2 * length + 1);        
    }
    return results;
}

} // namespace generic::math
