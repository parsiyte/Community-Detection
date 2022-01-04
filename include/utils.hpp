#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cstdlib>      // abs
#include <chrono>
#include <numeric>      // iota

#include "types.hpp"

inline bool is_equal(float value1, float value2, float eps=1e-10) {
    return std::abs(value1 - value2) < eps;
}

inline std::vector<community_t> fill_with_communities(size_t n_comms) {
    std::vector<community_t> arr(n_comms);
    std::iota(arr.begin(), arr.end(), 0);
    return arr;
}

inline std::chrono::system_clock::time_point get_time() { return  std::chrono::high_resolution_clock::now(); }
inline double get_elapsed_time(std::chrono::system_clock::time_point last_check_point, std::chrono::system_clock::time_point curr_check_point) {
    return std::chrono::duration_cast<std::chrono::microseconds>(curr_check_point - last_check_point).count() / 1'000'000.0f;
}

#endif