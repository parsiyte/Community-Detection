#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <unordered_map>
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

using time_point = std::chrono::system_clock::time_point;
inline time_point get_time() { 
    return  std::chrono::high_resolution_clock::now(); 
}
inline double get_elapsed_time(time_point last_check_point,  time_point curr_check_point) {
    return std::chrono::duration_cast<std::chrono::microseconds>(curr_check_point - last_check_point).count() / 1'000'000.0f;
}

template<typename K, typename V>
bool contains(const std::vector<V>& vec, K key) {
    return key < vec.size();
}

template<typename K, typename V>
bool contains(const std::unordered_map<K, V>& map, K key) {
    return map.find(key) != map.end();
}
#endif