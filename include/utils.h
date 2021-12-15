#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <chrono>

inline bool is_equal(float value1, float value2, float eps=1e-10) {
    return std::abs(value1 - value2) < eps;
}
inline bool is_not_equal(float value1, float value2, float eps=1e-10) {
    return !is_equal(value1, value2, eps);
}
inline float take_square(float value) { return value * value; }

using time_point = std::chrono::system_clock::time_point;
inline time_point get_time() { return  std::chrono::high_resolution_clock::now(); }
inline double get_elapsed_time(time_point last_check_point, time_point curr_check_point) {
    return std::chrono::duration_cast<std::chrono::microseconds>(curr_check_point - last_check_point).count() / 1'000'000.0f;
}

#endif