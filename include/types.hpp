#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstddef>  // size_t
#include <tuple>

using size_t = std::size_t;
using vertex_t = size_t;
using community_t = size_t;
using edge_t = std::tuple<vertex_t, vertex_t, float>;

#endif