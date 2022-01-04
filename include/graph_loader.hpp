#ifndef GRAPH_LOADER_HPP
#define GRAPH_LOADER_HPP

#include <vector>
#include <string>

#include "types.hpp"

std::vector<edge_t> load_edges(std::string path);
std::vector<vertex_t> get_vertices(const std::vector<edge_t>& edges);
#endif