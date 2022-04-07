#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <algorithm>    // max

#include "types.hpp"
#include "graph_loader.hpp"

std::vector<edge_t> load_edges(std::string path) {
    std::ifstream finput;
    finput.open(path, std::fstream::in);
    if(!finput.is_open() || finput.fail()) {
        std::cerr << "Error while opening " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<edge_t> edges;
    std::string line;
    while(std::getline(finput, line)) {
        int idx = -1;
        while(line[++idx] == ' ');
        if(line[idx] == '#') continue;
        
        vertex_t v_from, v_to;
        std::istringstream(line) >> v_from >> v_to;
        edges.emplace_back(v_from, v_to, 1.0f);
    }
    finput.close();

    return edges;
}

std::vector<vertex_t> get_vertices(const std::vector<edge_t>& edges) {
    vertex_t max_vertex = 0;
    for(auto [v_from, v_to, _]: edges) {
        max_vertex = std::max(max_vertex, std::max(v_from, v_to));
    }
    size_t n_vertices = max_vertex + 1;
    
    std::vector<bool> is_vertex_added(n_vertices);
    std::vector<vertex_t> vertices;
    for(auto [v_from, v_to, _]: edges) {
        if(!is_vertex_added[v_from]) {
            vertices.push_back(v_from);
            is_vertex_added[v_from] = true;
        }
        if(!is_vertex_added[v_to]) {
            vertices.push_back(v_to);
            is_vertex_added[v_to] = true;
        }
    }
    return vertices;
}