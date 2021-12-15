#include <algorithm>    // max
#include "Graph.hpp"
#include "utils.h"

Graph::Builder& Graph::Builder::insert(uint32_t v1, uint32_t v2, float weight) {
    edges.emplace_back(v1, v2, weight);
    return *this;
}

Graph Graph::Builder::build() {
    if(is_built) {
        throw "Graph has already been built! You cannot use build() method more than once!";
    }

    is_built = true;
    if(edges.size() == 0) return Graph();
    
    uint32_t n_vertices = 0;
    std::vector<uint32_t> renumbered_vs;
    for(auto [v1, v2, _]: edges) {
        uint32_t max = std::max(v1, v2) + 1;
        if(renumbered_vs.size() < max) 
            renumbered_vs.resize(max);
        if(renumbered_vs[v1] == 0)
            renumbered_vs[v1] = n_vertices++;
        if(renumbered_vs[v2] == 0)
            renumbered_vs[v2] = n_vertices++;
    }

    Graph graph;    
    graph.total_weight = 0;
    graph.v_to_w.resize(n_vertices);
    graph.v_to_self_w.resize(n_vertices);
    graph.rowptr.resize(n_vertices + 1); // 0 is inserted into beginning for convenience
    graph.col.resize(edges.size());

    for(auto [v1, _, __]: edges) {
        v1 = renumbered_vs[v1];
        graph.rowptr[v1 + 1]++;
    }
    for(uint32_t i = 2; i < graph.rowptr.size(); i++) {
        graph.rowptr[i] += graph.rowptr[i - 1];
    }

    std::vector<uint64_t> indexes = graph.rowptr;
    for(auto [v1, v2, w]: edges) {
        v1 = renumbered_vs[v1];
        v2 = renumbered_vs[v2];
        
        if(v1 == v2) graph.v_to_self_w[v1] = w;
        graph.total_weight += w;
        graph.v_to_w[v1] += w;
        graph.col[indexes[v1]] = std::make_pair(v2, w);
        indexes[v1]++;
    }

    return graph;
}