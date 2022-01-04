#include <vector>
#include <limits>
#include <algorithm>    // max

#include "types.hpp"
#include "Graph.hpp"

Graph::Builder& Graph::Builder::insert(vertex_t v1, vertex_t v2, float weight) {
    edges.emplace_back(v1, v2, weight);
    return *this;
}

Graph* Graph::Builder::build() {
    if(is_built) {
        throw "Graph has already been built! You cannot use build() method more than once!";
    }

    is_built = true;
    if(edges.size() == 0) return new Graph();
    
    vertex_t max_vertex = 0;
    for(auto [v_from, v_to, _]: edges) {
        max_vertex = std::max(max_vertex, std::max(v_from, v_to));
    }
    size_t n_vertices = max_vertex + 1;

    Graph* graph = new Graph();    
    graph->total_weight = 0;
    graph->vertex_weights.resize(n_vertices);
    graph->vertex_self_weights.resize(n_vertices);
    graph->rowptr.resize(n_vertices + 1); // 0 is inserted into beginning for convenience
    graph->col.resize(edges.size());

    for(auto [v, _, __]: edges) {
        graph->rowptr[v + 1]++;
    }
    for(vertex_t v = 2; v < graph->rowptr.size(); v++) {
        graph->rowptr[v] += graph->rowptr[v - 1];
    }

    std::vector<size_t> indexes = graph->rowptr;
    for(auto [v_from, v_to, w]: edges) {
        if(v_from == v_to) graph->vertex_self_weights[v_from] = w;
        graph->total_weight += w;
        graph->vertex_weights[v_from] += w;
        graph->col[indexes[v_from]] = std::make_pair(v_to, w);
        indexes[v_from]++;
    }
    edges.clear();
    
    return graph;
}