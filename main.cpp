
#include <iostream>
#include <fstream>
#include <algorithm>    // max_element
#include <string>
#include <memory>

#include "types.hpp"
#include "graph_loader.hpp"
#include "communities_finders.hpp"

int main(int argc, const char* argv[]) {
    std::vector<edge_t> edges = load_edges(std::string(argv[2]));
    std::vector<vertex_t> vertices = get_vertices(edges);
    size_t n_vertices = *std::max_element(vertices.begin(), vertices.end()) + 1;
    std::vector<vertex_t> renumbered_vertices(n_vertices);
    vertex_t vertex = 0;
    for(auto v: vertices) renumbered_vertices[v] = vertex++;

    Graph::Builder builder;
    for(auto [v1, v2, _]: edges) {
        v1 = renumbered_vertices[v1]; v2 = renumbered_vertices[v2];
        builder.insert(v1, v2); builder.insert(v2, v1);
    }
    std::shared_ptr<const Graph> graph_ptr(builder.build());
    edges.clear();

    std::vector<Partition> partitions;
    if(std::string(argv[1]) == "louvain") {
        partitions = run_louvain(graph_ptr, -1);
    }
    else if(std::string(argv[1]) == "leiden") {
        partitions = run_leiden(graph_ptr, -1);
    }
    std::cout << std::endl;
/*
    auto& partition = partitions.front();
    std::vector<std::vector<vertex_t>> comm_to_nodes(partition.get_ccount());
    for(auto v: vertices) {
        comm_to_nodes[partition.get_comm(renumbered_vertices[v])].push_back(v);
    }
    std::ofstream fout;
    fout.open("comms.txt", std::fstream::out);
    for(auto& vertices: comm_to_nodes){
        if(vertices.size() == 0) continue;
        for(auto v: vertices) fout << v << " ";
        fout << std::endl;
    }
    fout.close();
*/
}