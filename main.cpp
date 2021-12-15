
#include <iostream>
#include <string>

#include "GraphLoader.hpp"
#include "CommunitiesFinder.hpp"
#include <chrono>

using std::chrono::high_resolution_clock;
#define duration_cast(x) (std::chrono::duration_cast<std::chrono::microseconds>(x).count() / 1'000'000.0f)

std::string file_path;

void leiden() {
    std::cout << "LEIDEN:" << std::endl;
    size_t phase = 0;
    Graph curr_graph = load_graph(file_path);
    const float res = 1 / curr_graph.get_weight();
    uint32_t n_vertices_beg;
    do {
        n_vertices_beg = curr_graph.get_vcount();

        LeidenAlgorithm leiden(curr_graph, res);
        std::cout << phase << ") Nodes: " << curr_graph.get_vcount() << ", Edges: " << curr_graph.get_ecount() << ", Quality: " << leiden.quality();
        std::cout.flush();

        auto cp = high_resolution_clock::now();
        bool any_change_occured = leiden.move_vertices();
        std::cout << ", Moving Nodes: " << duration_cast(high_resolution_clock::now() - cp);
        std::cout.flush();

        if(!any_change_occured) break;
        
        // cp = high_resolution_clock::now();
        // leiden.refine_partition();
        // std::cout << ", Refine: " << duration_cast(high_resolution_clock::now() - cp);
        // std::cout.flush();

        cp = high_resolution_clock::now();
        curr_graph = leiden.aggregate_communities();
        std::cout << ", Aggregate Communities: " << duration_cast(high_resolution_clock::now() - cp) << std::endl;
        phase++;
    } while(curr_graph.get_vcount() < n_vertices_beg);
    std::cout << std::endl;
}

void louvain() {
    std::cout << "LOUVAIN:" << std::endl;
    size_t phase = 0;
    Graph curr_graph = load_graph(file_path);
    while(true) {
        LouvainAlgorithm louvain(curr_graph);
        std::cout << phase << ") Nodes: " << curr_graph.get_vcount() << ", Edges: " << curr_graph.get_ecount() << ", Quality: " << louvain.quality();
        std::cout.flush();

        auto cp = high_resolution_clock::now();
        bool any_change_occured = louvain.move_vertices();
        std::cout << ", Moving Nodes: " << duration_cast(high_resolution_clock::now() - cp);
        std::cout.flush();

        if(!any_change_occured) break;

        cp = high_resolution_clock::now();
        curr_graph = louvain.aggregate_communities();
        std::cout << ", Aggregate Communities: " << duration_cast(high_resolution_clock::now() - cp) << std::endl;
        phase++;
    }
    std::cout << std::endl;
}

int main(int argc, const char* argv[]) {
    file_path = std::string(argv[2]);
    if(std::string(argv[1]) == "louvain")
        louvain();
    else if(std::string(argv[1]) == "leiden")
        leiden();
}