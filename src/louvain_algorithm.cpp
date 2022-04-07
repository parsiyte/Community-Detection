#include <iostream>
#include <vector>

#include "types.hpp"
#include "communities_finders.hpp"
#include "Modularity.hpp"

namespace {
bool move_vertices(Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();

    std::vector<float> k_i_ins(partition.get_ccount());
    bool any_changes_occured = false;
    bool any_changes_occured_in_last_pass;
    do {
        any_changes_occured_in_last_pass = false;
        for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
            for(auto [u, w]: graph_ptr->adjecencies(v)) {
                if(v == u) continue;
                k_i_ins[partition.get_comm(u)] += w;
            }

            community_t initial_community = partition.get_comm(v);
            community_t best_community = Modularity::find_best_community(v, partition, k_i_ins);
            
            if(initial_community != best_community) {
                partition.change_comm(v, best_community, k_i_ins[initial_community], k_i_ins[best_community]);
                any_changes_occured = any_changes_occured_in_last_pass = true;
            }

            for(auto [u, _]: graph_ptr->adjecencies(v)) {
                k_i_ins[partition.get_comm(u)] = 0;
            }
        }
    } while(any_changes_occured_in_last_pass);
    return any_changes_occured;
}

void run(std::vector<Partition>& hierarchy, std::shared_ptr<const Graph> graph_ptr, int max_level) {
    std::cout << "Nodes: " << graph_ptr->get_vcount() << ", Edges: " << graph_ptr->get_ecount(); std::cout.flush();
    if (max_level == 0) return;
    
    hierarchy.emplace_back(graph_ptr);
    auto& partition = hierarchy.back();
    std::cout << ", Quality: " << Modularity::quality(partition); std::cout.flush();

    auto cp = get_time();
    bool any_change_occured = move_vertices(partition); partition.renumber();
    std::cout << ", Moving nodes: " << get_elapsed_time(cp, get_time()); std::cout.flush();
    if(!any_change_occured) return;

    cp = get_time();
    std::shared_ptr<const Graph> new_graph(aggregate_communities(partition));
    std::cout << ", Aggregate communities: " << get_elapsed_time(cp, get_time()) << std::endl;

    run(hierarchy, new_graph, max_level - 1);
}
}

std::vector<Partition> run_louvain(std::shared_ptr<const Graph> graph_ptr, int max_level) {
    std::vector<Partition> hierarchy;
    run(hierarchy, graph_ptr, max_level);
    return hierarchy;
}