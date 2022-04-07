#include <vector>
#include <unordered_set>

#include "types.hpp"
#include "communities_finders.hpp"

Graph* aggregate_communities(const Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();
    
    std::vector<std::vector<vertex_t>> comm_to_nodes(partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        comm_to_nodes[partition.get_comm(v)].push_back(v);
    }
    
    Graph::Builder builder;
    std::vector<float> weights_cache(partition.get_ccount());
    for(community_t c_from = 0; c_from < comm_to_nodes.size(); c_from++){        
        std::unordered_set<community_t> looked_comms;
        for(auto v_from: comm_to_nodes[c_from]) {
            for(auto [v_to, w]: graph_ptr->adjecencies(v_from)) {
                community_t c_to = partition.get_comm(v_to);
                weights_cache[c_to] += w;
                looked_comms.insert(c_to);
            }
        }

        for(auto c_to: looked_comms){
            builder.insert(c_from, c_to, weights_cache[c_to]);
            weights_cache[c_to] = 0;
        }
    }
    return builder.build();
}