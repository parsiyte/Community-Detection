#include "Modularity.hpp"

community_t Modularity::find_best_community(vertex_t v, const Partition& p, const std::vector<float>& k_i_ins) {
    auto g = p.get_graph_ptr();

    float _2m = g->get_weight();
    float k_i = g->get_weight(v);
    community_t best_comm = p.get_comm(v);
    
    // We need to isolate the vertex from its current community
    // since we want to regard the vertex as not in any community
    float isolated_sigma_tot = p.get_comm_weight(best_comm) - k_i;
    float max_gain = compute_delta_Q(_2m, isolated_sigma_tot, k_i, k_i_ins[best_comm]);

    for(auto [u, _]: g->adjecencies(v)) {
        community_t c = p.get_comm(u);
        if(p.get_comm(v) == c) continue;
        
        float gain = compute_delta_Q(_2m, p.get_comm_weight(c), k_i, k_i_ins[c]);
        if(gain > 0 && gain > max_gain) {
            best_comm = c;
            max_gain = gain;
        }
    }
    
    return best_comm;
}