#ifndef MODULARITY_HPP
#define MODULARITY_HPP

#include "types.hpp"
#include "Partition.hpp"

class Modularity {
public:
    /* 
        Returns the vertex's best community
        v: The vector
        p: The partition
        k_i_ins: The weights to the all community defined in the partition
    */
    static community_t find_best_community(vertex_t v, const Partition& p, const auto& k_i_ins) {
        auto g = p.get_graph_ptr();

        float _2m = g->get_weight();
        float k_i = g->get_weight(v);
        community_t best_comm = p.get_comm(v);
        
        // We need to isolate the vertex from its current community
        // since we want to regard the vertex as not in any community
        float isolated_sigma_tot = p.get_comm_weight(best_comm) - k_i;
        float max_gain = compute_delta_Q(_2m, isolated_sigma_tot, k_i, k_i_ins.at(best_comm));
        
        for(auto [u, _]: g->adjecencies(v)) {
            community_t c = p.get_comm(u);
            if(p.get_comm(v) == c || !contains(k_i_ins, c)) continue;

            float gain = compute_delta_Q(_2m, p.get_comm_weight(c), k_i, k_i_ins.at(c));
            if(gain > 0 && gain > max_gain) {
                best_comm = c;
                max_gain = gain;
            }
        }
        
        return best_comm;
    }

    static float quality(const Partition& p) noexcept {
        const float _2m = p.get_graph_ptr()->get_weight();
        double q = 0.;
        for (community_t c = 0; c < p.get_ccount(); c++) {
            if (p.get_comm_weight(c) > 0) {
                q += p.get_comm_inner_weight(c) / _2m - (p.get_comm_weight(c) / _2m) * (p.get_comm_weight(c) / _2m);
            }
        }
        return q;
    }

    static float compute_delta_Q(float _2m, float sigma_tot, float k_i, float k_i_in) noexcept {
        return (k_i_in - sigma_tot * k_i / _2m);
    }
};

#endif