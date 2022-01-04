#ifndef MODULARITY_HPP
#define MODULARITY_HPP

#include <vector>

#include "types.hpp"
#include "Partition.hpp"

class Modularity {
public:
    /* 
        Returns the vertex's best community
        v: The vector
        p: The partition
        k_i_ins: The weights from the vector to the all community defined in the partition
    */    
    static vertex_t find_best_community(vertex_t v, const Partition& p, const std::vector<float>& k_i_ins);

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