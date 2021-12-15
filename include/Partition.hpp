#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <cstdint>          // uint32_t, uint64_t
#include <vector>
#include "Graph.hpp"

class Partition {
    const Graph* graph;
    std::vector<uint32_t> v_to_comm;
    std::vector<float> comm_weights; // the sum of the weights of the edges incident to vertices in the community
    std::vector<float> comm_inner_weights; // the sum of the weights of the edges inside the community
public:
    Partition(const Graph& graph);

    void change_comm(uint32_t v, uint32_t c, float old_k_i_in, float new_k_i_in);
    
    uint32_t get_ccount() const noexcept { return v_to_comm.size(); }

    uint32_t get_comm(uint32_t v) const { return v_to_comm.at(v); }
    
    // Returns the sum of the weights of the edges incident to vertices in the community "c"
    float get_comm_weight(uint32_t c) const { return comm_weights.at(c); }
    
    // Returns the sum of the weights of the edges inside the community "c"
    float get_comm_inner_weight(uint32_t c) const { return comm_inner_weights.at(c); }

    const Graph& get_graph() const noexcept { return *graph; }
};

#endif