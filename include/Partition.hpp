#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <vector>
#include <memory>           // shared_ptr
#include <mutex>

#include "utils.hpp"
#include "Graph.hpp"

class Partition {
    size_t n_comms;
    std::shared_ptr<const Graph> graph_ptr;
    std::vector<community_t> vertex_comms;
    std::vector<size_t> comm_sizes;
    std::vector<std::mutex> comm_locks;
    std::vector<float> comm_weights; // the sum of the weights of the edges incident to vertices in the community
    std::vector<float> comm_inner_weights; // the sum of the weights of the edges inside the community
public:
    Partition(std::shared_ptr<const Graph> graph_ptr, std::vector<community_t> vertex_comms);
    Partition(std::shared_ptr<const Graph> graph_ptr) : Partition(graph_ptr, fill_with_communities(graph_ptr->get_vcount())) {};
    
    void change_comm(vertex_t v, community_t c, float old_k_i_in, float new_k_i_in, bool thread_safe=false);
    void renumber() noexcept;

    size_t get_ccount() const noexcept { return n_comms; }
    community_t get_comm(vertex_t v) const { return vertex_comms.at(v); }
    size_t get_comm_size(community_t c) const { return comm_sizes.at(c); }

    // Returns the sum of the weights of the edges incident to vertices in the community "c"
    float get_comm_weight(community_t c) const { return comm_weights.at(c); }
    
    // Returns the sum of the weights of the edges inside the community "c"
    float get_comm_inner_weight(community_t c) const { return comm_inner_weights.at(c); }

    std::shared_ptr<const Graph> get_graph_ptr() const noexcept { return graph_ptr; }

    const std::vector<community_t>& get_vertex_comms() const noexcept { return vertex_comms; }
};

#endif