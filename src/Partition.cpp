
#include "types.hpp"
#include "Partition.hpp"

Partition::Partition(std::shared_ptr<const Graph> graph_ptr, std::vector<community_t> vertex_comms) : graph_ptr(graph_ptr) {
    n_comms = 0;
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        community_t comm = vertex_comms[v];
        if(comm > n_comms) n_comms = comm;
    }
    n_comms++;

    comm_sizes.resize(n_comms);
    comm_weights.resize(n_comms);
    comm_inner_weights.resize(n_comms);
    comm_locks = std::vector<std::mutex>(n_comms);

    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        community_t comm = vertex_comms[v];
        comm_sizes[comm]++;
        comm_weights[comm] += graph_ptr->get_weight(v);
        comm_inner_weights[comm] += graph_ptr->get_self_weight(v);
    }
    this->vertex_comms = std::move(vertex_comms);
}

void Partition::change_comm(vertex_t v, community_t c, float old_k_i_in, float new_k_i_in, bool thread_safe) {      
    /* Remove from its old community*/
    community_t old_c = vertex_comms[v];
    if(thread_safe) comm_locks[old_c].lock();
    comm_sizes[old_c]--;
    comm_weights[old_c] -= graph_ptr->get_weight(v);
    comm_inner_weights[old_c] -= 2 * old_k_i_in + graph_ptr->get_self_weight(v);
    if(thread_safe) comm_locks[old_c].unlock();

    /* Insert to its new community */
    if(thread_safe) comm_locks[c].lock();
    comm_sizes[c]++;
    comm_weights[c] += graph_ptr->get_weight(v);
    comm_inner_weights[c] += 2 * new_k_i_in + graph_ptr->get_self_weight(v);
    if(thread_safe) comm_locks[c].unlock();
    
    vertex_comms[v] = c;
}

void Partition::renumber() noexcept {
    std::vector<community_t> renumbered_comms(get_ccount());
    std::vector<community_t> non_empty_comms;
    for(community_t c = 0; c < get_ccount(); c++) {
        if(get_comm_size(c) > 0) {
            renumbered_comms[c] = non_empty_comms.size();
            non_empty_comms.push_back(c);
        } 
    }
    
    std::vector<size_t> new_comm_sizes(non_empty_comms.size());
    std::vector<float> new_comm_weights(non_empty_comms.size());
    std::vector<float> new_comm_inner_weights(non_empty_comms.size());
    for(community_t c = 0; c < non_empty_comms.size(); c++) {
        community_t old_comm = non_empty_comms[c];
        new_comm_sizes[c] = comm_sizes[old_comm];
        new_comm_weights[c] = comm_weights[old_comm];
        new_comm_inner_weights[c] = comm_inner_weights[old_comm];
    }
    
    n_comms = non_empty_comms.size();
    comm_sizes = std::move(new_comm_sizes);
    comm_weights = std::move(new_comm_weights);
    comm_inner_weights = std::move(new_comm_inner_weights);
    comm_locks = std::vector<std::mutex>(n_comms);    
    
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        vertex_comms[v] = renumbered_comms[vertex_comms[v]];
    }
}