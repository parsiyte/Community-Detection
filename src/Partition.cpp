#include "Partition.hpp"

Partition::Partition(const Graph& graph) : graph(&graph) {
    v_to_comm.resize(graph.get_vcount());
    comm_weights.resize(graph.get_vcount());
    comm_inner_weights.resize(graph.get_vcount());

    for(uint32_t v = 0; v < graph.get_vcount(); v++) {
        v_to_comm[v] = v;
        comm_weights[v] = graph.get_weight(v);
        comm_inner_weights[v] = graph.get_self_weight(v);
    }
}

void Partition::change_comm(uint32_t v, uint32_t c, float old_k_i_in, float new_k_i_in) {      
    /* Remove from its old community*/
    uint32_t old_c = v_to_comm[v];
    comm_weights[old_c] -= graph->get_weight(v);
    comm_inner_weights[old_c] -= 2 * old_k_i_in + graph->get_self_weight(v);
    
    /* Insert to its new community */
    comm_weights[c] += graph->get_weight(v);
    comm_inner_weights[c] += 2 * new_k_i_in + graph->get_self_weight(v);

    v_to_comm[v] = c;
}