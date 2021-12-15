#include "CommunitiesFinder.hpp"
#include <cstddef>              //size_t
#include <algorithm>            //iota
#include <cmath>                //exp
#include <deque>
#include <random>
#include <unordered_set>
#include "utils.h"

Graph CommunitiesFinder::aggregate_communities() const noexcept {
    auto& graph = partition.get_graph();
    
    std::vector<std::vector<uint32_t>> comm_to_nodes(graph.get_vcount());
    for(uint32_t v = 0; v < graph.get_vcount(); v++) {
        comm_to_nodes[partition.get_comm(v)].push_back(v);
    }
    
    Graph::Builder builder;
    std::vector<float> weights_cache(graph.get_vcount());
    for(uint32_t c_from = 0; c_from < comm_to_nodes.size(); c_from++){        
        if(comm_to_nodes[c_from].size() == 0) continue;
        
        std::unordered_set<uint32_t> looked_comms;
        for(auto v_from: comm_to_nodes[c_from]) {
            for(auto [v_to, w]: graph.adjecencies(v_from)) {
                uint32_t c_to = partition.get_comm(v_to);
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

bool LouvainAlgorithm::move_vertices() noexcept {
    auto& graph = partition.get_graph();

    std::vector<float> k_i_ins(partition.get_ccount());
    bool any_changes_occured = false;
    bool any_changes_occured_in_last_pass;
    do {
        any_changes_occured_in_last_pass = false;
        for(uint32_t v = 0; v < graph.get_vcount(); v++) {
            for(auto [u, w]: graph.adjecencies(v)) {
                if(v == u) continue;
                k_i_ins[partition.get_comm(u)] += w;
            }

            uint32_t initial_community = partition.get_comm(v);
            uint32_t best_community = Modularity::find_best_community(v ,partition, k_i_ins);
            
            if(initial_community != best_community) {
                partition.change_comm(v, best_community, k_i_ins[initial_community], k_i_ins[best_community]);
                any_changes_occured = any_changes_occured_in_last_pass = true;
            }

            for(auto [u, _]: graph.adjecencies(v)) {
                if(v == u) continue;
                k_i_ins[partition.get_comm(u)] = 0;
            }
        }
    } while(any_changes_occured_in_last_pass);
    return any_changes_occured;
}

bool LeidenAlgorithm::move_vertices() noexcept {
    auto& graph = partition.get_graph();

    std::deque<uint32_t> v_queue(graph.get_vcount());
    std::iota(v_queue.begin(), v_queue.end(), 0);
    std::random_device rd;
    std::default_random_engine g(rd());
    std::shuffle(v_queue.begin(), v_queue.end(), g);
    std::vector<bool> v_in_queue(graph.get_vcount(), true);

    std::vector<float> k_i_ins(partition.get_ccount());
    bool any_changes_occured = false;
    while(!v_queue.empty()) {
        uint32_t v = v_queue.front(); v_queue.pop_front(); v_in_queue[v] = false;
        for(auto [u, w]: graph.adjecencies(v)) {
            if(v == u) continue;
            k_i_ins[partition.get_comm(u)] += w;
        }

        uint32_t initial_community = partition.get_comm(v);
        uint32_t best_community = Modularity::find_best_community(v ,partition, k_i_ins);
        
        if (initial_community != best_community) {
            partition.change_comm(v, best_community, k_i_ins[initial_community], k_i_ins[best_community]);
            any_changes_occured = true;
            
            for(auto [u, _]: graph.adjecencies(v)) {
                if(!v_in_queue[u] && partition.get_comm(u) != best_community){
                    v_in_queue[u] = true;
                    v_queue.push_back(u);
                } 
            }
        }

        for(auto [u, _]: graph.adjecencies(v)) {
            if(v == u) continue;
            k_i_ins[partition.get_comm(u)] = 0;
        }
    }
    return any_changes_occured;   
}

// inline float uni_rand() {
//     static std::random_device rd;
//     static std::default_random_engine g(rd());
//     static std::uniform_real_distribution<float> unif(0, 1);
//     return unif(g);
// }
    
// inline uint32_t find_index(const std::vector<double>& vec, double searched_val) {
//     int min_idx = 0, mid_idx, max_idx = vec.size() - 1;
//     while (min_idx <= max_idx){
//         mid_idx = (min_idx + max_idx) / 2;
//         double prev_val = mid_idx == 0 ? 0 : vec.at(mid_idx - 1);
//         double curr_val = vec.at(mid_idx);
//         if(prev_val <= searched_val && searched_val <= curr_val) break;
//         else if (prev_val <= searched_val)
//             min_idx = mid_idx + 1;
//         else
//             max_idx = mid_idx;
//     }
//     return mid_idx;
// }

// void LeidenAlgorithm::refine_partition() noexcept {
//     auto& graph = partition.get_graph();

//     Partition ref_p(graph);
//     std::vector<float> c_to_s(graph.get_vcount());
//     for(uint32_t v = 0; v < graph.get_vcount(); v++) {
//         for(auto [u, w]: graph.adjecencies(v)) {
//             if(v == u || partition.get_comm(v) != partition.get_comm(u)) continue;
//             c_to_s[v] += w;
//         }
//     }
//     // auto& v_to_s = c_to_s;

//     std::vector<uint64_t> comm_sizes(ref_p.get_ccount(), 1);
//     std::vector<float> k_i_ins(ref_p.get_ccount());
//     for(uint32_t v = 0; v < graph.get_vcount(); v++) {
//         if(comm_sizes[ref_p.get_comm(v)] != 1) continue;
//         // if(v_to_s[v] < graph.get_weight(v) * (partition.get_comm_weight(partition.get_comm(v)) - graph.get_weight(v)) * res) continue;
        
//         for(auto [u, w]: graph.adjecencies(v)) {
//             if(v == u) continue;
//             k_i_ins[ref_p.get_comm(u)] += w;
//         }
        
//         double cumm_exp, max_gain = 0;
//         uint32_t best_comm_idx;
//         std::vector<uint32_t> comms;
//         std::vector<double> qualities;
//         std::unordered_set<uint32_t> visited_comms;
//         for(auto [u, _]: graph.adjecencies(v)) {
//             if(partition.get_comm(v) != partition.get_comm(u)) continue;
            
//             uint32_t c = ref_p.get_comm(u);
//             if(visited_comms.find(c) != visited_comms.end()) continue;
//             visited_comms.insert(c);

//             if(c_to_s[c] >= ref_p.get_comm_weight(c) * (partition.get_comm_weight(partition.get_comm(v)) - ref_p.get_comm_weight(c)) * res) {
//                 // No need to isolate the vertex from its community when "c" equals to its community
//                 // since this function runs for only the vertex "v" whose community size is 1
//                 float gain = Modularity::compute_delta_Q(graph.get_weight(), ref_p.get_comm_weight(c), graph.get_weight(v), k_i_ins[c]);
//                 if(gain > 0) {
//                     if(gain > max_gain) {
//                         max_gain = gain;
//                         best_comm_idx = comms.size();
//                     }

//                     cumm_exp += std::exp(gain / 0.01); 
//                     comms.push_back(c); 
//                     qualities.push_back(cumm_exp);
//                 }
//             }
//         }

//         if(!comms.empty()) {
//             if(cumm_exp < std::numeric_limits<double>::infinity()) {
//                 best_comm_idx = find_index(qualities, cumm_exp * uni_rand());
//             }
//             uint32_t best_comm = comms.at(best_comm_idx);
            
//             c_to_s[ref_p.get_comm(v)] = 0;
//             for(auto [u, w]: graph.adjecencies(v)) {
//                 if(v == u || partition.get_comm(v) != partition.get_comm(u)) continue;
//                 c_to_s[best_comm] += best_comm == partition.get_comm(u) ? -w : w;
//             }

//             comm_sizes[ref_p.get_comm(v)]--;
//             comm_sizes[best_comm]++;
//             ref_p.change_comm(v, best_comm, k_i_ins[ref_p.get_comm(v)], k_i_ins[best_comm]);
//         }

//         for(auto [u, _]: graph.adjecencies(v)) {
//             if(v == u) continue;
//             k_i_ins[ref_p.get_comm(u)] = 0;
//         }
//     }

//     partition = std::move(ref_p);
// }