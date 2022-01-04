#include <iostream>
#include <vector>
#include <memory>
#include <deque>
#include <random>               // shuffle
#include <algorithm>            // iota

#include "types.hpp"
#include "communities_finders.hpp"
#include "Modularity.hpp"
namespace {
std::random_device RANDOM_DEVICE;
std::default_random_engine RANDOM_ENGINE(RANDOM_DEVICE());
std::uniform_real_distribution<float> UNIFORM_DIST(0, 1); 

inline size_t find_index(const std::vector<double>& vec, double searched_val) {
    size_t min_idx = 0, mid_idx, max_idx = vec.size() - 1;
    while (min_idx <= max_idx){
        mid_idx = (min_idx + max_idx) / 2;
        double prev_val = mid_idx == 0 ? 0 : vec.at(mid_idx - 1);
        double curr_val = vec.at(mid_idx);
        if(prev_val <= searched_val && searched_val <= curr_val) break;
        else if (prev_val <= searched_val)
            min_idx = mid_idx + 1;
        else
            max_idx = mid_idx;
    }
    return mid_idx;
}

bool move_vertices(Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();

    std::deque<vertex_t> v_queue(graph_ptr->get_vcount());
    std::iota(v_queue.begin(), v_queue.end(), 0);
    std::shuffle(v_queue.begin(), v_queue.end(), RANDOM_ENGINE);
    std::vector<bool> v_in_queue(graph_ptr->get_vcount(), true);

    std::vector<float> k_i_ins(partition.get_ccount());
    bool any_changes_occured = false;
    while(!v_queue.empty()) {
        vertex_t v = v_queue.front(); v_queue.pop_front(); v_in_queue[v] = false;
        for(auto [u, w]: graph_ptr->adjecencies(v)) {
            if(v == u) continue;
            k_i_ins[partition.get_comm(u)] += w;
        }

        community_t initial_community = partition.get_comm(v);
        community_t best_community = Modularity::find_best_community(v ,partition, k_i_ins);
        
        if (initial_community != best_community) {
            partition.change_comm(v, best_community, k_i_ins[initial_community], k_i_ins[best_community]);
            any_changes_occured = true;
            
            for(auto [u, _]: graph_ptr->adjecencies(v)) {
                if(!v_in_queue[u] && partition.get_comm(u) != best_community){
                    v_in_queue[u] = true;
                    v_queue.push_back(u);
                } 
            }
        }

        for(auto [u, _]: graph_ptr->adjecencies(v)) {
            k_i_ins[partition.get_comm(u)] = 0;
        }
    }
    return any_changes_occured;   
}

void refine_partition(Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();

    Partition refined_partition(graph_ptr);
    std::vector<float> ref_comm_to_comm_w(graph_ptr->get_vcount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        for(auto [u, w]: graph_ptr->adjecencies(v)) {
            if(partition.get_comm(v) != partition.get_comm(u)) {
                ref_comm_to_comm_w[refined_partition.get_comm(u)] += w;
            }
        }
    }

    std::vector<float> k_i_ins(refined_partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        community_t v_ref_comm = refined_partition.get_comm(v);
        if(refined_partition.get_comm_size(v_ref_comm) > 1) continue;

        float comm_weight = partition.get_comm_weight(partition.get_comm(v));
        if(ref_comm_to_comm_w[v_ref_comm] < refined_partition.get_comm_weight(v_ref_comm) * (comm_weight - refined_partition.get_comm_weight(v_ref_comm)) / graph_ptr->get_weight()) continue;
        ref_comm_to_comm_w[v_ref_comm] = 0;

        std::vector<community_t> target_ref_comms;
        target_ref_comms.push_back(v_ref_comm);
        for(auto [u, w]: graph_ptr->adjecencies(v)) {
            if(partition.get_comm(v) != partition.get_comm(u)) continue;
            community_t u_ref_comm = refined_partition.get_comm(u);
            if(ref_comm_to_comm_w[u_ref_comm] < refined_partition.get_comm_weight(u_ref_comm) * (comm_weight - refined_partition.get_comm_weight(u_ref_comm)) / graph_ptr->get_weight()) continue;
            if(is_equal(k_i_ins[u_ref_comm], 0)) {
                target_ref_comms.push_back(u_ref_comm);
            } 
            k_i_ins[u_ref_comm] += w;
        }

        community_t best_comm = v_ref_comm;
        double max_gain = 0, cumm_exp = 0;
        std::vector<double> qualities;
        for(auto ref_comm: target_ref_comms){
            float gain = Modularity::compute_delta_Q(graph_ptr->get_weight(), refined_partition.get_comm_weight(ref_comm), graph_ptr->get_weight(v), k_i_ins[ref_comm]);
            if(gain >= 0) {
                if(gain > max_gain) {
                    max_gain = gain;
                    best_comm = ref_comm;
                }
                cumm_exp += std::exp(gain / 0.01); 
            }
            qualities.push_back(cumm_exp);
        }

        if(qualities.back() < std::numeric_limits<double>::infinity()) {
            best_comm = target_ref_comms[find_index(qualities, qualities.back() * UNIFORM_DIST(RANDOM_ENGINE))];
        }

        for(auto [u, w]: graph_ptr->adjecencies(v)) {
            if(partition.get_comm(v) != partition.get_comm(u)) continue;
            ref_comm_to_comm_w[best_comm] += best_comm == refined_partition.get_comm(u) ? -w : w;
        }
        
        refined_partition.change_comm(v, best_comm, k_i_ins[refined_partition.get_comm(v)], k_i_ins[best_comm]);
        
        for(auto [u, _]: graph_ptr->adjecencies(v)) {
            k_i_ins[refined_partition.get_comm(u)] = 0;
        }
    }
    
    std::vector<std::vector<vertex_t>> comm_to_nodes(partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        comm_to_nodes[partition.get_comm(v)].push_back(v);
    }
    
    std::vector<community_t> merged_comms(graph_ptr->get_vcount());
    for(const auto& vertices: comm_to_nodes) {
        community_t merging_comm = 0;
        bool comm_merge_enabled = false;
        for(auto v: vertices) {
            community_t ref_comm = refined_partition.get_comm(v);
            merged_comms[v] = ref_comm;
            if(refined_partition.get_comm_size(ref_comm) == 1) {
                if(comm_merge_enabled) {
                    merged_comms[v] = merging_comm;
                } else {
                    merging_comm = ref_comm;
                    comm_merge_enabled = true;
                }
            }
        }
    }

    partition = Partition(graph_ptr, std::move(merged_comms));
}

void run(std::vector<Partition>& hierarchy, std::shared_ptr<const Graph> graph_ptr, std::vector<community_t> initial_vertex_comms, int max_level) {
    std::cout << "Nodes: " << graph_ptr->get_vcount() << ", Edges: " << graph_ptr->get_ecount(); std::cout.flush();
    if (max_level == 0) return;

    hierarchy.emplace_back(graph_ptr, std::move(initial_vertex_comms));
    auto& partition = hierarchy.back();
    std::cout << ", Quality: " << Modularity::quality(partition); std::cout.flush();

    auto cp = get_time();
    bool any_change_occured = move_vertices(partition); partition.renumber();
    std::cout << ", Moving nodes: " << get_elapsed_time(cp, get_time()); std::cout.flush();
    if(!any_change_occured) return;
    
    cp = get_time();
    std::vector<community_t> non_refined_comms = partition.get_vertex_comms();
    refine_partition(partition); partition.renumber();
    std::vector<community_t> comms(partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        comms[partition.get_comm(v)] = non_refined_comms[v];
    }
    std::cout << ", Refine partition: " << get_elapsed_time(cp, get_time()); std::cout.flush();

    cp = get_time();
    std::shared_ptr<const Graph> new_graph_ptr(aggregate_communities(partition));
    std::cout << ", Aggregate communities: " << get_elapsed_time(cp, get_time()) << std::endl;
    run(hierarchy, new_graph_ptr, std::move(comms), max_level - 1);
}
}

std::vector<Partition> run_leiden(std::shared_ptr<const Graph> graph_ptr, int max_level) {
    std::vector<Partition> hierarchy;
    run(hierarchy, graph_ptr, fill_with_communities(graph_ptr->get_vcount()), max_level);
    return hierarchy;
}