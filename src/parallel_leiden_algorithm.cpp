#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <deque>
#include <random>               // shuffle
#include <algorithm>            // iota
#include <omp.h>
#include <mutex>

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

    std::vector<vertex_t> vertices(graph_ptr->get_vcount());
    std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<bool> v_in_queue(graph_ptr->get_vcount());
    std::vector<vertex_t> new_vertices;

    bool any_changes_occured = false;
    #pragma omp parallel 
    while(!vertices.empty()) {
        using changing_elem = std::tuple<vertex_t, community_t, float, float>;
        std::vector<changing_elem> changing_elems;
        std::vector<vertex_t> adj_changed_vertices;
        
        #pragma omp for
        for(int i = 0; i < vertices.size(); i++) {
            const vertex_t v = vertices[i];
            
            std::unordered_map<community_t, float> k_i_ins;
            k_i_ins[partition.get_comm(v)] = 0;
            for(auto [u, _]: graph_ptr->adjecencies(v)) {
                k_i_ins[partition.get_comm(u)] = 0;    
            }

            for(auto [u, w]: graph_ptr->adjecencies(v)) {
                if(v == u) continue;
                k_i_ins[partition.get_comm(u)] += w;
            }

            community_t initial_community = partition.get_comm(v);
            community_t best_community = Modularity::find_best_community(v, partition, k_i_ins);
            if (initial_community != best_community) {
                partition.change_comm(v, best_community, k_i_ins[initial_community], k_i_ins[best_community], true);

                any_changes_occured = true;
                for(auto [u, _]: graph_ptr->adjecencies(v)) {
                    if(!v_in_queue[u] && partition.get_comm(u) != best_community){
                        v_in_queue[u] = true;
                        adj_changed_vertices.push_back(u);
                    } 
                }
            }
        }

        size_t new_vertices_beg_idx = 0;
        #pragma omp critical
        {
            new_vertices_beg_idx = new_vertices.size();
            new_vertices.resize(new_vertices_beg_idx + adj_changed_vertices.size());
        }
        std::copy(adj_changed_vertices.begin(), adj_changed_vertices.end(), new_vertices.begin() + new_vertices_beg_idx);

        #pragma omp barrier
        #pragma omp single // has implicit barrier
        {
            v_in_queue.resize(graph_ptr->get_vcount(), false);
            vertices = std::move(new_vertices);
            new_vertices.clear();
        }
    }

    return any_changes_occured;
}

void refine_partition(Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();

    std::vector<vertex_t> renumbered_vertex(graph_ptr->get_vcount());
    std::vector<std::vector<vertex_t>> comms(partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        renumbered_vertex[v] = comms[partition.get_comm(v)].size();
        comms[partition.get_comm(v)].push_back(v);
    }

    Partition refined_partition(graph_ptr);
    std::vector<float> ref_comm_to_comm_w(graph_ptr->get_vcount());
    #pragma omp parallel for
    for(const auto& vertices: comms) {
        for(vertex_t v: vertices) {
            for(auto [u, w]: graph_ptr->adjecencies(v)) {
                if(partition.get_comm(v) != partition.get_comm(u)) continue;
                ref_comm_to_comm_w[refined_partition.get_comm(v)] += w;
            }
        }

        std::vector<float> k_i_ins(vertices.size());
        for(vertex_t v: vertices) {
            community_t v_ref_comm = refined_partition.get_comm(v);
            if(refined_partition.get_comm_size(v_ref_comm) > 1) continue;
            
            float comm_weight = partition.get_comm_weight(partition.get_comm(v));
            if(ref_comm_to_comm_w[v_ref_comm] < graph_ptr->get_weight(v) * (comm_weight - graph_ptr->get_weight(v)) / graph_ptr->get_weight()) continue;
            ref_comm_to_comm_w[v_ref_comm] = 0;

            std::vector<community_t> target_ref_comms;
            target_ref_comms.push_back(v_ref_comm);
            for(auto [u, w]: graph_ptr->adjecencies(v)) {
                if(partition.get_comm(v) != partition.get_comm(u)) continue;
                community_t u_ref_comm = refined_partition.get_comm(u);
                if(ref_comm_to_comm_w[u_ref_comm] < refined_partition.get_comm_weight(u_ref_comm) * (comm_weight - refined_partition.get_comm_weight(u_ref_comm)) / graph_ptr->get_weight()) continue;
                if(is_equal(k_i_ins[renumbered_vertex[u_ref_comm]], 0)) {
                    target_ref_comms.push_back(u_ref_comm);
                } 
                k_i_ins[renumbered_vertex[u_ref_comm]] += w;
            }

            community_t best_comm = v_ref_comm;
            double max_gain = 0, cumm_exp = 0;
            std::vector<double> qualities;
            for(auto ref_comm: target_ref_comms){
                float gain = Modularity::compute_delta_Q(graph_ptr->get_weight(), refined_partition.get_comm_weight(ref_comm), graph_ptr->get_weight(v), k_i_ins[renumbered_vertex[ref_comm]]);
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
            
            refined_partition.change_comm(v, best_comm, k_i_ins[renumbered_vertex[v_ref_comm]], k_i_ins[renumbered_vertex[best_comm]]);

            for(int i = 0; i < k_i_ins.size(); i++) k_i_ins[i] = 0;
        }
    }

    partition = std::move(refined_partition);
}

Graph* aggregate_communities(const Partition& partition) {
    auto graph_ptr = partition.get_graph_ptr();
    
    std::vector<std::vector<vertex_t>> comm_to_nodes(partition.get_ccount());
    for(vertex_t v = 0; v < graph_ptr->get_vcount(); v++) {
        comm_to_nodes[partition.get_comm(v)].push_back(v);
    }
    
    Graph::Builder builder;
    #pragma omp parallel 
    {   
        using edge = std::tuple<vertex_t, vertex_t, float>;
        std::vector<edge> new_edges;

        #pragma omp for
        for(community_t c_from = 0; c_from < comm_to_nodes.size(); c_from++){        
            std::unordered_map<community_t, float> weights_cache;
            std::unordered_set<community_t> looked_comms;
            for(auto v_from: comm_to_nodes[c_from]) {
                for(auto [v_to, w]: graph_ptr->adjecencies(v_from)) {
                    community_t c_to = partition.get_comm(v_to);
                    looked_comms.insert(c_to);
                    weights_cache[c_to] = 0;
                }
            }

            for(auto v_from: comm_to_nodes[c_from]) {
                for(auto [v_to, w]: graph_ptr->adjecencies(v_from)) {
                    weights_cache[partition.get_comm(v_to)] += w;
                }
            }

            for(auto c_to: looked_comms){
                new_edges.emplace_back(c_from, c_to, weights_cache[c_to]);
                weights_cache[c_to] = 0;
            }
        }

        #pragma omp critical
        for(auto [from, to, weight]: new_edges) {
            builder.insert(from, to, weight);
        }
    }
    
    return builder.build();
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
    std::shared_ptr<const Graph> new_graph_ptr(::aggregate_communities(partition));
    std::cout << ", Aggregate communities: " << get_elapsed_time(cp, get_time()) << std::endl;
    run(hierarchy, new_graph_ptr, std::move(comms), max_level - 1);
}
}

std::vector<Partition> run_parallel_leiden(std::shared_ptr<const Graph> graph_ptr, int max_level) {
    std::vector<Partition> hierarchy;
    run(hierarchy, graph_ptr, fill_with_communities(graph_ptr->get_vcount()), max_level);
    return hierarchy;
}