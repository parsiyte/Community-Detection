#ifndef COMMUNITIES_FINDERS_HPP
#define COMMUNITIES_FINDERS_HPP

#include <vector>
#include <memory>           // shared_ptr

#include "Graph.hpp"
#include "Partition.hpp"

Graph* aggregate_communities(const Partition& partition);
std::vector<Partition> run_louvain(std::shared_ptr<const Graph> graph_ptr, int max_level=-1);
std::vector<Partition> run_leiden(std::shared_ptr<const Graph> graph_ptr, int max_level=-1);
#endif