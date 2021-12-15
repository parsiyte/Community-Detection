#ifndef COMMUNITIES_FINDER_HPP
#define COMMUNITIES_FINDER_HPP

#include <cstdint>          // uint32_t, uint64_t
#include "Graph.hpp"
#include "Partition.hpp"
#include "Modularity.hpp"

class CommunitiesFinder {
protected:
    Partition partition;
public:
    CommunitiesFinder(const Graph& graph) : partition(graph) {}
    Graph aggregate_communities() const noexcept;
    float quality() const noexcept { return Modularity::quality(partition); }
};


class LouvainAlgorithm : public CommunitiesFinder {
public:
    LouvainAlgorithm(const Graph& graph): CommunitiesFinder(graph) {};
    bool move_vertices() noexcept;
};
class LeidenAlgorithm : public CommunitiesFinder {
    const float res;
public:
    LeidenAlgorithm(const Graph& graph, float res): CommunitiesFinder(graph), res(res) {};
    bool move_vertices() noexcept;
    // void refine_partition() noexcept;
};
#endif