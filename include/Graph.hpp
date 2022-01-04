#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <utility>          // pair

#include "types.hpp"

class Graph {
    Graph() {};
    float total_weight = 0;
    std::vector<float> vertex_weights;
    std::vector<float> vertex_self_weights;
    /*
        Sparse Matrix Code Dependence Analysis Simplification at Compile Time - Scientific Figure on ResearchGate. Available from: 
        https://www.researchgate.net/figure/Compressed-Sparse-Row-CSR-sparse-matrix-format-The-val-array-stores-the-nonzeros-by_fig1_326696933 [accessed 12 Nov, 2021]
    */
    std::vector<size_t> rowptr;
    std::vector<std::pair<vertex_t, float>> col;
public:
    class Builder;
    friend class Builder;

    bool has_vertex(vertex_t v) const noexcept { return v < get_vcount(); }
    float get_self_weight(vertex_t v) const { return vertex_self_weights.at(v); }
    float get_weight(vertex_t v) const { return vertex_weights.at(v); };
    float get_weight() const noexcept { return total_weight; }

    size_t get_vcount() const noexcept { return vertex_weights.size(); }
    size_t get_ecount() const noexcept { return rowptr.back(); }

    class Adjecencies {
    public:
        Adjecencies(const Graph& g, vertex_t v) : beg_iter(g, g.rowptr[v]), end_iter(g, g.rowptr[v + 1]) { }
        auto begin() noexcept { return beg_iter; }
        auto end() noexcept { return end_iter; }
    private:
        class EdgeIter {
        public:
            EdgeIter(const Graph& g, size_t beg_indx) : g(g), curr_indx(beg_indx) {};
            std::pair<vertex_t, float> operator*() const noexcept { return g.col.at(curr_indx); }
            EdgeIter& operator++() noexcept { curr_indx++; return *this; }  
            EdgeIter operator++(int) noexcept { EdgeIter tmp = *this; ++(*this); return tmp; }
            friend bool operator== (const EdgeIter& a, const EdgeIter& b) { return a.curr_indx == b.curr_indx; };
            friend bool operator!= (const EdgeIter& a, const EdgeIter& b) { return a.curr_indx != b.curr_indx; };
        private:
            const Graph& g;
            size_t curr_indx;
        };
        EdgeIter beg_iter, end_iter;
    };

    /*
        can be used in range-based for loop: 
        for(auto [u, w]: grap.adjecencies(v)) {
            .......
        }
    */
    auto adjecencies(vertex_t v) const  {
        if(!has_vertex(v)) throw std::out_of_range("No vertex found!");
        return Adjecencies(*this, v);
    }
};

class Graph::Builder {
    bool is_built = false;
    std::vector<edge_t> edges;
public:
    Graph::Builder& insert(vertex_t v_from, vertex_t v_to, float weight=1.0f);
    Graph* build();
};
#endif
