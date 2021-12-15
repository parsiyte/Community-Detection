#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <cstdint>          // uint32_t, uint64_t
#include <utility>          // pair
#include <tuple>
#include <vector>

class Graph {
    float total_weight;
    std::vector<float> v_to_w;
    std::vector<float> v_to_self_w;
    /*
        Sparse Matrix Code Dependence Analysis Simplification at Compile Time - Scientific Figure on ResearchGate. Available from: 
        https://www.researchgate.net/figure/Compressed-Sparse-Row-CSR-sparse-matrix-format-The-val-array-stores-the-nonzeros-by_fig1_326696933 [accessed 12 Nov, 2021]
    */
    std::vector<uint64_t> rowptr;
    std::vector<std::pair<uint32_t, float>> col;
public:
    class Builder;
    friend class Builder;

    bool has_vertex(uint32_t v) const noexcept { return v < get_vcount(); }
    float get_self_weight(uint32_t v) const { return v_to_self_w.at(v); }
    float get_weight(uint32_t v) const { return v_to_w.at(v); };
    float get_weight() const noexcept { return total_weight; }

    uint32_t get_vcount() const noexcept { 
        return rowptr.size() - 1; /* First element isn't counted since we added it for convenience*/ 
    }
    uint64_t get_ecount() const noexcept { return rowptr.back(); }

    class Adjecencies {
    public:
        Adjecencies(const Graph& g, uint32_t v) : beg_iter(g, g.rowptr[v]), end_iter(g, g.rowptr[v + 1]) { }
        auto begin() noexcept { return beg_iter; }
        auto end() noexcept { return end_iter; }
    private:
        class EdgeIter {
        public:
            EdgeIter(const Graph& g, uint64_t beg_indx) : g(g), curr_indx(beg_indx) {};
            std::pair<uint32_t, float> operator*() const noexcept { return g.col.at(curr_indx); }
            EdgeIter& operator++() noexcept { curr_indx++; return *this; }  
            EdgeIter operator++(int) noexcept { EdgeIter tmp = *this; ++(*this); return tmp; }
            friend bool operator== (const EdgeIter& a, const EdgeIter& b) { return a.curr_indx == b.curr_indx; };
            friend bool operator!= (const EdgeIter& a, const EdgeIter& b) { return a.curr_indx != b.curr_indx; };
        private:
            const Graph& g;
            uint64_t curr_indx;
        };
        EdgeIter beg_iter, end_iter;
    };

    /*
        can be used in range-based for loop: 
        for(auto [u, w]: grap.adjecencies(v)) {
            .......
        }
    */
    auto adjecencies(uint32_t v) const  {
        if(!has_vertex(v)) throw std::out_of_range("No vertex found!");
        return Adjecencies(*this, v);
    }
};

class Graph::Builder {
    bool is_built = false;
    std::vector<std::tuple<uint32_t, uint32_t, float>> edges;
public:
    Graph::Builder& insert(uint32_t v1, uint32_t v2, float weight=1.0f);
    Graph build();
};
#endif
