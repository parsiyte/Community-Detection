#include "GraphLoader.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

Graph load_graph(std::string path) {
    std::ifstream finput;
    finput.open(path, std::fstream::in);
    if(!finput.is_open() || finput.fail()) {
        std::cerr << "Error while opening " << path << std::endl;
        exit(EXIT_FAILURE);
    }
    
    Graph::Builder builder;
    uint32_t v, u;
    while(finput >> v >> u) {
        builder.insert(v, u);
        builder.insert(u, v);
    }
    finput.close();
    return builder.build();
}