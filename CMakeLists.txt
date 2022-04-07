cmake_minimum_required(VERSION 3.16)
project(community-detection LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-fopenmp -fconcepts")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)
add_executable(${PROJECT_NAME}  src/Graph.cpp 
                                src/Partition.cpp
                                src/graph_loader.cpp
                                src/communities_finders.cpp
                                src/louvain_algorithm.cpp
                                src/leiden_algorithm.cpp
                                src/parallel_leiden_algorithm.cpp
                                main.cpp)
