cmake_minimum_required(VERSION 3.16)
project(community-detection LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-pg")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)
add_executable(${PROJECT_NAME} src/Graph.cpp src/Partition.cpp src/Modularity.cpp src/GraphLoader.cpp src/CommunitiesFinder.cpp main.cpp)