cmake_minimum_required(VERSION 3.10)
project(VectorFieldAnalysis VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(vector_field STATIC src/vector_field.cpp)
target_include_directories(vector_field PUBLIC include)

add_executable(example_1 examples/example_1.cpp)
target_link_libraries(example_1 vector_field)