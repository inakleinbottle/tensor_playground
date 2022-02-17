cmake_minimum_required(VERSION 3.16)

project(tensor_playground_helpers)




add_library(TensorPlayground_helpers INTERFACE  increasing_degree_walker.h decreasing_degree_walker.h index_word.h)
target_include_directories(TensorPlayground_helpers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


add_library(Playground::helpers ALIAS TensorPlayground_helpers)