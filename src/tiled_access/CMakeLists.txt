cmake_minimum_required(VERSION 3.16)

project(tensor_tiled_access VERSION 0.0.0)

add_library(TensorPlayground_tensor_tiled_access INTERFACE strided_rw.h)
target_sources(TensorPlayground_tensor_tiled_access INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}/implementation.h
        ${CMAKE_CURRENT_LIST_DIR}/tile.h
        )

target_include_directories(TensorPlayground_tensor_tiled_access INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TensorPlayground_tensor_tiled_access INTERFACE Playground::helpers)
add_library(Playground::tensor_tiled_access ALIAS TensorPlayground_tensor_tiled_access)