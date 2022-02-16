cmake_minimum_required(VERSION 3.16)

project(giles_templated_tensor VERSION 0.0.0)

add_library(TensorPlayground_giles_template_tensor INTERFACE mul_level_func.h)
target_sources(TensorPlayground_giles_template_tensor INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}/giles_template_tensor.h
        )

target_include_directories(TensorPlayground_giles_template_tensor INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TensorPlayground_giles_template_tensor INTERFACE Playground::helpers Playground::templated_tensor Playground::tensor_tiled_access)
add_library(Playground::giles_templated_tensor ALIAS TensorPlayground_giles_template_tensor)