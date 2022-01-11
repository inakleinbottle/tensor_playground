cmake_minimum_required(VERSION 3.16)

project(simple_templated_tensor VERSION 0.0.0)




add_library(TensorPlayground_templated_tensor INTERFACE)
target_sources(TensorPlayground_templated_tensor
        PUBLIC
        simple_template_tensor.h
        )

target_include_directories(TensorPlayground_templated_tensor INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TensorPlayground_templated_tensor INTERFACE Playground::helpers)
add_library(Playground::templated_tensor ALIAS TensorPlayground_templated_tensor)