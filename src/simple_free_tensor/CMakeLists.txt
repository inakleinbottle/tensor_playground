cmake_minimum_required(VERSION 3.16)


project(simple_free_tensor)


add_library(TensorPlayground_free_tensor STATIC)

find_package(OpenMP COMPONENTS CXX REQUIRED)

target_sources(TensorPlayground_free_tensor
        PUBLIC
            free_tensor.h
        PRIVATE
            free_tensor.cpp
        )

target_link_libraries(TensorPlayground_free_tensor
        PRIVATE
            Playground::helpers
        PUBLIC
            OpenMP::OpenMP_CXX)
target_include_directories(TensorPlayground_free_tensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(TensorPlayground_free_tensor PRIVATE -fno-omit-frame-pointer -p -pg)
add_library(Playground::free_tensor ALIAS TensorPlayground_free_tensor)