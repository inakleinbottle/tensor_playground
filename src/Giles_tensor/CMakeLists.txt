cmake_minimum_required(VERSION 3.16)


project(Giles_free_tensor)


add_library(TensorPlayground_Giles_tensor STATIC)

find_package(OpenMP COMPONENTS CXX REQUIRED)

target_sources(TensorPlayground_Giles_tensor
        PUBLIC
            Giles_tensor.h
        PRIVATE
            Giles_tensor.cpp
        )

target_link_libraries(TensorPlayground_Giles_tensor
        PRIVATE
            Playground::helpers
        PUBLIC
            OpenMP::OpenMP_CXX)
target_include_directories(TensorPlayground_Giles_tensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(TensorPlayground_Giles_tensor PRIVATE -fno-omit-frame-pointer -p -pg)

add_library(Playground::Giles_tensor ALIAS TensorPlayground_Giles_tensor)