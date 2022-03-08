cmake_minimum_required(VERSION 3.16)

project(tensor_playground_tests)

enable_testing()

find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

enable_testing()
add_executable(TensorPlayground_tests)

target_link_libraries(TensorPlayground_tests PRIVATE
        GTest::gtest
        GTest::gtest_main
        Playground::templated_tensor
        Playground::giles_templated_tensor
        )

target_sources(TensorPlayground_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/test_giles_tensor.cpp
        )

gtest_discover_tests(TensorPlayground_tests)