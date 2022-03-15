cmake_minimum_required(VERSION 3.19)

include(FetchContent)
include(GNUInstallDirs)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.8.0
)

include_directories(
        ${PROJECT_SOURCE_DIR}/src
)

FetchContent_MakeAvailable(googletest)
add_library(GTest::gtest ALIAS gtest)
add_library(GTest::Main ALIAS gtest_main)

include(GoogleTest)


set(TEST_FILES 
simple-test.cpp
)


add_executable(${CMAKE_PROJECT_NAME}_test ${TEST_FILES})
target_link_libraries(${CMAKE_PROJECT_NAME}_test gtest_main
        Playground::templated_tensor
        Playground::giles_templated_tensor
        )

include(GoogleTest)
gtest_discover_tests(${CMAKE_PROJECT_NAME}_test)