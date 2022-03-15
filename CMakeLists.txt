cmake_minimum_required(VERSION 3.20)
project(tensor_playground)

set(CMAKE_CXX_STANDARD 17)

add_executable(tensor_playground main.cpp src/Giles_tensor/Giles_tensor.cpp src/Giles_tensor/Giles_tensor.h src/template-Giles-tensor/giles_template_tensor.h)



add_subdirectory(src)

target_link_libraries(tensor_playground PUBLIC Playground::templated_tensor)
target_link_libraries(tensor_playground PRIVATE Playground::free_tensor)

enable_testing()
add_subdirectory(tests)