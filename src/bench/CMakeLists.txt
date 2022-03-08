cmake_minimum_required(VERSION 3.16)


project(tensor_playground_benchmarks)

find_package(benchmark REQUIRED)


add_executable(TensorPlayground_benchmarks params.h)

target_sources(TensorPlayground_benchmarks
        PRIVATE
        bench_template_tensor.cpp
        bench_free_tensor.cpp
        bench_Giles_template_tensor.cpp
        )

target_link_libraries(TensorPlayground_benchmarks PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        Playground::templated_tensor
        Playground::free_tensor
        Playground::giles_templated_tensor
        )
target_compile_options(TensorPlayground_benchmarks PRIVATE -fno-omit-frame-pointer -g -p -pg)