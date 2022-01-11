//
// Created by sam on 07/11/2021.
//

#include <benchmark/benchmark.h>
#include <free_tensor.h>



static void BM_free_tensor_mul(benchmark::State& state)
{
    free_tensor t1(5, 12, {1.0, 1.0, 2.0, 3.0, 4.0, 5.0});
    free_tensor t2(5, 12, {1.0, 0.0, 2.0, -3.0, 4.0, 5.0});

    for (auto _ : state) {
        auto result = t1 * t2;
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }

}

BENCHMARK(BM_free_tensor_mul);