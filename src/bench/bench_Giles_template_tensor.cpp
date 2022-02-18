//
// Created by sam on 14/02/2022.
//

#include <benchmark/benchmark.h>
#include <giles_template_tensor.h>

using namespace playground;

static void BM_Giles_templated_tensor_mul(benchmark::State &state) {
    giles_template_tensor<5, 8> t1{1.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    giles_template_tensor<5, 8> t2{1.0, 0.0, -2.0, 3.0, -4.0, -5.0};

    for (auto _ : state) {
        auto result = t1 * t2;
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK(BM_Giles_templated_tensor_mul);