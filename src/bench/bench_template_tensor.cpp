//
// Created by sam on 07/11/2021.
//

#include <benchmark/benchmark.h>
#include <simple_template_tensor.h>
#include "params.h"

using namespace playground;

static void BM_templated_tensor_mul(benchmark::State &state) {
    simple_template_tensor<WIDTH, DEPTH> t1{1.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    simple_template_tensor<WIDTH, DEPTH> t2{1.0, 0.0, -2.0, 3.0, -4.0, -5.0};

    for (auto _ : state) {
        auto result = t1 * t2;
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK(BM_templated_tensor_mul);