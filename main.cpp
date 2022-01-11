#include <iostream>

#include "src/simple_templated_tensor/simple_template_tensor.h"
#include "src/simple_free_tensor/free_tensor.h"

int main()
{
    free_tensor tensor1 (5, 5,  {1.0, 0.0, 2.0});
    free_tensor tensor2 (5, 5,  {1.0, 1.0, 0.0});

    auto result1 = tensor1*tensor2;


    simple_template_tensor<5, 5> tt1{1.0, 0.0, 2.0};
    simple_template_tensor<5, 5> tt2{1.0, 1.0, 0.0};

    auto result2 = tt1*tt2;

    std::cout << result1.size() << ' ' << result2.size() << '\n';

    return 0;
}
