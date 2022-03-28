#include"dust_cuda.h"
#include<math.h>
#include<vector.h>
#include"/Users/helenraichie/GitHub/cholla/src/global/global.h"
#include"/Users/helenraichie/GitHub/cholla/src/global/global_cuda.h"
#include"/Users/helenraichie/GitHub/cholla/src/global/gpu.hpp"

__global__ void cuda_hello() {
    printf("Hello World from GPU!\n");
}

int main() {
    cuda_hello<<<1, 1>>>();
    return 0;
}