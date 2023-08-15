#include <stdio.h>
#define WARPSIZE 8
static constexpr int maxWarpsPerBlock = 64 / WARPSIZE;

inline __device__ float atomic_add_bits(float* address, float val)
{
  return atomicAdd(address, val);
}

__inline__ __device__ float warp_reduce_sum(float val)
{
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    val += __shfl_down(val, offset);
  }
  return val;
}

__inline__ __device__ float block_reduce_sum(float val)
{
  // Shared memory for storing the results of each warp-wise partial
  // reduction
  __shared__ float shared[::maxWarpsPerBlock];

  int lane   = threadIdx.x % warpSize;  // thread ID within the warp,
  int warpId = threadIdx.x / warpSize;  // ID of the warp itself

  val = warp_reduce_sum(val);  // Each warp performs partial reduction

  if (lane == 0) {
    shared[warpId] = val;
  }  // Write reduced value to shared memory

  __syncthreads();  // Wait for all partial reductions

  // read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

  if (warpId == 0) {
    val = warp_reduce_sum(val);
  }  // Final reduce within first warp

  return val;
}

__inline__ __device__ void grid_reduce_sum(float val, float* out)
{
  // Reduce the entire block in parallel
  val = block_reduce_sum(val);

  // Write block level reduced value to the output scalar atomically
  if (threadIdx.x == 0) {
    atomic_add_bits(out, val);
  }
}

__global__ void kernel_reduce_sum(float* in, float* out, size_t N);

__global__ void kernel_reduce_sum(float* in, float* out, size_t N) {
  float __shared__ sum_stride[128]; // should be the size of a block

  printf("%d, %d, %d, %d\n", threadIdx.x, blockIdx.x, blockDim.x, gridDim.x);
  // Grid stride loop to perform as much of the reduction as possible
  float local_sum = 0.0;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x) {
    local_sum += in[i];
  }

  sum_stride[threadIdx.x] = local_sum;

  __syncthreads;

  grid_reduce_sum(, out);
  printf("value kernel %d: %f\n", blockIdx.x * blockDim.x + threadIdx.x, sum_val);

}


int main()
{
  const size_t N = pow(8, 2);
  size_t size = N * sizeof(float);

  printf("Checkpoint 1\n");

  // Allocate input vectors h_A and h_B in host memory
  float* host_in = (float*)malloc(size);
  float* host_out = (float*)malloc(1);

  printf("Checkpoint 2\n");

  // Initialize input vectors
  host_out = 0;
  for (int i = 0; i < N; i ++) {
    host_in[i] = 1.0;
  }

  printf("Checkpoint 3\n");

  // Allocate vectors in device memory
  float* dev_in;
  cudaMalloc(&dev_in, size);
  float* dev_out;
  cudaMalloc(&dev_out, 1);

  printf("Checkpoint 4\n");
  
  // Copy vectors from host memory to device memory
  cudaMemcpy(dev_in, host_in, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_out, host_out, 1, cudaMemcpyHostToDevice);

  printf("Checkpoint 5\n");

  size_t threadsPerBlock = 128;
  size_t blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  kernel_reduce_sum<<<blocksPerGrid, threadsPerBlock>>>(dev_in, dev_out, N);
  cudaDeviceSynchronize();

  printf("Checkpoint 6\n");

  cudaMemcpy(host_out, dev_out, 1, cudaMemcpyDeviceToHost);

  printf("Checkpoint 7\n");

  printf("value: %f\n", dev_out);
  
  return 0;
}