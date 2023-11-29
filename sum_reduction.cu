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

  // Write block level reduced value to the output scalar automically
  if (threadIdx.x == 0) {
    atomic_add_bits(out, val);
  }
}

__global__ void kernel_reduce_sum(float* in, float* out, size_t N);

__global__ void kernel_reduce_sum(float* in, float* out, size_t N) {

  // __shared__ float sum_stride[128];  // array shared between each block

  // Grid stride loop to read global array into shared block-wide array
  float sum_stride = 0;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x) {
    // sum_stride[threadIdx.x] += in[i];
    sum_stride += in[i];
  }

  __syncthreads();

  // grid_reduce_sum(sum_stride[threadIdx.x], out);
  grid_reduce_sum(sum_stride, out);

}


int main()
{
  const size_t N = 512;
  size_t size = N * sizeof(float);
  size_t fsize = 1 * sizeof(float);

  // Allocate input vectors h_A and h_B in host memory
  float* host_in = (float*)malloc(size);
  float* host_out = (float*)malloc(1);

  // Initialize input vectors
  *host_out = 0;
  for (int i = 0; i < N; i ++) {
    host_in[i] = 1.0;
  }

  // Allocate vectors in device memory
  float* dev_in;
  cudaMalloc((void**)&dev_in, size);
  float* dev_out;
  cudaMalloc((void**)&dev_out, 1*sizeof(float));
  
  // Copy vectors from host memory to device memory
  cudaMemcpy(dev_in, host_in, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_out, host_out, fsize, cudaMemcpyHostToDevice);

  size_t threadsPerBlock = 128;
  size_t blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  kernel_reduce_sum<<<blocksPerGrid, threadsPerBlock>>>(dev_in, dev_out, N);
  cudaDeviceSynchronize();

  cudaMemcpy(host_out, dev_out, 1*sizeof(float), cudaMemcpyDeviceToHost);

  printf("Result: %f\n", *host_out);

  free(host_in);
  free(host_out);
  cudaFree(dev_in);
  cudaFree(dev_out);
  
  return 0;
}


// cudaError_t err = cudaGetLastError();
// if (err != cudaSuccess) {
//     printf("CUDA error: %s\n", cudaGetErrorString(err));
//     exit(-1);
// }