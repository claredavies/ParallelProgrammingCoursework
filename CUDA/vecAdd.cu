#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(float *real_rasc, float *real_decl, float *rand_rasc, float *rand_decl,float *histogram_DD
    ,float *histogram_DR, float *histogram_RR, int n)
{
    // Get our global thread ID
    int id = blockIdx.x*blockDim.x+threadIdx.x;

    printf("id:   %d    N:  %d \n",id,n);
    // Make sure we do not go out of bounds
    if (id < n) {
        histogram_DD[id] = real_rasc[id] + real_decl[id] + rand_rasc[id] + rand_decl[id];
        histogram_DR[id] = real_rasc[id] + real_decl[id] + rand_rasc[id] + rand_decl[id];
        histogram_RR[id] = real_rasc[id] + real_decl[id] + rand_rasc[id] + rand_decl[id];
    }
}
 
int main( int argc, char* argv[] )
{
    // Size of vectors
    int n = 100000;
 
    // Host input vectors
    float *h_real_rasc;
    float *h_real_decl;
    float *h_rand_rasc;
    float *h_rand_decl;

    //Host output vector
    float *h_histogram_DD;
    float *h_histogram_DR;
    float *h_histogram_RR;

    // Device input vectors
    float *d_real_rasc;
    float *d_real_decl;
    float *d_rand_rasc;
    float *d_rand_decl;

    //Device output vector
     float *d_histogram_DD;
     float *d_histogram_DR;
     float *d_histogram_RR;
 
    // Size, in bytes, of each vector
    size_t bytes = n*sizeof(double);
 
    // Allocate memory for each vector on host
    h_real_rasc = (float*)malloc(bytes);
    h_real_decl = (float*)malloc(bytes);
    h_rand_rasc = (float*)malloc(bytes);
    h_rand_decl = (float*)malloc(bytes);

    h_histogram_DD = (float*)malloc(bytes);
    h_histogram_DR = (float*)malloc(bytes);
    h_histogram_RR = (float*)malloc(bytes);

    // Allocate memory for each vector on GPU
    cudaMalloc(&d_real_rasc, bytes);
    cudaMalloc(&d_real_decl, bytes);
    cudaMalloc(&d_rand_rasc, bytes);
    cudaMalloc(&d_rand_decl, bytes);
    cudaMalloc(&d_histogram_DD, bytes);
    cudaMalloc(&d_histogram_DR, bytes);
    cudaMalloc(&d_histogram_RR, bytes);
 
    int i;
    // Initialize vectors on host
    for( i = 0; i < n; i++ ) {
        h_real_rasc[i] = 0.25;
        h_real_decl[i] = 0.25;
        h_rand_rasc [i] = 0.25;
        h_rand_decl[i] = 0.25;
    }
 
    // Copy host vectors to device
    cudaMemcpy(d_real_rasc, h_real_rasc, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_real_decl, h_real_decl, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_rasc, h_rand_rasc, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_decl, h_rand_decl, bytes, cudaMemcpyHostToDevice);

    int blockSize, gridSize;
 
    // Number of threads in each thread block
    blockSize = 1024;
 
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)n/blockSize);
 
    // Execute the kernel
    vecAdd<<<gridSize, blockSize>>>(d_real_rasc, d_real_decl, d_rand_rasc, d_rand_decl, d_histogram_DD, d_histogram_DR
     , d_histogram_RR, n);
 
    // Copy array back to host
    cudaMemcpy( h_histogram_DD, d_histogram_DD, bytes, cudaMemcpyDeviceToHost );
    cudaMemcpy( h_histogram_DR, d_histogram_DR, bytes, cudaMemcpyDeviceToHost );
    cudaMemcpy( h_histogram_RR, d_histogram_RR, bytes, cudaMemcpyDeviceToHost );

    // Sum up vector c and print result divided by n, this should equal 1 within error
    double sum_DD = 0;
    double sum_DR = 0;
    double sum_RR = 0;
    for(i=0; i<n; i++) {
         sum_DD += h_histogram_DD[i];
         sum_DR += h_histogram_DR[i];
         sum_RR += h_histogram_RR[i];
    }
    printf("final result DD: %f  final result DR: %f  final result RR: %f\n", sum_DD,sum_DR,sum_RR);
 
    // Release device memory
    cudaFree(d_real_rasc);
    cudaFree(d_real_decl);
    cudaFree(d_rand_rasc);
    cudaFree(d_rand_decl);
    cudaFree(d_histogram_DD);
    cudaFree(d_histogram_DR);
    cudaFree(d_histogram_RR);

 
    // Release host memory
    free(h_real_rasc);
    free(h_real_decl);
    free(h_rand_rasc);
    free(h_rand_decl);
    free(h_histogram_DD);
    free(h_histogram_DR);
    free(h_histogram_RR);
 
    return 0;
}