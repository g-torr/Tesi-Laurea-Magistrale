/*
 * This program uses the device CURAND API to calculate what 
 * proportion of pseudo-random ints have low bit set.
 * It then generates uniform results to calculate how many
 * are greater than .5.
 * It then generates  normal results to calculate how many 
 * are within one standard deviation of the mean.
 */
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__); \
    return EXIT_FAILURE;}} while(0)

__global__ void setup_kernel(curandState *state)
{
    int id = threadIdx.x + blockIdx.x * 64;
    /* Each thread gets same seed, a different sequence 
       number, no offset */
    curand_init(1234, id, 0, &state[id]);
}

__global__ void setup_kernel(curandStatePhilox4_32_10_t *state)
{
    int id = threadIdx.x + blockIdx.x * 64;
    /* Each thread gets same seed, a different sequence 
       number, no offset */
    curand_init(1234, id, 0, &state[id]);
}

__global__ void setup_kernel(curandStateMRG32k3a *state)
{
    int id = threadIdx.x + blockIdx.x * 64;
    /* Each thread gets same seed, a different sequence 
       number, no offset */
    curand_init(0, id, 0, &state[id]);
}

__global__ void generate_kernel(curandState *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    int count = 0;
    unsigned int x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[id];
    /* Generate pseudo-random unsigned ints */
    for(int i = 0; i < n; i++) {
        x = curand(&localState);
        /* Check if low bit set */
        if(x & 1) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_kernel(curandStatePhilox4_32_10_t *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    int count = 0;
    unsigned int x;
    /* Copy state to local memory for efficiency */
    curandStatePhilox4_32_10_t localState = state[id];
    /* Generate pseudo-random unsigned ints */
    for(int i = 0; i < n; i++) {
        x = curand(&localState);
        /* Check if low bit set */
        if(x & 1) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_uniform_kernel(curandState *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    float x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[id];
    /* Generate pseudo-random uniforms */
    for(int i = 0; i < n; i++) {
        x = curand_uniform(&localState);
        /* Check if > .5 */
        if(x > .5) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_uniform_kernel(curandStatePhilox4_32_10_t *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    float x;
    /* Copy state to local memory for efficiency */
    curandStatePhilox4_32_10_t localState = state[id];
    /* Generate pseudo-random uniforms */
    for(int i = 0; i < n; i++) {
        x = curand_uniform(&localState);
        /* Check if > .5 */
        if(x > .5) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_normal_kernel(curandState *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    float2 x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[id];
    /* Generate pseudo-random normals */
    for(int i = 0; i < n/2; i++) {
        x = curand_normal2(&localState);
        /* Check if within one standard deviaton */
        if((x.x > -1.0) && (x.x < 1.0)) {
            count++;
        }
        if((x.y > -1.0) && (x.y < 1.0)) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_normal_kernel(curandStatePhilox4_32_10_t *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    float2 x;
    /* Copy state to local memory for efficiency */
    curandStatePhilox4_32_10_t localState = state[id];
    /* Generate pseudo-random normals */
    for(int i = 0; i < n/2; i++) {
        x = curand_normal2(&localState);
        /* Check if within one standard deviaton */
        if((x.x > -1.0) && (x.x < 1.0)) {
            count++;
        }
        if((x.y > -1.0) && (x.y < 1.0)) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_kernel(curandStateMRG32k3a *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    unsigned int x;
    /* Copy state to local memory for efficiency */
    curandStateMRG32k3a localState = state[id];
    /* Generate pseudo-random unsigned ints */
    for(int i = 0; i < n; i++) {
        x = curand(&localState);
        /* Check if low bit set */
        if(x & 1) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_uniform_kernel(curandStateMRG32k3a *state,
                                int n, 
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    double x;
    /* Copy state to local memory for efficiency */
    curandStateMRG32k3a localState = state[id];
    /* Generate pseudo-random uniforms */
    for(int i = 0; i < n; i++) {
        x = curand_uniform_double(&localState);
        /* Check if > .5 */
        if(x > .5) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

__global__ void generate_normal_kernel(curandStateMRG32k3a *state, 
                                int n,
                                unsigned int *result)
{
    int id = threadIdx.x + blockIdx.x * 64;
    unsigned int count = 0;
    double2 x;
    /* Copy state to local memory for efficiency */
    curandStateMRG32k3a localState = state[id];
    /* Generate pseudo-random normals */
    for(int i = 0; i < n/2; i++) {
        x = curand_normal2_double(&localState);
        /* Check if within one standard deviaton */
        if((x.x > -1.0) && (x.x < 1.0)) {
            count++;
        }
        if((x.y > -1.0) && (x.y < 1.0)) {
            count++;
        }
    }
    /* Copy state back to global memory */
    state[id] = localState;
    /* Store results */
    result[id] += count;
}

int main(int argc, char *argv[])
{

    int i;
    unsigned int total;
    curandState *devStates;
    curandStateMRG32k3a *devMRGStates;
    curandStatePhilox4_32_10_t *devPHILOXStates;
    unsigned int *devResults, *hostResults;
    bool useMRG = 0;
    bool usePHILOX = 0;
    int sampleCount = 10000;
    bool doubleSupported = 0;
    int device;
    struct cudaDeviceProp properties;  

    /* check for double precision support */
    CUDA_CALL(cudaGetDevice(&device));
    CUDA_CALL(cudaGetDeviceProperties(&properties,device));
    if ( properties.major >= 2 || (properties.major == 1 && properties.minor >= 3) ) {
        doubleSupported = 1;
    }
    
    /* Check for MRG32k3a option (default is XORWOW) */
    if (argc >= 2)  {
        if (strcmp(argv[1],"-m") == 0) {
            useMRG = 1;
            if (!doubleSupported){
                printf("MRG32k3a requires double precision\n");
                printf("^^^^ test WAIVED due to lack of double precision\n");
                return EXIT_SUCCESS;
            }
        }else if (strcmp(argv[1],"-p") == 0) {
		usePHILOX = 1;
	} 
        /* Allow over-ride of sample count */    
        sscanf(argv[argc-1],"%d",&sampleCount); 
    }

    /* Allocate space for results on host */
    hostResults = (unsigned int *)calloc(64 * 64, sizeof(int));

    /* Allocate space for results on device */
    CUDA_CALL(cudaMalloc((void **)&devResults, 64 * 64 * 
              sizeof(unsigned int)));

    /* Set results to 0 */
    CUDA_CALL(cudaMemset(devResults, 0, 64 * 64 * 
              sizeof(unsigned int)));

    /* Allocate space for prng states on device */
    if (useMRG) {
        CUDA_CALL(cudaMalloc((void **)&devMRGStates, 64 * 64 * 
                  sizeof(curandStateMRG32k3a)));
    }else if(usePHILOX) {
        CUDA_CALL(cudaMalloc((void **)&devPHILOXStates, 64 * 64 * 
                  sizeof(curandStatePhilox4_32_10_t)));
    }else {
        CUDA_CALL(cudaMalloc((void **)&devStates, 64 * 64 * 
                  sizeof(curandState)));
    }
    
    /* Setup prng states */
    if (useMRG) {
        setup_kernel<<<64, 64>>>(devMRGStates);
    }else if(usePHILOX)
    {
        setup_kernel<<<64, 64>>>(devPHILOXStates);
    }else {
        setup_kernel<<<64, 64>>>(devStates);
    }
    
    /* Generate and use pseudo-random  */
    for(i = 0; i < 50; i++) {
        if (useMRG) {
            generate_kernel<<<64, 64>>>(devMRGStates, sampleCount, devResults);
        }else if (usePHILOX){
            generate_kernel<<<64, 64>>>(devPHILOXStates, sampleCount, devResults);
	}else {
            generate_kernel<<<64, 64>>>(devStates, sampleCount, devResults);
        }
    }
    
    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostResults, devResults, 64 * 64 * 
        sizeof(unsigned int), cudaMemcpyDeviceToHost));

    /* Show result */
    total = 0;
    for(i = 0; i < 64 * 64; i++) {
        total += hostResults[i];
    }
    printf("Fraction with low bit set was %10.13f\n", 
        (float)total / (64.0f * 64.0f * sampleCount * 50.0f));
        
    /* Set results to 0 */
    CUDA_CALL(cudaMemset(devResults, 0, 64 * 64 * 
              sizeof(unsigned int)));

    /* Generate and use uniform pseudo-random  */
    for(i = 0; i < 50; i++) {
        if (useMRG) {
            generate_uniform_kernel<<<64, 64>>>(devMRGStates, sampleCount, devResults);
        }else if(usePHILOX) {
            generate_uniform_kernel<<<64, 64>>>(devPHILOXStates, sampleCount, devResults);
	}else {
            generate_uniform_kernel<<<64, 64>>>(devStates, sampleCount, devResults);
        }
    }

    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostResults, devResults, 64 * 64 * 
        sizeof(unsigned int), cudaMemcpyDeviceToHost));

    /* Show result */
    total = 0;
    for(i = 0; i < 64 * 64; i++) {
        total += hostResults[i];
    }
    printf("Fraction of uniforms > 0.5 was %10.13f\n", 
        (float)total / (64.0f * 64.0f * sampleCount * 50.0f));
    /* Set results to 0 */
    CUDA_CALL(cudaMemset(devResults, 0, 64 * 64 * 
              sizeof(unsigned int)));

    /* Generate and use normal pseudo-random  */
    for(i = 0; i < 50; i++) {
        if (useMRG) {
            generate_normal_kernel<<<64, 64>>>(devMRGStates, sampleCount, devResults);
        }else if(usePHILOX) {
            generate_normal_kernel<<<64, 64>>>(devPHILOXStates, sampleCount, devResults);
	}else {
            generate_normal_kernel<<<64, 64>>>(devStates, sampleCount, devResults);
        }
    }

    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostResults, devResults, 64 * 64 * 
        sizeof(unsigned int), cudaMemcpyDeviceToHost));

    /* Show result */
    total = 0;
    for(i = 0; i < 64 * 64; i++) {
        total += hostResults[i];
    }
    printf("Fraction of normals within 1 standard deviation was %10.13f\n", 
        (float)total / (64.0f * 64.0f * sampleCount * 50.0f));

    /* Cleanup */
    if (useMRG) {
        CUDA_CALL(cudaFree(devMRGStates));
    }else if(usePHILOX)
    {
        CUDA_CALL(cudaFree(devPHILOXStates));
    }else {
        CUDA_CALL(cudaFree(devStates));
    }    
    CUDA_CALL(cudaFree(devResults));
    free(hostResults);
    printf("^^^^ kernel_example PASSED\n");
    return EXIT_SUCCESS;
}

The following example uses the cuRAND host MTGP setup API, and the cuRAND device API, to generate integers using the MTGP32 generator, and calculates the proportion that have the low bit set.

/*
 * This program uses the device CURAND API to calculate what 
 * proportion of pseudo-random ints have low bit set.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>


#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__); \
    return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__); \
    return EXIT_FAILURE;}} while(0)

__global__ void generate_kernel(curandStateMtgp32 *state, 
                                int n,
                                int *result)
{
    int id = threadIdx.x + blockIdx.x * 256;
    int count = 0;
    unsigned int x;
    /* Generate pseudo-random unsigned ints */
    for(int i = 0; i < n; i++) {
        x = curand(&state[blockIdx.x]);
        /* Check if low bit set */
        if(x & 1) {
            count++;
        }
    }
    /* Store results */
    result[id] += count;
}

int main(int argc, char *argv[])
{
    int i;
    long long total;
    curandStateMtgp32 *devMTGPStates;
    mtgp32_kernel_params *devKernelParams;
    int *devResults, *hostResults;
    int sampleCount = 10000;
    
    /* Allow over-ride of sample count */    
    if (argc == 2) {
        sscanf(argv[1],"%d",&sampleCount);
    }
        
    /* Allocate space for results on host */
    hostResults = (int *)calloc(64 * 256, sizeof(int));

    /* Allocate space for results on device */
    CUDA_CALL(cudaMalloc((void **)&devResults, 64 * 256 * 
              sizeof(int)));

    /* Set results to 0 */
    CUDA_CALL(cudaMemset(devResults, 0, 64 * 256 * 
              sizeof(int)));

    /* Allocate space for prng states on device */
    CUDA_CALL(cudaMalloc((void **)&devMTGPStates, 64 * 
              sizeof(curandStateMtgp32)));
    
    /* Setup MTGP prng states */
    
    /* Allocate space for MTGP kernel parameters */
    CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
    
    /* Reformat from predefined parameter sets to kernel format, */
    /* and copy kernel parameters to device memory               */
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
    
    /* Initialize one state per thread block */
    CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, 
                mtgp32dc_params_fast_11213, devKernelParams, 64, 1234));
    
    /* State setup is complete */
    
    /* Generate and use pseudo-random  */
    for(i = 0; i < 10; i++) {
        generate_kernel<<<64, 256>>>(devMTGPStates, sampleCount, devResults);
    }

    /* Copy device memory to host */
    CUDA_CALL(cudaMemcpy(hostResults, devResults, 64 * 256 * 
        sizeof(int), cudaMemcpyDeviceToHost));

    /* Show result */
    total = 0;
    for(i = 0; i < 64 * 256; i++) {
        total += hostResults[i];
    }
    
    
    printf("Fraction with low bit set was %10.13g\n", 
        (double)total / (64.0f * 256.0f * sampleCount * 10.0f));

    /* Cleanup */
    CUDA_CALL(cudaFree(devMTGPStates));
    CUDA_CALL(cudaFree(devResults));
    free(hostResults);
    printf("^^^^ kernel_mtgp_example PASSED\n");
    return EXIT_SUCCESS;
}


Read more at: http://docs.nvidia.com/cuda/curand/index.html#ixzz4MIRQ90YK
Follow us: @GPUComputing on Twitter | NVIDIA on Facebook

