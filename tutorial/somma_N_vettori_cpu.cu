#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

/*__global__ void somma(float *a,float *b, float *c, int N){
		int i=threadIdx.x+ blockIdx.x * blockDim.x;
			while(i<N){*(c+i)=*(a+i)+*(b+i);
				i=i+(gridDim.x*blockDim.x);}}*/
main(){ 
clock_t t1 = clock();
srand(time(NULL));
	int N,i;
printf("inserire il numero di elementi del vettore");
scanf("%d",&N);
float *a; float *b; float *c;

a=(float*)malloc(N*sizeof(float));
b=(float*)malloc(N*sizeof(float));
c=(float*)malloc(N*sizeof(float));
for(i=0;i<N;i++){*(a+i)=(float)rand()/RAND_MAX;
		*(b+i)=(float)rand()/RAND_MAX;}
float *dev_a; float * dev_b; float *dev_c;
/*cudaMalloc((void**)&dev_a, N*sizeof(float));
cudaMalloc((void**)&dev_b, N*sizeof(float));
cudaMalloc((void**)&dev_c, N*sizeof(float));
cudaMemcpy(dev_a,a,N*sizeof(float),cudaMemcpyHostToDevice);
cudaMemcpy(dev_b,b,N*sizeof(float),cudaMemcpyHostToDevice);
somma<<<12,12>>>(dev_a,dev_b,dev_c,N);
cudaMemcpy(c,dev_c,N*sizeof(float),cudaMemcpyDeviceToHost);
*/
for(i=0;i<N;i++){*(c+i)=*(a+i)+*(b+i);}

//for(i=0;i<N;i++){
//printf("%f + %f = %f \n",*(a+i),*(b+i),*(c+i));}

 clock_t t2 = clock();
double time_sec = 
       (double)(t2-t1)/(double)(CLOCKS_PER_SEC); 
 
    printf("Time (sec): %lf\n",time_sec); 
return 0;
}

