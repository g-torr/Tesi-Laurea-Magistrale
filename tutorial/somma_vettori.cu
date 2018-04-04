#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#define N 4 
__global__ void somma(int* a, int* b, int* c,int *j){
int i=blockIdx.x; *j=i;
if (i<N) c[i]=a[i]+b[i];
}


int main(){
printf("è l'esercizio di pag 41 del libro");
int* a;int*b;int *c;
int *dev_a; int* dev_b; int* dev_c;
cudaMalloc((void**)&dev_a, N*sizeof(int));
cudaMalloc((int**)&dev_b, N*sizeof(int));
cudaMalloc((void**)&dev_c,N*sizeof(int));

a=(int*)malloc(N*sizeof(int));
b=(int*)malloc(N*sizeof(int));
c=(int*)malloc(N*sizeof(int));
int i;
for(i=0;i<N;i++){
//printf("acquisisci a, b \n");
*(a+i)=i;
*(b+i)=i*i;}
cudaMemcpy(dev_a,a,N*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dev_b,b,N*sizeof(int),cudaMemcpyHostToDevice);
int j; int* dev_j;                     	//j è il  blockIdx.x che copierò dal device all'host
cudaMalloc((void**)&dev_j,sizeof(int));
somma<<<N,1>>>(dev_a,dev_b,dev_c,dev_j);
cudaMemcpy(c,dev_c,N*sizeof(int),cudaMemcpyDeviceToHost);
cudaMemcpy(&j,dev_j,sizeof(int),cudaMemcpyDeviceToHost);
printf("il blockIdx.x =  %d",j);

for(i=0;i<N;i++){
printf("il tuo numero è %d \n",*(c+i));}

return 0;}
