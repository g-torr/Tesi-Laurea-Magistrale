#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#define N 20*50000
#define imin(a,b) (a<b?a:b)
const int threadsperblock= 256;
const int numberofblocks=imin(32,(N+ threadsperblock-1)/threadsperblock);

__global__ void somma(float *a, float *somme){int i= threadIdx.x+ blockIdx.x *blockDim.x;
	__shared__ float parziale[threadsperblock]; //si sfrutta la proprietà della shared memory, infatti vengono creati numberofblocks copie della variabile  parziale, cosicchè l'indice di parziale gira solo sui treads di un singolo blocco, 
		int cacheindex=threadIdx.x;
	float temp=0;
	while(i<N){temp=temp + *(a+i);
			i=i+ blockDim.x*gridDim.x;}
	parziale[cacheindex]=temp; // in questo passaggio l'indice di parziale gira solo sui treads di un singolo blocco, ma concretamente noi abbiamo   numberofblocks  di questi threads che mi ricoprono tutti gli N elementi (in realtà l'ultimo blocco potrebbe rimanere incompleto se N non è un multiplo di threadsperblock)

	__syncthreads();
	int k=blockDim.x/2;  //riduzione, funziona solo per un threadsperblock pari a una potenza di 2
	while(k!=0){
	if(cacheindex<k) parziale[cacheindex]+=parziale[cacheindex+k];
			__syncthreads();	k=k/2;}
	if (cacheindex==0) somme[blockIdx.x] =parziale[0];
						}


void stampa(float *a,int n){int i;
		for(i=0;i<n;i++)printf("%f\n",a[i]);}

__global__ void somma2(float *a,float media,float* somme2){
			int i=threadIdx.x+blockIdx.x *blockDim.x ;
			float temp=0;
			while(i<N){	temp=temp+a[i]*a[i];	 i+=blockDim.x*gridDim.x;}
		__shared__ float parziale[threadsperblock];
			parziale[threadIdx.x]=temp;
			__syncthreads();
			int j= blockDim.x/2;
			while (j!=0){
					if(threadIdx.x< j) parziale[threadIdx.x]+= parziale[threadIdx.x+j];
					__syncthreads(); j=j/2; }
			if(threadIdx.x==0)	somme2[blockIdx.x]= parziale[0];
						}

main(){ float* a; float* dev_a;
float *somme; float* dev_somme;
int i;
srand(time(NULL));
a=(float*)malloc(N*sizeof(float));
somme=(float*)malloc(numberofblocks*sizeof(float));
cudaMalloc((void**)&dev_somme,numberofblocks*sizeof(float));
cudaMalloc((void**)&dev_a, N*sizeof(float));
for(i=0;i<N;i++) *(a+i)=(float)rand()/RAND_MAX;
cudaMemcpy(dev_a,a,N*sizeof(float),cudaMemcpyHostToDevice);
cudaMemcpy(dev_somme,somme,numberofblocks*sizeof(float),cudaMemcpyHostToDevice);
somma<<<numberofblocks,threadsperblock>>>(dev_a,dev_somme);
cudaMemcpy(somme,dev_somme,numberofblocks*sizeof(float),cudaMemcpyDeviceToHost);
//stampa(somme,numberofblocks);
float media=0;
for(i=0;i<numberofblocks;i++){
media=media+ *(somme+i);
}
int n=N;
media= media /n;
printf("la media è %f \n",media);

// ora calcolo la varianza
float * somme2; float * dev_somme2;
somme2=(float*)malloc(numberofblocks*sizeof(float));
cudaMalloc((void**)&dev_somme2,numberofblocks*sizeof(float));
somma2<<<numberofblocks,threadsperblock>>>(dev_a,media,dev_somme2);
cudaMemcpy(somme2,dev_somme2,numberofblocks*sizeof(float),cudaMemcpyDeviceToHost);
float varianza=0;
stampa (somme2,numberofblocks);
for(i=0;i<numberofblocks;i++){varianza+=somme2[i];}		
			varianza=varianza/n- media*media;
	printf("varianza = %f\n",varianza);
return 0;
}
