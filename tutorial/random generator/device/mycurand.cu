#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#define n 1024

__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			
			
				}

__global__ void fissa_lo_stato(unsigned int*states){
			unsigned int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			states[id]=id;}			
__device__ float rng_uni(unsigned int *state)
{
        // generates uniform ran num, 
        // keeps state val pointed by *state updated 

        unsigned int x = *state;

        x = x ^ (x >> 13);
        x = x ^ (x << 17);
        x = x ^ (x >> 5);

        *state = x;

        return (float) x / 4294967296;
}

__global__ void dyn_kernel(unsigned int *states,float*x)
{
	
        int id = threadIdx.x+blockIdx.x*blockDim.x;
	x[id]=rng_uni(&states[id]);
}

__global__ void genera_random(curandState* stato,float *x){int id=threadIdx.x+ blockIdx.x*blockDim.x;
								
							x[id]=curand(&stato[id]);}
void stampa(float* x,int N){int i=2;
//for(i=0;i<N;i++)
 printf("%1.4f \n",x[i]);}

main(){
float *x;float*dev_x;
int t=100;
x=(float*)malloc(n*n*sizeof(float));
cudaMalloc((float**)&dev_x,n*n*sizeof(float));
unsigned int* states;
cudaMalloc((unsigned int**)&states,n*n*sizeof(unsigned int));
fissa_lo_stato<<<n,n>>>(states);
//ora alloco la variabile di tipo curandState
curandState * stato;
cudaMalloc((void**)&stato,n*n*sizeof(curandState));
setup_random_kernel<<<n,n>>>(stato);
for(int j=0;j<t;j++){
genera_random<<<n,n>>>(stato,dev_x);
dyn_kernel<<<n,n>>>(states,dev_x);
cudaMemcpy(x,dev_x,n*n*sizeof(float),cudaMemcpyDeviceToHost);
stampa(x,n*n);}
cudaFree(dev_x);
cudaFree(stato);
free(x);}
