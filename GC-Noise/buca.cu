#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#define tau 0.6 // costante nel processo di O-U
#define k 2.0	//convessitù parabola potenziale armonico
#define D 0.1
#define dt 0.001
#define durata 3000
#define termalizzazione 70

const int blocks=40;
const int threads=1024;
struct configurazione{
		float* x;
		float* eta;
		};
__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}
__global__ void inizializza(float *x,float *eta){int id=threadIdx.x+ blockIdx.x*blockDim.x;
					x[id]=0;eta[id]=0;}

__global__ void evolvi(curandState* stato,float*x, float * eta){float w;int id=threadIdx.x+ blockIdx.x*blockDim.x;	
							w=curand_normal(&stato[id])*sqrt(2.);	
							x[id]=x[id]-(1/tau)*x[id]*dt+sqrt(D)*w*sqrt(dt)/tau; 
							eta[id]=eta[id]-k*eta[id]*dt+x[id]*dt;
							
}
void 			stampa(float* x,int m)	{int i;for(i=0;i<m;i++) printf("%1.4f \n",x[i]);}
void		stampa_file(float*x, int m){FILE*f; f=fopen("traettoria.txt","w");int i;
					for(i=0;i<m;i++){fprintf(f,"%f\n",x[i]);}
							
					fclose(f);}

main(){configurazione sistema; int N=blocks*threads;int t; // x è il sistema  dinamico, n è il rumore
cudaEvent_t start,stop; 
float *traettoria; float * storage; 
traettoria=(float*)malloc(durata*sizeof(float));
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
storage=(float*)malloc(N*sizeof(float));

cudaMalloc((float**)&sistema.x,N*sizeof(float));
cudaMalloc((float**)&sistema.eta,N*sizeof(float));
inizializza<<<blocks,threads>>>(sistema.x,sistema.eta);
curandState * stato;
cudaMalloc((void**)&stato,N*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(stato);
for(t=0;t<durata;t++){	evolvi<<<blocks,threads>>>(stato,sistema.x,sistema.eta);
			if(t>=termalizzazione)cudaMemcpy(&traettoria[t-termalizzazione],sistema.eta,sizeof(float),cudaMemcpyDeviceToHost);}
cudaMemcpy(storage,sistema.eta,N*sizeof(float),cudaMemcpyDeviceToHost);

stampa(storage,N);
stampa_file(traettoria,durata-termalizzazione);
free(storage);
cudaFree(sistema.x);
cudaFree(sistema.eta);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
//printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}
