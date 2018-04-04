#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#define n 40 //descrivo n*n particelle
#define k 0.6 // costante nel processo di O-U
#define D 0.3
#define dt 0.1
#define durata 3000
#define termalizzazione 70

__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}
__global__ void inizializza(float *x){int id=threadIdx.x+ blockIdx.x*blockDim.x;
					x[id]=0;}

__global__ void evolvi(curandState* stato,float* x){float w;int id=threadIdx.x+ blockIdx.x*blockDim.x;	
							w=curand_normal(&stato[id]);	
							x[id]=x[id]-k*x[id]*dt+D*w*sqrt(dt); 
							
}
void 			stampa(float* x,int m)	{int i;for(i=0;i<m;i++) printf("%1.4f \n",x[i]);}
void		stampa_file(float*x, int m){FILE*f; f=fopen("traettoria.txt","w");int i;
					for(i=0;i<m;i++){fprintf(f,"%f\n",x[i]);}
							
					fclose(f);}

main(){float * x;float* dev_x;int N=n*n;int t;
cudaEvent_t start,stop;
float *traettoria;
traettoria=(float*)malloc(durata*sizeof(float));
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
x=(float*)malloc(N*sizeof(float));
//cudaMemset((float*)&dev_x,0,N*sizeof(float));
cudaMalloc((float**)&dev_x,N*sizeof(float));
cudaMalloc((float**)&dev_x,N*sizeof(float));
inizializza<<<n,n>>>(dev_x);
curandState * stato;
cudaMalloc((void**)&stato,N*N*sizeof(curandState));
setup_random_kernel<<<n,n>>>(stato);
for(t=0;t<durata;t++){	evolvi<<<n,n>>>(stato,dev_x);
			if(t>=termalizzazione)cudaMemcpy(&traettoria[t-termalizzazione],dev_x,sizeof(float),cudaMemcpyDeviceToHost);}//sto dicendo che aspetto t=100 per la termalizzazione
cudaMemcpy(x,dev_x,N*sizeof(float),cudaMemcpyDeviceToHost);

stampa(x,N);
stampa_file(traettoria,durata-termalizzazione);
free(x);
cudaFree(dev_x);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
//printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}
