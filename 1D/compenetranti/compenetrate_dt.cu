#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define tau 0.06 // costante nel processo di O-U
#define D 0.1
//#define durata 9000
//#define termalizzazione 70
#define raggio 0.3
#define mu 1
const int blocks=10;
const int threads=1024;
const float durata=4./(2.*D)*10.;
const float tsalva= durata/10.;
struct configurazione{
		float* eta;
		float* x;
		};
__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}
__global__ void inizializza(float *x,float *eta){int id=threadIdx.x+ blockIdx.x*blockDim.x;
					x[id]=0;eta[id]=0;}

__global__ void evolvi(curandState* stato,float*x, float * eta,float * forza,float dt){float w;int id=threadIdx.x+ blockIdx.x*blockDim.x;	
							w=curand_normal(&stato[id])*sqrt(2.);	
							eta[id]=eta[id]-(1/tau)*eta[id]*dt+sqrt(D)*w*sqrt(dt)/tau; 
							x[id]=x[id]+eta[id]*dt;
							float delta;forza[id]=0;
							if (x[id]+raggio>1){delta=x[id]+raggio-1;x[id]=1-raggio;forza[id]=delta/(dt*mu);}
						else if (x[id]-raggio<-1) {delta=x[id]-raggio +1;x[id]=-1+raggio;forza[id]=delta/(dt*mu);} 
							
										}
void 			stampa(float* x,float dt)	{FILE*f;char indirizzo [50];int i,j;int N=blocks*threads;
						sprintf(indirizzo,"./dati_%f",dt);
						f=fopen(indirizzo,"w");	
						for(i=0;i<N;i++) fprintf(f,"%1.4f \n",x[i]);
							fclose(f);}

void		stampa_file(float*x, int m,float dt){FILE*f;char indirizzo [50];int i,j;int N=blocks*threads;
						
						sprintf(indirizzo,"forza%f",dt);
            mkdir(indirizzo,0700);
            char indirizzo2[50]; sprintf(indirizzo2,"./forza%f/forza_%d",dt,m);
						f=fopen(indirizzo2,"w");
							for(i=0;i<N;i++)fprintf(f,"%f\n",x[i]);
						fclose(f);}
					

int main(void){
configurazione sistema; int numero_passi;
int N=blocks*threads;int t; // x è il sistema  dinamico, n è il rumore
cudaEvent_t start,stop; 
float *traettoria; float * storage;float* forza; float* dev_forza; 
traettoria=(float*)malloc(durata*sizeof(float));
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
storage=(float*)malloc(N*sizeof(float));
forza=(float*)malloc(N*sizeof(float));

cudaMalloc((float**)&sistema.x,N*sizeof(float));
cudaMalloc((float**)&sistema.eta,N*sizeof(float));
cudaMalloc((float**)&dev_forza,N*sizeof(float));
inizializza<<<blocks,threads>>>(sistema.x,sistema.eta);
curandState * stato;
cudaMalloc((void**)&stato,N*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(stato);

double dt;int j;dt=2.;


for(j=1;j<17;j++){

  dt=dt/(2.);
  printf("%f\n",dt);
  
  inizializza<<<blocks,threads>>>(sistema.x,sistema.eta);
  
  int numero_passi=(int)(durata/dt);
  int passi_salvataggio =(int)(tsalva/dt);
  int i=0;
  
  for(t=0;t<numero_passi;t++){
    evolvi<<<blocks,threads>>>(stato,sistema.x,sistema.eta,dev_forza,dt);
//			if(t>=termalizzazione)cudaMemcpy(&traettoria[t-termalizzazione],sistema.x,sizeof(float),cudaMemcpyDeviceToHost);}
    if((t% passi_salvataggio==0)&&(t>0)){
      cudaMemcpy(forza,dev_forza,N*sizeof(float),cudaMemcpyDeviceToHost);
      stampa_file(forza,i,dt);i++;
    }
  }

  cudaMemcpy(storage,sistema.x,N*sizeof(float),cudaMemcpyDeviceToHost);
  stampa(storage,dt);
 }


free(storage);
cudaFree(sistema.x);
cudaFree(sistema.eta);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
//printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}
