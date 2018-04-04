#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define tau 1. // costante nel processo di O-U
#define k 0.001	//convessitù parabola potenziale armonico
#define D 0.1
#define size 3.
#define half_size size/2.
const int termalizzazione= half_size*half_size/(2.*D);
const int blocks=600;
const int threads=1024;
const float durata=size*size/(2.*D)*4.;
const float tsalva= durata/60.;
struct configurazione{
		double* x;
		double* eta;
		};
__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}
__global__ void inizializza(double *x,double *eta){int id=threadIdx.x+ blockIdx.x*blockDim.x;
					x[id]=0.75;eta[id]=0.5;}

__global__ void evolvi(curandState* stato,double*x, double * eta,double* forza,double dt){double w;int id=threadIdx.x+ blockIdx.x*blockDim.x;	
							w=curand_normal(&stato[id])*sqrt(2.0);	
							eta[id]=eta[id]-(1/tau)*eta[id]*dt+sqrt(D)*w*sqrt(dt)/tau; //Ornstein Oulenbeck
              
	forza[id]=12*k*pow(x[id],-13);          
	x[id] = x[id] +forza[id]*dt + eta[id]*dt;// processo dinamico 
	  
         
					
          if (x[id]>=half_size) x[id] = x[id]-size;
					else if(x[id]< -half_size) x[id] = x[id]+size;		      
}

void 			stampa(double* x,double dt)	{FILE*f;char indirizzo [50];int i,j;int N=blocks*threads;
						sprintf(indirizzo,"./dati_%f",dt);
						f=fopen(indirizzo,"w");	
						for(i=0;i<N;i++) fprintf(f,"%1.4f \n",x[i]);
							fclose(f);}

void		stampa_file(double*x, int m,double dt){FILE*f;char indirizzo [50];int i,j;int N=blocks*threads;
						
						sprintf(indirizzo,"forza%f",dt);
            mkdir(indirizzo,0700);
            char indirizzo2[50]; sprintf(indirizzo2,"./forza%f/forza_%d",dt,m);
						f=fopen(indirizzo2,"w");
							for(i=0;i<N;i++)fprintf(f,"%f\n",x[i]);
						fclose(f);}

main(){
  configurazione sistema; int N=blocks*threads;int t; // x è il sistema  dinamico, n è il rumore
  cudaEvent_t start,stop; 
  double *traettoria; double * storage; double* forza; double* dev_forza; 
  traettoria=(double*)malloc(durata*sizeof(double));
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
  storage=(double*)malloc(N*sizeof(double));
	forza=(double*)malloc(N*sizeof(double));
cudaMalloc((double**)&sistema.x,N*sizeof(double));
cudaMalloc((double**)&sistema.eta,N*sizeof(double));
cudaMalloc((double**)&dev_forza,N*sizeof(double));

curandState * stato;
cudaMalloc((void**)&stato,N*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(stato);

double dt;int j;dt=0.2;
	int numero_passi,passi_salvataggio,passi_termalizzazione;

for(j=0;j<10;j++){
	dt=dt/2;
  printf("%f\n",dt);
 inizializza<<<blocks,threads>>>(sistema.x,sistema.eta);
  
   numero_passi=(int)(durata/dt);
   passi_salvataggio =(int)(tsalva/dt);
  passi_termalizzazione=(int)(termalizzazione/dt);
  int i=0;
  
   for(t=0;t<numero_passi;t++){
    evolvi<<<blocks,threads>>>(stato,sistema.x,sistema.eta,dev_forza,dt);
//			if(t>=termalizzazione)cudaMemcpy(&traettoria[t-termalizzazione],sistema.x,sizeof(float),cudaMemcpyDeviceToHost);}
    if((t% passi_salvataggio==0)&&(t>passi_termalizzazione)){
      cudaMemcpy(forza,dev_forza,N*sizeof(double),cudaMemcpyDeviceToHost);
      stampa_file(forza,i,dt);i++;
    }
  }

  cudaMemcpy(storage,sistema.x,N*sizeof(double),cudaMemcpyDeviceToHost);
  stampa(storage,dt);
 }


free(storage);
cudaFree(sistema.x);
cudaFree(sistema.eta);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}
