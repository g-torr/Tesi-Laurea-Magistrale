#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define tau 0.6 // costante nel processo di O-U
#define D 0.1
#define dt 0.1
#define durata size*size*100/(2*D)
#define termalizzazione 0
#define raggio 0.3
#define mu 1
const int blocks=40;
const int threads=1024;
const float size=10.;
const float half_size=size/2.;
struct configurazione{
		float* eta;
		float* x;
		};


void creo_cartelle()					{system( "rm -rf ./forza" );
							mkdir("forza",0700);}
__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che sta utilizza thread richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}
__global__ void inizializza(float *x,float *eta){int id=threadIdx.x+ blockIdx.x*blockDim.x;
					x[id]=0;eta[id]=0.;}

__global__ void evolvi(curandState* stato,float*x, float * eta,float * forza,float* v){float w;int id=threadIdx.x+ blockIdx.x*blockDim.x;	
							w=curand_normal(&stato[id])*sqrt(2.);	
							//w=0;
							eta[id]=eta[id]-(1/tau)*eta[id]*dt+sqrt(D)*w*sqrt(dt)/tau; 
							x[id]=x[id]+eta[id]*dt;
							float delta;forza[id]=0;
						if (x[id]+raggio>half_size){v[id]=(half_size-raggio)/dt-x[id]/dt+eta[id];
									delta=x[id]+raggio-half_size;x[id]=half_size-raggio;forza[id]=delta/(dt*mu);}
						else if (x[id]-raggio<-half_size) {v[id]=(-half_size+raggio)/dt-x[id]/dt+eta[id];
									delta=x[id]-raggio +half_size;x[id]=-half_size+raggio;forza[id]=delta/(dt*mu);} 					else v[id]=eta[id];
							
										}
void 			stampa(float* x,int m)	{int i;for(i=0;i<m;i++) printf("%1.4f %1.4f \n",x[i],x[i+m]);}
void		stampa_file(float*x, int m){FILE*f;char indirizzo [50];int i,j;int N=blocks*threads;
						
						sprintf(indirizzo,"./forza/forza_%d.txt",m);
						f=fopen(indirizzo,"w");
							for(i=0;i<N;i++)fprintf(f,"%f\n",x[i]);
						fclose(f);}
void		stampa_traettoria(float*x,float*v, int m){FILE*f; f=fopen("traettoria.txt","w");int i;
					for(i=0;i<m;i++){fprintf(f,"%f	%f\n",x[i],v[i]);}
							
					fclose(f);}

					

main(){
creo_cartelle();configurazione sistema; int N=blocks*threads;int t; // x è il sistema  dinamico, n è il rumore
cudaEvent_t start,stop; 
float *traettoria; float * storage;float* forza; float* dev_forza;float* velocit;float * v; 
traettoria=(float*)malloc(durata*sizeof(float));
velocit=(float*)malloc(durata*sizeof(float));//sempre per la singola traettoria

cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
storage=(float*)malloc(2*N*sizeof(float));
forza=(float*)malloc(N*sizeof(float));

cudaMalloc((float**)&sistema.x,N*sizeof(float));
cudaMalloc((float**)&sistema.eta,N*sizeof(float));
cudaMalloc((float**)&dev_forza,N*sizeof(float));
cudaMalloc((float**)&v,N*sizeof(float));

inizializza<<<blocks,threads>>>(sistema.x,sistema.eta);
curandState * stato;
cudaMalloc((void**)&stato,N*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(stato);
int i=0;
for(t=0;t<durata;t++){	
			evolvi<<<blocks,threads>>>(stato,sistema.x,sistema.eta,dev_forza,v);
			if(t>=termalizzazione){cudaMemcpy(&traettoria[t-termalizzazione],sistema.x,sizeof(float),cudaMemcpyDeviceToHost);
						cudaMemcpy(&velocit[t-termalizzazione],sistema.eta,sizeof(float),cudaMemcpyDeviceToHost);}
if(t%1000==999)	{cudaMemcpy(forza,dev_forza,N*sizeof(float),cudaMemcpyDeviceToHost);stampa_file(forza,i);i++;}
cudaMemcpy(storage,sistema.x,N*sizeof(float),cudaMemcpyDeviceToHost);
cudaMemcpy(storage+N,v,N*sizeof(float),cudaMemcpyDeviceToHost);}
stampa(storage,N);
free(storage);
stampa_traettoria(traettoria,velocit,durata-termalizzazione);
cudaFree(sistema.x);
cudaFree(sistema.eta);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
//printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}
