#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
//#include <helper_cuda.h>
#include <cuda_runtime_api.h>

const float size=3.;
const float half_size=size/2.;
#define tau 256 // costante nel processo di O-U
#define D 1.
#define dt 0.001
#define mobility 1
const int durata=(int)size*size*400/D;
const float tsalva= durata/20.;
const float raggio=1.; //ricordati che deve essere minore di half_size
const int blocks=256;
const int threads=1024;
const int particelle=blocks*threads;
__constant__ float costanti[2];
struct point{	float x;
		float y;
			};



struct configurazione{	point*	eta;
			point*	r;
			point*	forza;};
void	stampa(configurazione stato, int i){int id;FILE*f;FILE*g;char indirizzo_posizione[50];char indirizzo_forza[50];
					mkdir("posizione",0700);mkdir("forza",0700);
					sprintf(indirizzo_posizione,"./posizione/dati_%d",i);
					sprintf(indirizzo_forza,"./forza/forza_%d",i);
					f=fopen(indirizzo_posizione,"w");g=fopen(indirizzo_forza,"w");
					for(id=0;id<particelle;id++){fprintf(f,"%f		%f \n",stato.r[id].x,stato.r[id].y); 
					fprintf(g,"%f		%f \n",stato.forza[id].x,stato.forza[id].y);}
					fclose(f);fclose(g);}

__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che si utilizza richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}

__device__ point nuova_posizione(float raggio,point old){point nuovo;float d=sqrt((old.x*old.x)+(old.y*old.y));
						nuovo.x=raggio*old.x/d;
						nuovo.y=raggio*old.y/d;
						
						return nuovo;}







__device__ bool is_in(float raggio,point P){float d=sqrt(P.x*P.x+P.y*P.y);
							if (d<raggio) return true;
							else return false;
						}



__global__ void	inizializza(configurazione stato,curandState* gen_random){point temp; int id=threadIdx.x+blockIdx.x*blockDim.x;	
								do {
							temp.x=-half_size+size*curand_uniform(&gen_random[id]);
							temp.y=-half_size+size*curand_uniform(&gen_random[id]);}	
							
								while (!is_in(raggio,temp));
						stato.r[id].x=temp.x;
						stato.r[id].y=temp.y;
						stato.eta[id].x=0;	
						stato.eta[id].y=0;
									 		 		
}
__global__ void	evolvi(configurazione stato,curandState* gen_random){	
							 point r,v,f;
						int id=threadIdx.x+blockIdx.x*blockDim.x; 	
						/*v.x=__ldg(&stato.eta[id].x);	v.y=__ldg(&stato.eta[id].y);
						r.x=__ldg(&stato.r[id].x);	r.x=__ldg(&stato.r[id].x);*/
						v=stato.eta[id];	r=stato.r[id];/*
						v.x=v.x-(1.f/tau)*v.x*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau;   
						v.y=v.y-(1.f/tau)*v.y*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau; 	*/
						v.x=v.x-v.x*costanti[0]+curand_normal(&gen_random[id])*costanti[1];   	
						v.y=v.y-v.y*costanti[0]+curand_normal(&gen_random[id])*costanti[1]; 	
						r.x=r.x+v.x*dt;
						r.y=r.y+v.y*dt;	
						f.x=0.f; f.y=0.f;
						
					bool in=is_in(raggio,r);
					if (!in){
										f=r;
										r=nuova_posizione(raggio,r);
										
										f.x=(f.x-r.x)/(dt*mobility);
										f.y=(f.y-r.y)/(dt*mobility); 
							}
					

																	
			if(r.x>half_size)	r.x=r.x-size;
				else if (r.x<-half_size) r.x=r.x+size;
			if	(r.y>half_size)	r.y=r.y-size;
				else if (r.y<-half_size)	r.y=r.y+size;
			
						stato.eta[id]=v;	stato.forza[id]=f;	stato.r[id]=r;

						
//if(id==5)printf(" la particella 5 si trova in %f \n",stato.r[id].x);
}


main(){
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
int numero_passi=(int)durata/dt;
int passi_salvataggio =(int)(tsalva/dt);
float costanti_host[2];
costanti_host[0]=dt/tau;
costanti_host[1]=(float)sqrt(2*D* dt)/tau;
cudaMemcpyToSymbol(costanti,costanti_host,2*sizeof(float));
if(raggio>half_size){printf("attenzione!! raggio deve essere maggiore di half size");return 1;}

cudaStream_t stream0,stream1,stream2,stream3;
cudaStreamCreate(&stream0);cudaStreamCreate(&stream1);cudaStreamCreate(&stream2);cudaStreamCreate(&stream3);//****creo gli stream
configurazione dev_stato_0,dev_stato_1,dev_stato_2,dev_stato_3,stato;//alloco lo stato del sistema sul device e sull'host tramite cudaHostAlloc (pinned memory)

cudaMalloc((point**)&dev_stato_0.eta,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_0.r,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_0.forza,(particelle/4)*sizeof(point));

cudaMalloc((point**)&dev_stato_1.eta,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_1.r,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_1.forza,(particelle/4)*sizeof(point));


cudaMalloc((point**)&dev_stato_2.eta,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_2.r,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_2.forza,(particelle/4)*sizeof(point));


cudaMalloc((point**)&dev_stato_3.eta,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_3.r,(particelle/4)*sizeof(point));
cudaMalloc((point**)&dev_stato_3.forza,(particelle/4)*sizeof(point));

cudaHostAlloc((point**)&stato.eta,particelle*sizeof(point),cudaHostAllocDefault);
cudaHostAlloc((point**)&stato.r,particelle*sizeof(point),cudaHostAllocDefault);
cudaHostAlloc((point**)&stato.forza,particelle*sizeof(point),cudaHostAllocDefault);

curandState * generatori_random;//alloco il generatore dei numeri random sul device
if (numero_passi> pow(2,67)) printf("ATTENZIONE! ricontrollare il generatore di numeri casuali"); 
cudaMalloc((void**)&generatori_random,particelle*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(generatori_random);

inizializza<<<blocks/4,threads,0,stream0>>>(dev_stato_0,generatori_random);
inizializza<<<blocks/4,threads,0,stream1>>>(dev_stato_1,generatori_random);
inizializza<<<blocks/4,threads,0,stream2>>>(dev_stato_2,generatori_random);
inizializza<<<blocks/4,threads,0,stream3>>>(dev_stato_3,generatori_random);

int i=0;
for(int t=0;t<numero_passi;t++){
//for(t=0;t<10;t++){
			evolvi<<<blocks/4,threads,0,stream0>>>(dev_stato_0,generatori_random);
			evolvi<<<blocks/4,threads,0,stream1>>>(dev_stato_1,generatori_random);
			evolvi<<<blocks/4,threads,0,stream2>>>(dev_stato_2,generatori_random);
			evolvi<<<blocks/4,threads,0,stream3>>>(dev_stato_3,generatori_random);
	
	   if((t% passi_salvataggio==0)&&(t>0)&&(t>10*passi_salvataggio)){printf("siamo a %d/9 \n",i);
			cudaMemcpyAsync(stato.r,dev_stato_0.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream0);
			cudaMemcpyAsync(stato.r + (particelle/4),dev_stato_1.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream1);
			cudaMemcpyAsync(stato.r + (particelle/2),dev_stato_2.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream2);
			cudaMemcpyAsync(stato.r + (3*particelle/4),dev_stato_3.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream3);
	
			
			cudaMemcpyAsync(stato.forza,dev_stato_0.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream0); 			
			cudaMemcpyAsync(stato.forza + (particelle/4),dev_stato_1.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream1);
			cudaMemcpyAsync(stato.forza + (particelle/2),dev_stato_2.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream2);
			cudaMemcpyAsync(stato.forza + (3*particelle/4),dev_stato_3.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream3);
	  		cudaDeviceSynchronize();
			stampa(stato,i);i++;}

			}
 

cudaMemcpyAsync(stato.r,dev_stato_0.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream0);
cudaMemcpyAsync(stato.r + (particelle/4),dev_stato_1.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream1);
cudaMemcpyAsync(stato.r + (particelle/2),dev_stato_2.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream2);
cudaMemcpyAsync(stato.r + (3*particelle/4),dev_stato_3.r,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream3);
	
			
cudaMemcpyAsync(stato.forza,dev_stato_0.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost,stream0); 			
cudaMemcpyAsync(stato.forza + (particelle/4),dev_stato_1.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream1);
cudaMemcpyAsync(stato.forza + (particelle/2),dev_stato_2.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream2);
cudaMemcpyAsync(stato.forza + (3*particelle/4),dev_stato_3.forza,(particelle/4)*sizeof(point),cudaMemcpyDeviceToHost, stream3);
 
cudaDeviceSynchronize();
stampa(stato,i);
cudaFree(dev_stato_0.r);
cudaFree(dev_stato_0.eta);
cudaFree(dev_stato_0.forza);
cudaFree(dev_stato_1.r);
cudaFree(dev_stato_1.eta);
cudaFree(dev_stato_1.forza);
cudaFree(dev_stato_2.r);
cudaFree(dev_stato_2.eta);
cudaFree(dev_stato_2.forza);
cudaFree(dev_stato_3.r);
cudaFree(dev_stato_3.eta);
cudaFree(dev_stato_3.forza);

cudaFreeHost(stato.r);
cudaFreeHost(stato.forza);
cudaFreeHost(stato.eta);
cudaStreamSynchronize(stream0);		
cudaStreamSynchronize(stream1);
cudaStreamSynchronize(stream2);		
cudaStreamSynchronize(stream3);
 cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
cudaStreamDestroy(stream0);cudaStreamDestroy(stream1);cudaStreamDestroy(stream2);cudaStreamDestroy(stream3);
}



