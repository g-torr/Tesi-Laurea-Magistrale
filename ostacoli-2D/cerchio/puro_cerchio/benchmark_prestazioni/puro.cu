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
#define tau 0.5 // costante nel processo di O-U
#define D 0.1
#define dt 0.01
#define mobility 1
const int durata=(int)size*size*10/D;
int numero_passi=(int)(durata/dt);
int passi_salvataggio =(int)(durata/(10.*dt));

const float raggio=1.; //ricordati che deve essere minore di half_size
const int blocks=100;
const int threads=1024;
const int particelle=blocks*threads;
struct point{	float x;
		float y;
			};

__constant__ float costanti[2];

struct configurazione{	point*	eta;
			point*	r;
			point*	forza;
			bool* 	inside;};
void creo_cartelle()				{system( "rm -rf ./posizione" );system( "rm -rf ./forza" );system( "rm -rf ./velocity" );
							mkdir("posizione",0700);mkdir("forza",0700);mkdir("velocity",0700);}

void	stampa(configurazione stato, int i){int id;FILE*f;FILE*g;FILE*h;char indirizzo_posizione[50];char indirizzo_forza[50];char indirizzo_velocity[50];
					sprintf(indirizzo_posizione,"./posizione/dati_%d",i);
					sprintf(indirizzo_forza,"./forza/forza_%d",i);
					sprintf(indirizzo_velocity,"./velocity/velocity_%d",i);
	
					f=fopen(indirizzo_posizione,"w");g=fopen(indirizzo_forza,"w");h=fopen(indirizzo_velocity,"w");
		for(id=0;id<particelle;id++){fprintf(f,"%f		%f \n",stato.r[id].x,stato.r[id].y); //printf("sto scrivendo %d\n",id);
					fprintf(g,"%f		%f \n",stato.forza[id].x,stato.forza[id].y);
					fprintf(h,"%f		%f \n",stato.eta[id].x -stato.forza[id].x,stato.eta[id].y-stato.forza[id].y);
						}
					fclose(f);fclose(g);fclose(h);}

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
						stato.r[id].x=(float)temp.x;
						stato.r[id].y=(float)temp.y;
						stato.eta[id].x=0.f;	
						stato.eta[id].y=0.f;
									 		 		
}
__global__ void	evolvi(configurazione stato,curandState* gen_random){	
			point r,v,f;
			int id=threadIdx.x+blockIdx.x*blockDim.x; 	
			v=stato.eta[id];	r=stato.r[id];
			v.x=v.x-v.x*costanti[0]+curand_normal(&gen_random[id])*costanti[1];   	
			v.y=v.y-v.y*costanti[0]+curand_normal(&gen_random[id])*costanti[1]; 	
			r.x=r.x+v.x*dt;
			r.y=r.y+v.y*dt;	
			bool in =is_in(raggio,stato.r[id]);
					if (!in) stato.inside[id]=true;
					else {stato.inside[id]=false;f.x=0.f; f.y=0.f;stato.forza[id]=f;	

						/*if(r.x>half_size)	r.x=r.x-size;
						else if (r.x<-half_size) r.x=r.x+size;
						if	(r.y>half_size)	r.y=r.y-size;
						else if (r.y<-half_size)	r.y=r.y+size;*/
						}
												
			stato.eta[id]=v;		stato.r[id]=r;		
				/*	if (in){
										f=r;
										r=segmento_vicino(r);
										
										f.x=(f.x-r.x)/(dt*mobility);
										f.y=(f.y-r.y)/(dt*mobility); 
							}
					

																}	
		else {	if(r.x>half_size)	r.x=r.x-size;
				else if (r.x<-half_size) r.x=r.x+size;
			if	(r.y>half_size)	r.y=r.y-size;
				else if (r.y<-half_size)	r.y=r.y+size;
			}
						stato.eta[id]=v;	stato.forza[id]=f;	stato.r[id]=r;

						*/
}

__global__ void correct(configurazione stato,curandState* gen_random,int * ids,int count){
				point r,f;
				int k=threadIdx.x+blockIdx.x*blockDim.x;
				while(k<count){
				
				int id=ids[k]; 	
				r=stato.r[id];
				f=r;
				r=nuova_posizione(raggio,r);
				f.x=(f.x-r.x)/(dt*mobility);
				f.y=(f.y-r.y)/(dt*mobility); 
				stato.forza[id]=f;	stato.r[id]=r;		
				k=k+(gridDim.x*blockDim.x);}
					}

			


main(){
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
creo_cartelle();
int i,t;

float costanti_host[2];
costanti_host[0]=dt/tau;
costanti_host[1]=(float)sqrt(2*D* dt)/tau;
cudaMemcpyToSymbol(costanti,costanti_host,2*sizeof(float));

if(raggio>half_size){printf("attenzione!! raggio deve essere maggiore di half size");return 1;}

configurazione dev_stato,stato;//alloco lo stato del sistema
cudaMalloc((point**)&dev_stato.eta,particelle*sizeof(point));
cudaMalloc((point**)&dev_stato.r,particelle*sizeof(point));
cudaMalloc((point**)&dev_stato.forza,particelle*sizeof(point));
cudaMalloc((bool**)&dev_stato.inside,particelle*sizeof(bool));

curandState * generatori_random;//alloco il generatore dei numeri random sul device
if (numero_passi> pow(2,67)) printf("ATTENZIONE! ricontrollare il generatore di numeri casuali"); 
cudaMalloc((void**)&generatori_random,particelle*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(generatori_random);



inizializza<<<blocks,threads>>>(dev_stato,generatori_random);
stato.eta=(point*)malloc(particelle*sizeof(point));
stato.r  =(point*)malloc(particelle*sizeof(point));
stato.forza=(point*)calloc(particelle,sizeof(point));
stato.inside=(bool*)malloc(particelle*sizeof(bool));

i=0;int count;
int*  ids;int*  dev_ids;
cudaMalloc((int**)&dev_ids,particelle*sizeof(int));
ids=(int*)malloc(particelle*sizeof(int));
for(t=0;t<numero_passi;t++){
			evolvi<<<blocks,threads>>>(dev_stato,generatori_random);
			cudaMemcpy(stato.inside,dev_stato.inside,particelle*sizeof(bool),cudaMemcpyDeviceToHost);		
			count=0;
			for(int j=0;j<particelle;j++){
					if(stato.inside[j]==true){ ids[count]=j;count++;}
							}
		//	printf("numero di pallette che entrano %d \n", count);
			
			cudaMemcpy(dev_ids,ids,count*sizeof(int),cudaMemcpyHostToDevice);
			correct<<<(count+1023)/1024,1024>>>(dev_stato,generatori_random,dev_ids,count);					
							
			   if((t% passi_salvataggio==0)&&(t>0)){printf("siamo a %d/9 \n",i);
			cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
			cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
			cudaMemcpy(stato.eta,dev_stato.eta,particelle*sizeof(point),cudaMemcpyDeviceToHost);
				stampa(stato,i);i++;
	
}
//cudaDeviceSynchronize();
			}

 

cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
cudaMemcpy(stato.eta,dev_stato.eta,particelle*sizeof(point),cudaMemcpyDeviceToHost);
stampa(stato,i);
cudaFree(dev_stato.r);
cudaFree(dev_stato.eta);
cudaFree(dev_stato.forza);
free(stato.r);
free(stato.forza);
free(stato.eta);
 cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}


 

