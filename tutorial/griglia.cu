#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#define imin(a,b) (a<b?a:b)
const int DIM=4;
const int threadsperblock=2;
const int durata=10;
#define rate 0.03f
/*
texture <float,2> text_Told;
texture <float,2> text_Tnew;
texture <float,2> text_Tfix;*/
struct dati{
		float* old;
		float* upd;
		};


void	stampa(float*T, int n){int i,j;for(i=0;i<n;i++){for(j=0;j<n;j++){printf("%f	",T[i+(j*n)]);}
							printf("\n");}}
	
void	stampa_file(float*T, int n){FILE*f; f=fopen("output.txt","w");int i,j;
							fprintf(f," Vi sono n+1 matrici, la matrice 0 è la condizione iniziale, mentre l'elemento 0 del vettore di matrici updated rappresenta il primo timestep\n copio la matrice %d\n",0);
					for(i=0;i<n;i++){for(j=0;j<n;j++){fprintf(f,"%f	",T[i+(j*n)]);}
							fprintf(f,"\n");}
					fclose(f);}
void	stampa_file2(float ** storage){FILE*f; f=fopen("output.txt","a");int i,n,j;
for(n=0;n<durata;n++){
fprintf(f," copio la matrice %d\n",n+1);
					for(i=0;i<DIM;i++){for(j=0;j<DIM;j++){fprintf(f,"%f	",storage[n][i+(j*DIM)]);}
							fprintf(f,"\n");}}
fclose(f);}


void	 constrains(float* T){	int i;
				 float p;
				for(i=0;i<(DIM*DIM);i++){p=(float)rand()/RAND_MAX;
									if(p>0.7) T[i]=1;
									else T[i]=0;	
					}
				}


void     inizializza(float* T_fix,float *T)	{ int i;
                                float p;
                                for(i=0;i<(DIM*DIM);i++)	{if(T_fix[i]==0)	{p=(float)rand()/RAND_MAX;T[i]=p;}
                                                                        	else T[i]=1;
                                        		}			}
                               		    	 

__global__ void apply_constrains(float *fix, float *old){int x,y,offset;
						 x=threadIdx.x+blockIdx.x*blockDim.x;
						 y=threadIdx.y+blockIdx.y*blockDim.y;
						offset= x+y*DIM;
						if(fix[offset]!=0){old[offset]=1;} 
							}
							
__global__ void execute(float* old, float* upd){int x,y,offset;
							x=threadIdx.x+blockIdx.x*blockDim.x;
							y=threadIdx.y+blockIdx.y*blockDim.y;
							offset=x+y*DIM;
							int u,b,l,r;//interazione a primi vicini:up,bottom,left,right
							//periodic boundary condition
							if(x==0) 	{u=offset+(DIM-1);b=offset+1;}
							else if(x==DIM-1){b=offset-(DIM-1);u=offset-1;}
							else		{u=offset-1;b=offset+1;}
							if (y==0)	{l=offset+(DIM-1)*DIM;r=offset+DIM;}
							else if(y==DIM-1){r=offset-(DIM-1)*DIM;l=offset-DIM;}
							else		{r=offset+DIM;l=offset-DIM;}
							upd[offset]=old[offset]+(rate*(old[u]+old[b]+old[l]+old[r]-(4*old[offset])));


								}


void	time_step(dati* T,float * T_fissata)	{
			//	for(n=0;n<durata;n++){
				dim3 blocks(DIM/threadsperblock,DIM/threadsperblock);
				dim3 threads (threadsperblock,threadsperblock);
					
				apply_constrains<<<blocks,threads>>>(T_fissata,T->old);
				execute<<<blocks,threads>>>(T->old,T->upd);
				T->old=T->upd;				
			}
				




main(){srand(5);
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
float* T_fissata;float* T_iniziale;
T_fissata=(float*)malloc(DIM*DIM*sizeof(float));
constrains(T_fissata);	//pongo i vincoli su alcune celle scelte casualmente con probabilità = 0.3

T_iniziale=(float*)malloc(DIM*DIM*sizeof(float));
inizializza(T_fissata,T_iniziale);
dati T_dev;
//T_dev=(dati*)malloc(durata*sizeof(dati));


int n=0;			//definisco matrici nella GPU

cudaMalloc((void**)&T_dev.old,DIM*DIM*sizeof(float));//}
	
cudaMemcpy(T_dev.old,T_iniziale,DIM*DIM*sizeof(float),cudaMemcpyHostToDevice);
float *dev_fissata;

cudaMalloc((void**)&dev_fissata,DIM*DIM*sizeof(float));
cudaMemcpy(dev_fissata,T_fissata,DIM*DIM*sizeof(float),cudaMemcpyHostToDevice);

cudaMalloc((void**)&T_dev.upd,DIM*DIM*sizeof(float));
stampa_file (T_iniziale,DIM); 
float ** storage;
storage=(float**)malloc(durata*sizeof(float*));


for(n=0;n<durata;n++){
time_step(&T_dev,dev_fissata);
//cudaMemcpy(T_dev.old,T_dev.upd,DIM*DIM*sizeof(float),cudaMemcpyDeviceToDevice);
storage[n]=(float*)malloc(DIM*DIM*sizeof(float));
cudaMemcpy(storage[n],T_dev.upd,DIM*DIM*sizeof(float),cudaMemcpyDeviceToHost);
			}
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
stampa_file2(storage);
printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);


/*				serviva a controllaregli indici della griglia
int indici[DIM*DIM];int j,offset;printf("stampo gli indici \n");
for(i=0;i<DIM;i++){for(j=0;j<DIM;j++){offset=i+(j*DIM); printf("%d	",offset);}
		printf("\n");}
*/
cudaFree(T_dev.upd);
cudaFree(dev_fissata);
cudaFree(T_dev.old);



return 0;
}
