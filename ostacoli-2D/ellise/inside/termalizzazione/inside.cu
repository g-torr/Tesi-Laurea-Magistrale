#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cuda_runtime_api.h>
/*#include <helper_cuda.h>
#include <helper_functions.h>*/  
const float size=(float)3.;
const float half_size=size/(float)2.;
#define tau (float)0.06 // costante nel processo di O-U
#define D (float)0.1
#define dt (float)0.01
#define mobility (float)1
#define N 100
const int durata=size*size*2/D;
const float tsalva= durata/10;
const int blocks=128;
const int threads=1024;
const int particelle=blocks*threads;
struct point{	float x;
		float y;};


struct configurazione{	point*	eta;
			point*	r;
			point*	forza;
			bool* 	inside;};
__constant__ point vertice[N+1];
__constant__ point ottimizza[4+1];
__constant__ float costanti[2];

void creo_cartelle()					{system( "rm -rf ./posizione" );system( "rm -rf ./forza" );
							mkdir("posizione",0700);mkdir("forza",0700);}

void	stampa(configurazione stato, int i){int id;FILE*f;FILE*g;char indirizzo_posizione[50];char indirizzo_forza[50];
					sprintf(indirizzo_posizione,"./posizione/dati_%d",i);
					sprintf(indirizzo_forza,"./forza/forza_%d",i);
	
					f=fopen(indirizzo_posizione,"w");g=fopen(indirizzo_forza,"w");
		for(id=0;id<particelle;id++){fprintf(f,"%f		%f \n",stato.r[id].x,stato.r[id].y); //printf("sto scrivendo %d\n",id);
					fprintf(g,"%f		%f \n",stato.forza[id].x,stato.forza[id].y);}
					fclose(f);fclose(g);}

__global__ void setup_random_kernel(curandState* stato){
			int id=threadIdx.x+ blockIdx.x*blockDim.x; 
			curand_init(1234,id,0,&stato[id]); //inizializzo lo stato che mi genera i numeri casuali curand_init (seed,sequence,offset, curandState_t *state). per evitare ogni tipo di problemi qui si è scelto lo stesso seme per tutti i thread, ogni thread avrà sequenza differente: ovvero partendo dallo stesso seme si scartano i primi id*2^67 numeri casuali) a meno che il successivo  codice che si utilizza richiede 2^67 numeri casuali, noi siamo tranquilli che non ci siano overlap di numeri casuali. Si potrebbe cambiare il seme per ogni indice, fissare la sequenza a 0: si  potrebbero avere problemi di correlazione tra le sequenze di numeri con semi diversi(cosa molto rara). Dato che effettivamente l'idea di sprecare 2*67 numeri casuali  per ogni thread mi sembra una follia, forse la soluzione più semplice potrebbe essere quella di giocare sull'offset, in effetti se si fissano seme e sequence, mentre si  fissa il parametro offset= k*id, con k un qualsiasi numero > #numeri casuali che uso in ogni  tread  dovrei essere tranquillo che non si abbiano sovrapposizioni
			}


__device__ point segmento_vicino(point P){ //P è il punto rispetto al quale viene cercato il segmento più vicino,  la funzione restituisce closest, che è il punto sul perimetro del poligono più vicino a P, d è tale distanza minima.
	int i=0;float t;float d; 
	 point closest;
	t = ((P.x-vertice[0].x)*(vertice[1].x-vertice[0].x)+(P.y-vertice[0].y)*(vertice[0+1].y-vertice[0].y))/
								((vertice[0+1].x-vertice[0].x)*(vertice[0+1].x-vertice[0].x)+(vertice[0+1].y-vertice[0].y)*(vertice[0+1].y-vertice[0].y));

	 if(t<(float)0.0){t=(float)0.0;}
 	 if(t>(float)1.0){t=(float)1.0;} 


    closest.x = vertice[0].x+ (vertice[0+1].x-vertice[0].x)*t; 
    closest.y = vertice[0].y+ (vertice[0+1].y-vertice[0].y)*t;  
	d=pow(P.x- closest.x,2)+pow(P.y- closest.y,2);

	point temp;
	float d_temp;
 	for(i=1;i<N;i++){
    t = ((P.x-vertice[i].x)*(vertice[i+1].x-vertice[i].x)+(P.y-vertice[i].y)*(vertice[i+1].y-vertice[i].y))/
								((vertice[i+1].x-vertice[i].x)*(vertice[i+1].x-vertice[i].x)+(vertice[i+1].y-vertice[i].y)*(vertice[i+1].y-vertice[i].y));

    
    if(t<(float)0.0){t=(float)0.0;}
    if(t>(float)1.0){t=(float)1.0;} 

	
    temp.x = vertice[i].x+ (vertice[i+1].x-vertice[i].x)*t; 
    temp.y =vertice[i].y+ (vertice[i+1].y-vertice[i].y)*t;  
	d_temp=(P.x-temp.x)*(P.x-temp.x)+(P.y-temp.y)*(P.y-temp.y);
			if(d_temp<d){ closest.x=temp.x; closest.y=temp.y;d=d_temp;}
			}

	return closest;
     }

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
__device__ inline float	isLeft( point P0, point P1, point P2 ){
    							return ( (P1.x - P0.x) * (P2.y - P0.y)- (P2.x -  P0.x) * (P1.y - P0.y) );
				}


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
__device__ int wn_PnPoly( point P, int n )
{
    int    wn = 0;    // the  winding number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
        if (vertice[i].y <= P.y) {          // start y <= P.y
            if (vertice[i+1].y  > P.y)      // an upward crossing
                 if (isLeft( vertice[i], vertice[i+1], P) > 0)  // P left of  edge
                     ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if (vertice[i+1].y  <= P.y)     // a downward crossing
                 if (isLeft( vertice[i], vertice[i+1], P) < 0)  // P right of  edge
                     --wn;            // have  a valid down intersect
        }
    }
    return wn;
}




__device__ bool is_in(point P){int wn=wn_PnPoly(P,N);
					//printf("il wn=%d \n",wn);
					if (wn>0) return true;
							else return false;
						}



__global__ void	inizializza(configurazione stato,curandState* gen_random){point temp;
							int id=threadIdx.x+blockIdx.x*blockDim.x;	
						
							do {
							temp.x=-half_size+size*curand_uniform(&gen_random[id]);
							temp.y=-half_size+size*curand_uniform(&gen_random[id]);}
							
								while(is_in(temp)); 		 	     							stato.r[id].x=temp.x;
						stato.r[id].y=temp.y;
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
			//if((r.x>ottimizza[0].x)&&(r.x<ottimizza[1].x)&&(r.y<ottimizza[2].x)&&(r.y>ottimizza[0].y)){		
					bool in=is_in(r);
					if(!in) stato.inside[id]=true;
					else {stato.inside[id]=false;f.x=0.f; f.y=0.f;stato.forza[id]=f;	
						/*if(r.x>half_size)	r.x=r.x-size;
						else if (r.x<-half_size) r.x=r.x+size;
						if	(r.y>half_size)	r.y=r.y-size;
						else if (r.y<-half_size)	r.y=r.y+size;*/}
												
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
				point r,v,f;
				int k=threadIdx.x+blockIdx.x*blockDim.x;
				while(k<count){
				
				int id=ids[k]; 	
				v=stato.eta[id];	r=stato.r[id];
				f=r;
				r=segmento_vicino(r);
				f.x=(f.x-r.x)/(dt*mobility);
				f.y=(f.y-r.y)/(dt*mobility); 
				stato.eta[id]=v;	stato.forza[id]=f;	stato.r[id]=r;		
				k=k+(gridDim.x*blockDim.x);}
					}

void ottimizza_geometria(point* scheletro, point* temp_ottimizza){int i;
				
					float min_x,min_y,max_x,max_y;
					temp_ottimizza[0].x= scheletro[0].x;
					temp_ottimizza[0].y= scheletro[0].y;
					temp_ottimizza[1].x= scheletro[0].x;
					temp_ottimizza[2].y= scheletro[0].y;
					for(i=1;i<N;i++){
						min_x=scheletro[i].x;
						min_y=scheletro[i].y;
						max_x=scheletro[i].x;
						max_y=scheletro[i].y;
						if (temp_ottimizza[0].x>min_x) temp_ottimizza[0].x=min_x;
						if (temp_ottimizza[1].x<max_x) temp_ottimizza[1].x=max_x; 						
						if (temp_ottimizza[0].y>min_y) temp_ottimizza[0].y=min_y;
						else if (temp_ottimizza[2].y<max_y) temp_ottimizza[2].y=max_y;}
						
						temp_ottimizza[1].y=temp_ottimizza[0].y;
						temp_ottimizza[2].x=temp_ottimizza[1].x;
						temp_ottimizza[3].x=temp_ottimizza[0].x;
						temp_ottimizza[3].y=temp_ottimizza[2].y;
						temp_ottimizza[4].x=temp_ottimizza[0].x;
						temp_ottimizza[4].y=temp_ottimizza[0].y;

}
void alloco_punti(point* scheletro){  int i=0; float a,b;
					FILE*f;f=fopen("input.dat","r");
					while(fscanf(f,"%f" "%f",&a,&b)>0)i++;
					rewind(f); 
					if(N!=i){printf("il tuo file di input non è consistente con il numero di vertici aspettato"); 
							exit(-1);}	
					
					for(i=0;i<N;i++){
							fscanf(f,"%f %f",&scheletro[i].x,&scheletro[i].y);
							}
					fclose(f);
					scheletro[N].x=scheletro[0].x;
					scheletro[N].y=scheletro[0].y;//sto creando una geometria con un vertice in più che coincide con il primo		
					}

main(){
creo_cartelle();
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);

int i,t;
point* temp_scheletro;
point* temp_ottimizza;
temp_scheletro=(point*)malloc((N +1)*sizeof(point));
alloco_punti(temp_scheletro);
temp_ottimizza=(point*)malloc((4+1)*sizeof(point));
ottimizza_geometria(temp_scheletro,temp_ottimizza);

printf("per usare la constant memory il numero di vertici  deve essere fissato al tempo di compilazione, sicuro di avere %d vertici?\n",N);
cudaMemcpyToSymbol(ottimizza,temp_ottimizza,(4+1)*sizeof(point));
cudaMemcpyToSymbol(vertice,temp_scheletro,(N+1)*sizeof(point));

float costanti_host[2];
costanti_host[0]=dt/tau;
costanti_host[1]=(float)sqrt(2*D* dt)/tau;
cudaMemcpyToSymbol(costanti,costanti_host,2*sizeof(float));
int numero_passi=(int)durata/dt;
int passi_salvataggio =(int)(tsalva/dt);


/*
for(i=0;i<4;i++){
temp_ottimizza[i].x=0;
temp_ottimizza[i].y=0;}
srand(10);
cudaMemcpyFromSymbol(temp_ottimizza,ottimizza,(4+1)*sizeof(point));
for(i=0;i<4;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,temp_ottimizza[i].x,temp_ottimizza[i].y,temp_ottimizza[i+1].x,temp_ottimizza[i+1].y);}
*/
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
			correct<<<count/1024,1024>>>(dev_stato,generatori_random,dev_ids,count);					
							
			   if((t> 2**i)&&(t>0)){printf("siamo a %d/9 \n",i);
			cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
			cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
				stampa(stato,i);i++;
	
}
//cudaDeviceSynchronize();
			}

 

cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
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


