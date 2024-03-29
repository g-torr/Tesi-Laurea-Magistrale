#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
#include <cuda_runtime_api.h>
//#include <helper_cuda.h>

const float size=(float)3.;
const float half_size=size/(float)2.;
#define tau (float)0.06 // costante nel processo di O-U
#define D (float)0.1
#define dt (float)0.01
#define mobility (float)1
const int durata=size*size*2/D;
const float tsalva= durata/10.;
const int blocks=128;
const int threads=1024;
const int particelle=blocks*threads;
struct point{	float x;
		float y;
			};

struct geometria{	point*	vertici;};

struct configurazione{	point*	eta;
			point*	r;
			point*	forza;};
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

__device__ void segmento_vicino(point * V, int N, point P,point* closest,int id){ //P è il punto rispetto al quale viene cercato il segmento più vicino,  la funzione restituisce closest, che è il punto sul perimetro del poligono più vicino a P, d è tale distanza minima.
	int i=0;float t;float d; 
	 
	t = ((P.x-V[0].x)*(V[0+1].x-V[0].x)+(P.y-V[0].y)*(V[0+1].y-V[0].y))/
								((V[0+1].x-V[0].x)*(V[0+1].x-V[0].x)+(V[0+1].y-V[0].y)*(V[0+1].y-V[0].y));

	 if(t<0.0){t=0.0;}
 	 if(t>1.0){t=1.0;} 


    closest[id].x = V[0].x+ (V[0+1].x-V[0].x)*t; 
    closest[id].y = V[0].y+ (V[0+1].y-V[0].y)*t;  
	d=(P.x- closest[id].x)*(P.x- closest[id].x)+(P.y- closest[id].y)*(P.y- closest[id].y);

	point temp;
	float d_temp;
 	for(i=1;i<N;i++){
    t = ((P.x-V[i].x)*(V[i+1].x-V[i].x)+(P.y-V[i].y)*(V[i+1].y-V[i].y))/
								((V[i+1].x-V[i].x)*(V[i+1].x-V[i].x)+(V[i+1].y-V[i].y)*(V[i+1].y-V[i].y));

    
    if(t<0.0){t=0.0;}
    if(t>1.0){t=1.0;} 

	
    temp.x = V[i].x+ (V[i+1].x-V[i].x)*t; 
    temp.y =V[i].y+ (V[i+1].y-V[i].y)*t;  
	d_temp=(P.x-temp.x)*(P.x-temp.x)+(P.y-temp.y)*(P.y-temp.y);
			if(d_temp<d){ closest[id].x=temp.x; closest[id].y=temp.y;d=d_temp;}
				}

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
__device__ int wn_PnPoly( point P, point* V, int n )
{
    int    wn = 0;    // the  winding number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
        if (V[i].y <= P.y) {          // start y <= P.y
            if (V[i+1].y  > P.y)      // an upward crossing
                 if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
                     ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if (V[i+1].y  <= P.y)     // a downward crossing
                 if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
                     --wn;            // have  a valid down intersect
        }
    }
    return wn;
}




__device__ bool is_in(geometria scheletro,point P,int N){int wn=wn_PnPoly(P,scheletro.vertici,N);
					//printf("il wn=%d \n",wn);
					if (wn>0) return true;
							else return false;
						}



__global__ void	inizializza(configurazione stato,geometria scheletro,int N,curandState* gen_random){point temp;
							int id=threadIdx.x+blockIdx.x*blockDim.x;	
						
							do {
							temp.x=-half_size+size*curand_uniform(&gen_random[id]);
							temp.y=-half_size+size*curand_uniform(&gen_random[id]);}
							
								while(is_in(scheletro,temp,N)); 						stato.r[id].x=temp.x;
						stato.r[id].y=temp.y;
						stato.eta[id].x=0;	
						stato.eta[id].y=0;
										
 															 		
}
__global__ void	evolvi(configurazione stato,geometria scheletro, int N,geometria ottimizza,curandState* gen_random){	
						 
						int id=threadIdx.x+blockIdx.x*blockDim.x; 	
						stato.eta[id].x=stato.eta[id].x-(1/tau)*stato.eta[id].x*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau;      			
						stato.eta[id].y=stato.eta[id].y-(1/tau)*stato.eta[id].y*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau;			
						stato.r[id].x=stato.r[id].x+stato.eta[id].x*dt;
						stato.r[id].y=stato.r[id].y+stato.eta[id].y*dt;	
				if((stato.r[id].x>ottimizza.vertici[0].x)&&(stato.r[id].x<ottimizza.vertici[1].x)&&(stato.r[id].y<ottimizza.vertici[2].x)&&(stato.r[id].y>ottimizza.vertici[0].y)) {		
					if (is_in(scheletro,stato.r[id],N)){
										stato.forza[id]=stato.r[id];
										segmento_vicino(scheletro.vertici,N,stato.r[id],stato.r,id);
										__syncthreads();
										stato.forza[id].x=(stato.forza[id].x-stato.r[id].x)/(dt*mobility);
										stato.forza[id].y=(stato.forza[id].y-stato.r[id].y)/(dt*mobility); 
													}
					else {stato.forza[id].x=0.; stato.forza[id].y=0.;}

									}
		else {	if(stato.r[id].x>half_size)	stato.r[id].x=stato.r[id].x-size;
				else if (stato.r[id].x<-half_size) stato.r[id].x=stato.r[id].x+size;
			if	(stato.r[id].y>half_size)	stato.r[id].y=stato.r[id].y-size;
				else if (stato.r[id].y<-half_size)	stato.r[id].y=stato.r[id].y+size;
			stato.forza[id].x=0.; stato.forza[id].y=0.;	}

						
//if(id==5)printf(" la particella 5 si trova in %f \n",stato.r[id].x);
}

void ottimizza_geometria(geometria scheletro, geometria* ottimizza, int N){int i;
					(*ottimizza).vertici=(point*)malloc((4+1)*sizeof(point));
					float min_x,min_y,max_x,max_y;
					(*ottimizza).vertici[0].x= scheletro.vertici[0].x;
					(*ottimizza).vertici[0].y= scheletro.vertici[0].y;
					(*ottimizza).vertici[1].x= scheletro.vertici[0].x;
					(*ottimizza).vertici[2].y= scheletro.vertici[0].y;
					for(i=1;i<N;i++){
						min_x=scheletro.vertici[i].x;
						min_y=scheletro.vertici[i].y;
						max_x=scheletro.vertici[i].x;
						max_y=scheletro.vertici[i].y;
						if ((*ottimizza).vertici[0].x>min_x) (*ottimizza).vertici[0].x=min_x;
						if ((*ottimizza).vertici[1].x<max_x) (*ottimizza).vertici[1].x=max_x; 						
						if ((*ottimizza).vertici[0].y>min_y) (*ottimizza).vertici[0].y=min_y;
						else if ((*ottimizza).vertici[2].y<max_y) (*ottimizza).vertici[2].y=max_y;}
						
						(*ottimizza).vertici[1].y=(*ottimizza).vertici[0].y;
						(*ottimizza).vertici[2].x=(*ottimizza).vertici[1].x;
						(*ottimizza).vertici[3].x=(*ottimizza).vertici[0].x;
						(*ottimizza).vertici[3].y=(*ottimizza).vertici[2].y;
						(*ottimizza).vertici[4].x=(*ottimizza).vertici[0].x;
						(*ottimizza).vertici[4].y=(*ottimizza).vertici[0].y;

}
void alloco_punti(geometria* scheletro,int *N){  int i=0; float a,b;
					FILE*f;f=fopen("input2.dat","r");
					while(fscanf(f,"%f" "%f",&a,&b)>0)i++;
					*N=i;rewind(f);
					(*scheletro).vertici=(point*)malloc((*N +1)*sizeof(point));
					for(i=0;i<*N;i++){
							fscanf(f,"%f %f",&(*scheletro).vertici[i].x,&(*scheletro).vertici[i].y);
							}
					fclose(f);
					(*scheletro).vertici[*N].x=(*scheletro).vertici[0].x;
					(*scheletro).vertici[*N].y=(*scheletro).vertici[0].y;//sto creando una geometria con un vertice in più che coincide con il primo		
					}

main(){
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
int N,i,t;
geometria scheletro,dev_scheletro;
geometria ottimizza,dev_ottimizza;
alloco_punti(&scheletro,&N);
ottimizza_geometria(scheletro,&ottimizza,N);

cudaMalloc((point**)&dev_scheletro.vertici,(N+1)*sizeof(point));
cudaMalloc((point**)&dev_ottimizza.vertici,(4+1)*sizeof(point));
cudaMemcpy(dev_scheletro.vertici,scheletro.vertici,(N+1)*sizeof(point),cudaMemcpyHostToDevice);
cudaMemcpy(dev_ottimizza.vertici,ottimizza.vertici,(4+1)*sizeof(point),cudaMemcpyHostToDevice);
int numero_passi=(int)durata/dt;
int passi_salvataggio =(int)(tsalva/dt);

/*

for(i=0;i<N;i++){
scheletro.vertici[i].x=0;
scheletro.vertici[i].y=0;}
srand(10);
cudaMemcpy(scheletro.vertici,dev_scheletro.vertici,(N+1)*sizeof(point),cudaMemcpyDeviceToHost);*/

/*for(i=0;i<N;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,scheletro.vertici[i].x,scheletro.vertici[i].y,scheletro.vertici[i+1].x,scheletro.vertici[i+1].y);}*/
configurazione dev_stato,stato;//alloco lo stato del sistema
cudaMalloc((point**)&dev_stato.eta,particelle*sizeof(point));
cudaMalloc((point**)&dev_stato.r,particelle*sizeof(point));
cudaMalloc((point**)&dev_stato.forza,particelle*sizeof(point));

curandState * generatori_random;//alloco il generatore dei numeri random sul device
if (numero_passi> pow(2,67)) printf("ATTENZIONE! ricontrollare il generatore di numeri casuali"); 
cudaMalloc((void**)&generatori_random,particelle*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(generatori_random);



inizializza<<<blocks,threads>>>(dev_stato,dev_scheletro,N,generatori_random);
stato.eta=(point*)malloc(particelle*sizeof(point));
stato.r  =(point*)malloc(particelle*sizeof(point));
stato.forza=(point*)calloc(particelle,sizeof(point));


i=0;
for(t=0;t<numero_passi;t++){
			evolvi<<<blocks,threads>>>(dev_stato,dev_scheletro,N,dev_ottimizza,generatori_random);
			   if((t% passi_salvataggio==0)&&(t>0)){printf("siamo a %d/9 \n",i);
			cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
			cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
				stampa(stato,i);i++;
	
}
//	cudaDeviceSynchronize();	
			}
 
 

cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);
stampa(stato,i);

free(stato.r);
free(stato.forza);
free(stato.eta);
 cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
float tempo_trascorso; cudaEventElapsedTime(&tempo_trascorso,start,stop);
printf("il tempo necessario per eseguire il programma è %f ms\n",tempo_trascorso);
}


