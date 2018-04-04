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
#include <ftw.h>
#include <unistd.h>
/*#include <helper_cuda.h>
#include <helper_functions.h>*/  
const float size=3.0;
const float half_size=size/2.;
#define tau 0.06 // costante nel processo di O-U
#define D 0.1
#define dt 0.01
#define mobility 1
#define N 4
const int durata=size*size*2/D;
const float tsalva= durata/10.;
const int blocks=100;
const int threads=1024;
const int particelle=blocks*threads;
struct point{	float x;
		float y;};


struct configurazione{	point*	eta;
			point*	r;
			point*	forza;};
__constant__ point vertice[N+1];
__constant__ point ottimizza[4+1];

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

//				non la uso più!!
__device__ void segmento_vicino(point P,point* closest,int id){ //P è il punto rispetto al quale viene cercato il segmento più vicino,  la funzione restituisce closest, che è il punto sul perimetro del poligono più vicino a P, d è tale distanza minima.
	int i=0;float t;float d; 
	 
	t = ((P.x-vertice[0].x)*(vertice[0+1].x-vertice[0].x)+(P.y-vertice[0].y)*(vertice[0+1].y-vertice[0].y))/
								((vertice[0+1].x-vertice[0].x)*(vertice[0+1].x-vertice[0].x)+(vertice[0+1].y-vertice[0].y)*(vertice[0+1].y-vertice[0].y));

	 if(t<0.0){t=0.0;}
 	 if(t>1.0){t=1.0;} 


    closest[id].x = vertice[0].x+ (vertice[0+1].x-vertice[0].x)*t; 
    closest[id].y = vertice[0].y+ (vertice[0+1].y-vertice[0].y)*t;  
	d=(P.x- closest[id].x)*(P.x- closest[id].x)+(P.y- closest[id].y)*(P.y- closest[id].y);

	point temp;
	float d_temp;
 	for(i=1;i<N;i++){
    t = ((P.x-vertice[i].x)*(vertice[i+1].x-vertice[i].x)+(P.y-vertice[i].y)*(vertice[i+1].y-vertice[i].y))/
								((vertice[i+1].x-vertice[i].x)*(vertice[i+1].x-vertice[i].x)+(vertice[i+1].y-vertice[i].y)*(vertice[i+1].y-vertice[i].y));

    
    if(t<0.0){t=0.0;}
    if(t>1.0){t=1.0;} 

	
    temp.x = vertice[i].x+ (vertice[i+1].x-vertice[i].x)*t; 
    temp.y =vertice[i].y+ (vertice[i+1].y-vertice[i].y)*t;  
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


__device__ bool interseca(float*t, float*T,point P_old,point P_nuovo,point V1,point V2){float x1,x2,X1,X2,y1,y2,Y1,Y2;
						X1=P_old.x;	X2=P_nuovo.x;	x1=V1.x;	x2=V2.x;
						Y1=P_old.y;	Y2=P_nuovo.y;	y1=V1.y;	y2=V2.y;
					*T=(x2*(-y1 + Y1) + X1*(y1 - y2) + x1*(-Y1 + y2))/((X1 - X2)*(y1 - y2) - (x1 - x2)*(Y1 - Y2));
					*t=(X2*(y1 - Y1) + x1*(Y1 - Y2) + X1*(-y1 + Y2))/(-((X1 - X2)*(y1 - y2)) + (x1 - x2)*(Y1 - Y2));					
					if((*t>=0.0)&&(*t<=1.0)&&(*T>=0.)&&(*T<=1.)) return true;
     						else return false;	}

__device__ bool closest(point P_old,point P, point* nuovo){float t,temp_t,T,temp_T;int i;bool flag=false;
				T=1;				
				for(i=0;i<N;i++){ 
				if (interseca(&temp_t,&temp_T,P_old,P,vertice[i],vertice[i+1])){flag=true;if (temp_T<T) T=temp_T;}
						}
				
				(*nuovo).x=(P.x-P_old.x)*(T-0.00000001)+P_old.x;
				(*nuovo).y=(P.y-P_old.y)*(T-0.00000001)+P_old.y;
				return flag;	
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
						stato.eta[id].x=0;	
						stato.eta[id].y=0;
										
 															 		
}
__global__ void	evolvi(configurazione stato,curandState* gen_random){	
						 
						int id=threadIdx.x+blockIdx.x*blockDim.x; point proposta;	
						stato.eta[id].x=stato.eta[id].x-(1/tau)*stato.eta[id].x*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau;      			
						stato.eta[id].y=stato.eta[id].y-(1/tau)*stato.eta[id].y*dt+sqrt(D)*curand_normal(&gen_random[id])*sqrt(2.)*sqrt(dt)/tau;			
						proposta.x=stato.r[id].x+stato.eta[id].x*dt;
						proposta.y=stato.r[id].y+stato.eta[id].y*dt;	
						stato.forza[id]=proposta;						
						bool flag=closest(stato.r[id],proposta,&stato.r[id]);
		if((stato.r[id].x>ottimizza[0].x)&&(stato.r[id].x<ottimizza[1].x)&&(stato.r[id].y<ottimizza[2].x)&&(stato.r[id].y>ottimizza[0].y))
		{ 
					if(flag){
				//	if (is_in(proposta)){segmento_vicino(proposta,&proposta,0); 				
							
						//	stato.r[id]=proposta;
							stato.forza[id].x=(stato.forza[id].x-stato.r[id].x)/(dt*mobility);
							stato.forza[id].y=(stato.forza[id].y-stato.r[id].y)/(dt*mobility); }
					else{stato.forza[id].x=0.; stato.forza[id].y=0.;}
		}	
		else{
			stato.forza[id].x=0.; stato.forza[id].y=0.; 
			if(stato.r[id].x>half_size) 	stato.r[id].x=stato.r[id].x-size;
				else if (stato.r[id].x<-half_size) stato.r[id].x=stato.r[id].x+size;
			if	(stato.r[id].y>half_size)	stato.r[id].y=stato.r[id].y-size;
				else if (stato.r[id].y<-half_size)	stato.r[id].y=stato.r[id].y+size;
					}
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
					FILE*f;f=fopen("input2.dat","r");
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

main(){int deviceCount;
cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess)
    {
        printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
        printf("Result = FAIL\n");
        exit(EXIT_FAILURE);
    }
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
/*
stato.eta=(point*)malloc(particelle*sizeof(point));
stato.r  =(point*)malloc(particelle*sizeof(point));
stato.forza=(point*)calloc(particelle,sizeof(point));*/


curandState * generatori_random;//alloco il generatore dei numeri random sul device
if (numero_passi> pow(2,67)) printf("ATTENZIONE! ricontrollare il generatore di numeri casuali"); 
cudaMalloc((void**)&generatori_random,particelle*sizeof(curandState));
setup_random_kernel<<<blocks,threads>>>(generatori_random);


inizializza<<<blocks/4,threads,0,stream0>>>(dev_stato_0,generatori_random);
inizializza<<<blocks/4,threads,0,stream1>>>(dev_stato_1,generatori_random);
inizializza<<<blocks/4,threads,0,stream2>>>(dev_stato_2,generatori_random);
inizializza<<<blocks/4,threads,0,stream3>>>(dev_stato_3,generatori_random);
i=0;
for(t=0;t<numero_passi;t++){
//for(t=0;t<10;t++){
			evolvi<<<blocks/4,threads,0,stream0>>>(dev_stato_0,generatori_random);
			evolvi<<<blocks/4,threads,0,stream1>>>(dev_stato_1,generatori_random);
			evolvi<<<blocks/4,threads,0,stream2>>>(dev_stato_2,generatori_random);
			evolvi<<<blocks/4,threads,0,stream3>>>(dev_stato_3,generatori_random);
	
	   if((t% passi_salvataggio==0)&&(t>0)){printf("siamo a %d/9 \n",i);
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
/*
cudaMemcpy(stato.r,dev_stato.r,particelle*sizeof(point),cudaMemcpyDeviceToHost);
cudaMemcpy(stato.forza,dev_stato.forza,particelle*sizeof(point),cudaMemcpyDeviceToHost);*/
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


