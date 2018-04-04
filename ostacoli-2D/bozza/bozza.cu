#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>

const int particelle=100;
const int N=4;
const float size=3.0;

struct segmenti	{	float x0;
			float x1;
			float y0;
			float y1;		//poichè il segmento è definito a partire da 2 punti, devo passare 4 coordinate per segmento
			bool up;	// questo parametro deve essere assegnato a mano quando si costruisce la geometria, mi dice se devo 						guardare sopra o sotto rispetto a ogni segmento
		};


struct vertici{	float x;
		float y;
			};

struct geometria{	vertici*	punti;
			segmenti*	lati;};

struct stato{		float*	x;
			float*	y;
			float*	eta_x;
			float* 	eta_y;};
void	stampa(stato configurazione){int i;FILE*f;
					f=fopen("output.dat","w");
		for(i=0;i<particelle;i++){fprintf(f,"%f		%f \n",configurazione.x[i],configurazione.y[i]);}
					fclose(f);}

float  		vincolo(segmenti lato,float x, float y){
						float f=y-lato.y1-(((lato.y1-lato.y0)/(lato.x1-lato.x0))*(x-lato.x1));
							return f;}

bool is_present(geometria scheletro,float x, float y){int i;bool flag;
				for(i=0;i<N;i++){
					if ((vincolo(scheletro.lati[i],x,y)>0)==scheletro.lati[i].up) flag=true;
					else return false;
							}return true;}



void	inizializza(stato* configurazione,geometria scheletro){int i=0;float temp_x,temp_y;		
						while(i<particelle){
							temp_x=size*(float)rand()/RAND_MAX;
							temp_y=size*(float)rand()/RAND_MAX;
							if (is_present(scheletro,temp_x,temp_y)){
										(*configurazione).x[i]=temp_x;(*configurazione).y[i]=temp_y;
											i++;}
 										}					 									 															
					
											}



void	crea_body(geometria* scheletro){ int i;int flag;
						(*scheletro).lati=(segmenti*)malloc(N*sizeof(segmenti));
						for(i=0;i<N-1;i++)	{(*scheletro).lati[i].x0=(*scheletro).punti[i].x;
									(*scheletro).lati[i].y0=(*scheletro).punti[i].y;
									(*scheletro).lati[i].x1=(*scheletro).punti[i+1].x;
									(*scheletro).lati[i].y1=(*scheletro).punti[i+1].y;
									flag=2;	
									while((flag!=1)&&(flag!=-1))	{
									printf("se il vincolo del lato %d è >0 premi 1; altrimenti se è <0 premi -1\n",i);												
									scanf("%d",&flag); }
									if (flag==1)	(*scheletro).lati[i].up=true;
									else if (flag==-1) (*scheletro).lati[i].up=false;
 									}
						(*scheletro).lati[N-1].x0=(*scheletro).punti[N-1].x;
						(*scheletro).lati[N-1].y0=(*scheletro).punti[N-1].y;
						(*scheletro).lati[N-1].x1=(*scheletro).punti[0].x;
						(*scheletro).lati[N-1].y1=(*scheletro).punti[0].y;
					flag=2;	
					while((flag!=1)&&(flag!=-1))	{
					printf("se il vincolo del lato %d è >0 premi 1; altrimenti se è <0 premi -1\n",i);
					scanf("%d",&flag);		}
					if (flag==1)	(*scheletro).lati[i].up=true;
					else if (flag==-1) (*scheletro).lati[i].up=false;	}


void alloco_vertici(geometria* scheletro){  int i; 
					(*scheletro).punti=(vertici*)malloc(N*sizeof(vertici));
					for(i=0;i<N;i++){
					printf("inserisci ascissa e ordinata del %d° vertice\n",i);
					scanf("%f",&(*scheletro).punti[i].x);scanf("%f",&(*scheletro).punti[i].y);
							}
						}

main(){
geometria scheletro;
//numero di vertici dell'oggetto
alloco_vertici(&scheletro);
crea_body(&scheletro);
srand(10);
int i;

for(i=0;i<N;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,scheletro.lati[i].x0,scheletro.lati[i].y0,scheletro.lati[i].x1,scheletro.lati[i].y1);}

stato configurazione;

configurazione.x=(float*)malloc(particelle*sizeof(float));
configurazione.y=(float*)malloc(particelle*sizeof(float));

inizializza(&configurazione,scheletro);

stampa(configurazione);
}


