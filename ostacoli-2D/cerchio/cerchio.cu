#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

const int particelle=10000;
//const int N=4;
const float size=3.0;
const float half_size=size/2.;
#define tau 0.06 // costante nel processo di O-U
#define D 0.1
#define dt 0.01
#define mobility 1
const int durata=100;
const float tsalva= durata/10.;
const float raggio=1.; //ricordati che deve essere minore di half_size

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
		for(id=0;id<particelle;id++){fprintf(f,"%f		%f \n",stato.r[id].x,stato.r[id].y); //printf("sto scrivendo %d\n",id);
					fprintf(g,"%f		%f \n",stato.forza[id].x,stato.forza[id].y);}
					fclose(f);fclose(g);}
point nuova_posizione(float raggio,point old){point nuovo;float d=sqrt((old.x*old.x)+(old.y*old.y));
						nuovo.x=raggio*old.x/d;
						nuovo.y=raggio*old.y/d;
						
						return nuovo;}

double randn (double mu, double sigma){
  					double U1, U2, W, mult;
  					static double X1, X2;
  					static int call = 0;
 
  					if (call == 1)
    					{
      						call = !call;
      						return (mu + sigma * (double) X2);
   					 }
 
  					do
   					 {
     						 U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      						U2 = -1 + ((double) rand () / RAND_MAX) * 2;
     						 W = pow (U1, 2) + pow (U2, 2);
   					 }
  						while (W >= 1 || W == 0);
 
				mult = sqrt ((-2 * log (W)) / W);
				X1 = U1 * mult;
  				X2 = U2 * mult;
 
  				call = !call;
 
  return (mu + sigma * (double) X1);
}





bool is_in(float raggio,point P){float d=sqrt(P.x*P.x+P.y*P.y);
							if (d<raggio) return true;
							else return false;
						}



void	inizializza(configurazione* stato,float raggio){int i=0;point temp;	
						(*stato).eta=(point*)malloc(particelle*sizeof(point));
						(*stato).r  =(point*)malloc(particelle*sizeof(point));
						(*stato).forza=(point*)calloc(particelle,sizeof(point));
						while(i<particelle){
							temp.x=-half_size+size*(float)rand()/RAND_MAX;
							temp.y=-half_size+size*(float)rand()/RAND_MAX;
							if (!is_in(raggio,temp)){
										(*stato).r[i].x=temp.x;
										(*stato).r[i].y=temp.y;
										(*stato).eta[i].x=0;	
										(*stato).eta[i].y=0;
										i++;}
 										}					 		
											}
void	evolvi(configurazione stato,float raggio){	
			for(int i=0;i<particelle;i++){
						stato.eta[i].x=stato.eta[i].x-(1/tau)*stato.eta[i].x*dt+sqrt(D)*randn(0,sqrt(2.))*sqrt(dt)/    tau;      				//	printf("%f\n",stato.eta[i].x);
						stato.eta[i].y=stato.eta[i].y-(1/tau)*stato.eta[i].y*dt+sqrt(D)*randn(0,sqrt(2.))*sqrt(dt)/tau;
						stato.r[i].x=stato.r[i].x+stato.eta[i].x*dt;
						stato.r[i].y=stato.r[i].y+stato.eta[i].y*dt;		
					if (is_in(raggio,stato.r[i])){
										stato.forza[i]=stato.r[i];
										stato.r[i]=nuova_posizione(raggio,stato.r[i]);
										stato.forza[i].x=(stato.forza[i].x-stato.r[i].x)/(dt*mobility);
										stato.forza[i].y=(stato.forza[i].y-stato.r[i].y)/(dt*mobility);
													}


					else {	if(stato.r[i].x>half_size)	stato.r[i].x=stato.r[i].x-size;
							else if (stato.r[i].x<-half_size) stato.r[i].x=stato.r[i].x+size;
						if	(stato.r[i].y>half_size)	stato.r[i].y=stato.r[i].y-size;
							else if (stato.r[i].y<-half_size)	stato.r[i].y=stato.r[i].y+size;
						stato.forza[i].x=0.; stato.forza[i].y=0.;
						}

				}	
}


main(){
clock_t t1 = clock();
if(raggio>half_size){printf("attenzione!! raggio deve essere maggiore di half size");return 1;}


srand(10);
int t,i=0;
/*
for(i=0;i<4;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,ottimizza.vertici[i].x,ottimizza.vertici[i].y,ottimizza.vertici[i+1].x,ottimizza.vertici[i+1].y);}*/
configurazione stato;
inizializza(&stato,raggio);
int numero_passi=(int)durata/dt;
int passi_salvataggio =(int)(tsalva/dt);
for(t=0;t<numero_passi;t++){
			evolvi(stato,raggio);
			   if((t% passi_salvataggio==0)&&(t>0)){stampa(stato,i);i++;}
			}
point a;a.x=.25;a.y=-1.5;
point b=nuova_posizione(raggio,a);
point c;c.x=a.x-b.x; c.y=a.y-b.y;

printf("il punto (%f,%f) si trova a distanza%f dal punto (%f,%f)\n",a.x,a.y,sqrt(c.x*c.x+c.y*c.y),b.x,b.y);
//stampa(stato,1);
 clock_t t2 = clock();
double time_sec = 
       (double)(t2-t1)/(double)(CLOCKS_PER_SEC); 
 
    printf("Time (sec): %lf\n",time_sec); 
}


