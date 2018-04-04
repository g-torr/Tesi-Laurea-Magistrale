#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <ftw.h>
#include <unistd.h>

const int particelle=100000;
//const int N=4;
const double size=3.0;
const double half_size=size/2.;
#define tau 0.06 // costante nel processo di O-U
#define D 0.1
#define dt 0.01
#define mobility 1
const int durata=size*size*1/(2*D);
const double tsalva= durata/10.;

struct point{	double x;
		double y;
			};

struct geometria{	point*	vertici;};

struct configurazione{	point*	eta;
			point*	r;
			point*	forza;};
void creo_cartelle()					{system( "rm -rf ./posizione" );system( "rm -rf ./forza" );
							mkdir("posizione",0700);mkdir("forza",0700);}
void	stampa(configurazione stato, int i){int id;FILE*f;FILE*g;char indirizzo_posizione[50];char indirizzo_forza[50];
					sprintf(indirizzo_posizione,"./posizione/dati_%d",i);
					sprintf(indirizzo_forza,"./forza/forza_%d",i);
					
					f=fopen(indirizzo_posizione,"w");g=fopen(indirizzo_forza,"w");
		for(id=0;id<particelle;id++){fprintf(f,"%f		%f \n",stato.r[id].x,stato.r[id].y); //printf("sto scrivendo %d\n",id);
					fprintf(g,"%f		%f \n",stato.forza[id].x,stato.forza[id].y);}
					fclose(f);fclose(g);}

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


/*point segmento_vicino(point * V, int N, point P){ //P è il punto rispetto al quale viene cercato il segmento più vicino,  la funzione restituisce closest, che è il punto sul perimetro del poligono più vicino a P, d è tale distanza minima.
	int i=0;double t;double d; point closest;
	// t=(double*)malloc(N*sizeof(double));
	t = ((P.x-V[0].x)*(V[0+1].x-V[0].x)+(P.y-V[0].y)*(V[0+1].y-V[0].y))/
								((V[0+1].x-V[0].x)*(V[0+1].x-V[0].x)+(V[0+1].y-V[0].y)*(V[0+1].y-V[0].y));

	 if(t<0.0){t=0.0;}
 	 if(t>1.0){t=1.0;} 


    (closest).x = V[0].x+ (V[0+1].x-V[0].x)*t; 
    (closest).y = V[0].y+ (V[0+1].y-V[0].y)*t;  
	d=(P.x- (closest).x)*(P.x- (closest).x)+(P.y- (closest).y)*(P.y- (closest).y);

	point temp;
	double d_temp;
 	for(i=1;i<N;i++){
    t = ((P.x-V[i].x)*(V[i+1].x-V[i].x)+(P.y-V[i].y)*(V[i+1].y-V[i].y))/
								((V[i+1].x-V[i].x)*(V[i+1].x-V[i].x)+(V[i+1].y-V[i].y)*(V[i+1].y-V[i].y));

    
    if(t<0.0){t=0.0;}
    if(t>1.0){t=1.0;} 


    temp.x = V[i].x+ (V[i+1].x-V[i].x)*t; 
    temp.y =V[i].y+ (V[i+1].y-V[i].y)*t;  
	d_temp=(P.x-temp.x)*(P.x-temp.x)+(P.y-temp.y)*(P.y-temp.y);
			if(d_temp<d){ closest.x=temp.x; closest.y=temp.y;d=d_temp;}
				}

    return closest; }*/

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
inline double	isLeft( point P0, point P1, point P2 ){
    							return ( (P1.x - P0.x) * (P2.y - P0.y)- (P2.x -  P0.x) * (P1.y - P0.y) );
				}


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int
wn_PnPoly( point P, point* V, int n )
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

double interseca(point P_old,point P_nuovo,point V1,point V2){double x1,x2,X1,X2,y1,y2,Y1,Y2;double t; double T;
						X1=P_old.x;	X2=P_nuovo.x;	x1=V1.x;	x2=V2.x;
						Y1=P_old.y;	Y2=P_nuovo.y;	y1=V1.y;	y2=V2.y;
					T=(x2*(-y1 + Y1) + X1*(y1 - y2) + x1*(-Y1 + y2))/((X1 - X2)*(y1 - y2) - (x1 - x2)*(Y1 - Y2));
					t=(X2*(y1 - Y1) + x1*(Y1 - Y2) + X1*(-y1 + Y2))/(-((X1 - X2)*(y1 - y2)) + (x1 - x2)*(Y1 - Y2));		 
			
					if((t>=0.0)&&(t<=1.0)&&(T>=0.)&&(T<1.)) return T;
     						else return 2;	}

bool closest(point P_old,point* P,geometria scheletro,int N){double t,temp_t,T,temp_T;int i;bool flag=false;
				T=2;				
				for(i=0;i<N;i++){ 
				temp_T=interseca(P_old,*P,scheletro.vertici[i],scheletro.vertici[i+1]);
										if (temp_T<T){flag=true;T=temp_T;}

						}
				if(flag==true){
				(*P).x=((*P).x-P_old.x)*(T-0.000001)+P_old.x;
				(*P).y=((*P).y-P_old.y)*(T-0.000001)+P_old.y;
						}
			/*	else{(*P).x=0;
				(*P).y=;}*/
				return flag;	
					} 


bool is_in(geometria scheletro,point P,int N){int wn=wn_PnPoly(P,scheletro.vertici,N);
					//printf("il wn=%d \n",wn);
					if (wn>0) return true;
							else return false;
						}



void	inizializza(configurazione* stato,geometria scheletro,int N){int i=0;point temp;	
						(*stato).eta=(point*)malloc(particelle*sizeof(point));
						(*stato).r  =(point*)malloc(particelle*sizeof(point));
						(*stato).forza=(point*)calloc(particelle,sizeof(point));
						while(i<particelle){
							temp.x=-half_size+size*(double)rand()/RAND_MAX;
							temp.y=-half_size+size*(double)rand()/RAND_MAX;
							if (!is_in(scheletro,temp,N)){
										(*stato).r[i].x=temp.x;
										(*stato).r[i].y=temp.y;
										(*stato).eta[i].x=0;	
										(*stato).eta[i].y=0;
										i++;}
 										}					 		
											}
void	evolvi(configurazione stato,geometria scheletro, int N,geometria ottimizza){	point proposta;
			for(int i=0;i<particelle;i++){  
						stato.eta[i].x=stato.eta[i].x-(1/tau)*stato.eta[i].x*dt+sqrt(D)*randn(0,sqrt(2.))*sqrt(dt)/    tau;      				//	printf("%f\n",stato.eta[i].x);
						stato.eta[i].y=stato.eta[i].y-(1/tau)*stato.eta[i].y*dt+sqrt(D)*randn(0,sqrt(2.))*sqrt(dt)/tau;
						proposta.x=stato.r[i].x+stato.eta[i].x*dt;
						proposta.y=stato.r[i].y+stato.eta[i].y*dt;
				stato.forza[i].x=proposta.x;
				stato.forza[i].y=proposta.y;
				bool flag=closest(stato.r[i],&proposta,scheletro,N);
				stato.r[i]=proposta;	
				if((proposta.x>ottimizza.vertici[0].x)&&(proposta.x<ottimizza.vertici[1].x)&&(proposta.y<ottimizza.vertici[2].x)&&(proposta.y>ottimizza.vertici[0].y)) 	{		
						if (flag){
								stato.forza[i].x=(stato.forza[i].x-stato.r[i].x)/(dt*mobility);
								stato.forza[i].y=(stato.forza[i].y-stato.r[i].y)/(dt*mobility);}
													
					else {stato.forza[i].x=0.; stato.forza[i].y=0.;}

						}
		else {if(stato.r[i].x>half_size)	stato.r[i].x=stato.r[i].x-size;
				else if (stato.r[i].x<-half_size) stato.r[i].x=stato.r[i].x+size;
			if	(stato.r[i].y>half_size)	stato.r[i].y=stato.r[i].y-size;
				else if (stato.r[i].y<-half_size)	stato.r[i].y=stato.r[i].y+size;
			stato.forza[i].x=0.; stato.forza[i].y=0.;	}
/*					

					if (stato.r[i].x>ottimizza[0].x){
						if(stato.r[i].x<ottimizza[1].x){
							if(stato.r[i].y>ottimizza[0].y){
								if(stato.r[i].y<ottimizza[2].y{
									if (is_in(scheletro,stato.r[i],N)){
										stato.forza[i]=stato.r[i];
										stato.r[i]=segmento_vicino(scheletro.vertici,N,stato.r[i]);
										stato.forza[i].x=(stato.forza[i].x-stato.r[i].x)/(dt*mobility);
										stato.forza[i].y=(stato.forza[i].y-stato.r[i].y)/(dt*mobility);
													}
									else	{stato.forza[i].x=0.; stato.forza[i].y=0;}
												}
								else if(stato.r[i].y>half_size)	stato.r[i].y=stato.r[i].y-size;
									stato.forza[i].x=0.; stato.forza[i].y=0;
											}
							else if	(stato.r[i].y<-half_size) stato.r[i].y=stato.r[i].y+size;
									stato.forza[i].x=0.; stato.forza[i].y=0; 											}
					else if (stato.r[i].x<-half_size) stato.r[i].x=stato.r[i].x+size;	
						 						}}
				
					else if (stato.r[i].x<-half_size) stato.r[i].x=stato.r[i].x+size;
					if	(stato.r[i].y>half_size)	stato.r[i].y=stato.r[i].y-size;
					else if (stato.r[i].y<-half_size)	stato.r[i].y=stato.r[i].y+size;
						if (is_in(scheletro,stato.r[i],N)){ //printf("la particella %d è caduta dentro \n",i);
							stato.forza[i]=stato.r[i];
							stato.r[i]=segmento_vicino(scheletro.vertici,N,stato.r[i]);
							stato.forza[i].x=(stato.forza[i].x-stato.r[i].x)/(dt*mobility);
							stato.forza[i].y=(stato.forza[i].y-stato.r[i].y)/(dt*mobility);
}

						else {stato.forza[i].x=0.; stato.forza[i].y=0;}*/
						}
}

void ottimizza_geometria(geometria scheletro, geometria* ottimizza, int N){int i;
					(*ottimizza).vertici=(point*)malloc((4+1)*sizeof(point));
					double min_x,min_y,max_x,max_y;
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
void alloco_punti(geometria* scheletro,int *N){  int i=0; double a,b;
					FILE*f;f=fopen("input2.dat","r");
					while(fscanf(f,"%f" "%f",&a,&b)>0)i++;
					*N=i;rewind(f);
					(*scheletro).vertici=(point*)malloc((*N +1)*sizeof(point));
					for(i=0;i<*N;i++){
							fscanf(f,"%lf %lf",&(*scheletro).vertici[i].x,&(*scheletro).vertici[i].y);
							}
					fclose(f);
					(*scheletro).vertici[*N].x=(*scheletro).vertici[0].x;
					(*scheletro).vertici[*N].y=(*scheletro).vertici[0].y;//sto creando una geometria con un vertice in più che coincide con il primo		
					}

main(){
creo_cartelle();
clock_t t1 = clock();
int N;
geometria scheletro;
geometria ottimizza;
alloco_punti(&scheletro,&N);
ottimizza_geometria(scheletro,&ottimizza,N);
srand(10);
int t,i=0;
/*
for(i=0;i<4;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,ottimizza.vertici[i].x,ottimizza.vertici[i].y,ottimizza.vertici[i+1].x,ottimizza.vertici[i+1].y);}*/
configurazione stato;
inizializza(&stato,scheletro,N);
int numero_passi=(int)durata/dt;
int passi_salvataggio =(int)(tsalva/dt);
for(t=0;t<numero_passi;t++){
//for(t=0;t<8;t++){
			evolvi(stato,scheletro,N,ottimizza);
			   if((t% passi_salvataggio==0)&&(t>0)){stampa(stato,i);i++;}
			//stampa(stato,t);}
			}
point a,b,c;a.x=.25;a.y=-1.5; b.x=0.25;b.y=0;

c.x=b.x; c.y=b.y;
closest(a,&b,scheletro,N);


printf("tra il punto (%f,%f) e il punto (%f,%f) si trova (%f,%f)\n",a.x,a.y,c.x,c.y,b.x,b.y);
//stampa(stato,0);
 clock_t t2 = clock();
double time_sec = 
       (double)(t2-t1)/(double)(CLOCKS_PER_SEC); 
 
    printf("Time (sec): %lf\n",time_sec); 
}


