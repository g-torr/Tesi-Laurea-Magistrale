#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>

const int particelle=1000;
//const int N=4;
const float size=6.0;




struct point{	float x;
		float y;
			};

struct geometria{	point*	vertici;};

struct stato{		float*	x;
			float*	y;
			float*	eta_x;
			float* 	eta_y;};
void	stampa(stato configurazione){int i;FILE*f;
					f=fopen("output.dat","w");
		for(i=0;i<particelle;i++){fprintf(f,"%f		%f \n",configurazione.x[i],configurazione.y[i]);}
					fclose(f);}




float segmento_vicino(point * V, int N, point P,point * closest){ //P è il punto di cui viene calcolata la distanza rispetto al segmento di estremi P1-P2 ; t= indice  
	int i=0;float* t;float d;
	 t=(float*)malloc(N*sizeof(float));
	t[0] = ((P.x-V[0].x)*(V[0+1].x-V[0].x)+(P.y-V[0].y)*(V[0+1].y-V[0].y))/
								((V[0+1].x-V[0].x)*(V[0+1].x-V[0].x)+(V[0+1].y-V[0].y)*(V[0+1].y-V[0].y));

	 if(t[0]<0.0){t[0]=0.0;}
 	 if(t[0]>1.0){t[0]=1.0;} 


    (*closest).x = V[0].x+ (V[0+1].x-V[0].x)*t[0]; 
    (*closest).y = V[0].y+ (V[0+1].y-V[0].y)*t[0];  
	d=(P.x- (*closest).x)*(P.x- (*closest).x)+(P.y- (*closest).y)*(P.y- (*closest).y);

	point temp;
	float d_temp;
 	for(i=1;i<N;i++){
    t[i] = ((P.x-V[i].x)*(V[i+1].x-V[i].x)+(P.y-V[i].y)*(V[i+1].y-V[i].y))/
								((V[i+1].x-V[i].x)*(V[i+1].x-V[i].x)+(V[i+1].y-V[i].y)*(V[i+1].y-V[i].y));

    
    if(t[i]<0.0){t[i]=0.0;}
    if(t[i]>1.0){t[i]=1.0;} 


    temp.x = V[i].x+ (V[i+1].x-V[i].x)*t[i]; 
    temp.y =V[i].y+ (V[i+1].y-V[i].y)*t[i];  
	d_temp=(P.x-temp.x)*(P.x-temp.x)+(P.y-temp.y)*(P.y-temp.y);
			if(d_temp<d){ (*closest).x=temp.x; (*closest).y=temp.y;d=d_temp;}
				}
    
    return d; }

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
inline float	isLeft( point P0, point P1, point P2 ){
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




bool is_present(geometria scheletro,point P,int N){int wn=wn_PnPoly(P,scheletro.vertici,N);
					//printf("il wn=%d \n",wn);
					if (wn>0) return true;
							else return false;
						}



void	inizializza(stato* configurazione,geometria scheletro,int N){int i=0;point temp;		
						while(i<particelle){
							temp.x=-1+size*(float)rand()/RAND_MAX;
							temp.y=-1+size*(float)rand()/RAND_MAX;
							if (!is_present(scheletro,temp,N)){
										(*configurazione).x[i]=temp.x;(*configurazione).y[i]=temp.y;
											i++;}
 										}					 									 															
					
											}





void alloco_punti(geometria* scheletro,int *N){  int i=0; float a,b;
					FILE*f;f=fopen("input","r");
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
int N;
geometria scheletro;
//numero di punti dell'oggetto
alloco_punti(&scheletro,&N);
//crea_body(&scheletro);
srand(10);
int i;
/*
for(i=0;i<N;i++){
					printf("gli estremi del segmento %d° vertice sono (%f,%f) 	(%f,%f)\n",i,scheletro.lati[i].x0,scheletro.lati[i].y0,scheletro.lati[i].x1,scheletro.lati[i].y1);}*/

stato configurazione;

configurazione.x=(float*)malloc(particelle*sizeof(float));
configurazione.y=(float*)malloc(particelle*sizeof(float));

inizializza(&configurazione,scheletro,N);
point a,b;a.x=4.;a.y=.5;
float d=segmento_vicino(scheletro.vertici,N,a,&b);
printf("il punto (%f,%f) si trova a distanza%f dal punto (%f,%f)\n",a.x,a.y,sqrt(d),b.x,b.y);
stampa(configurazione);
}


