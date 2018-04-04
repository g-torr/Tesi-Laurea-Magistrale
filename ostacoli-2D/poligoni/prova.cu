#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>

const int threads=1024;
const int blocks=30;
const int particelle=threads*blocks;
const int N=4;
const float size=6.0;


struct vertici{	float x;
		float y;
			};

struct geometria{	vertici*	punti;};

struct stato{		float*	x;
			float*	y;
			float*	eta_x;
			float* 	eta_y;};
void	stampa(stato configurazione){int i;FILE*f;
					f=fopen("output.dat","w");
		for(i=0;i<particelle;i++){fprintf(f,"%f		%f \n",configurazione.x[i],configurazione.y[i]);}
					fclose(f);}




// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
__global__ inline float isLeft( vertici P0, vertici P1, vertici P2 )
{
    return ( (P1.x - P0.x) * (P2.y - P0.y)
            - (P2.x -  P0.x) * (P1.y - P0.y) );
}


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
__global__ int wn_PnPoly( vertici P, vertici* V, int n ){
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
//************non viene usato in questa simulazione perchè con il winding number all edges that are totally above or totally below P get rejected after only two (2) inequality tests. However, currently popular implementations of the cn algorithm ) use at least three (3) inequality tests for each rejected edge. 

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int	cn_PnPoly( vertici P, vertici* V, int n ){
   					 int    cn = 0;    // the  crossing number counter

   					 // loop through all edges of the polygon
   				 for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
     					if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
        				|| ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
      	      								// compute  the actual edge-ray intersect x-coordinate
            								float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
            								if (P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
               									  ++cn;   // a valid crossing of y=P.y right of P.x
        									}
    							}
   	 return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}

__global__ bool is_present(geometria scheletro,vertici P){int wn=wn_PnPoly(P,scheletro.punti,N);
					//printf("il wn=%d \n",wn);
					if (wn>0) return true;
							else return false;
						}



__global__ void	inizializza(stato* configurazione,geometria scheletro){int id=threadIdx.x+ blockIdx.x*blockDim.x;
						vertici temp;bool flag=false;
						while(flag){
							temp.x=-1+size*(float)curand()/RAND_MAX;
							temp.y=-1+size*(float)curand()/RAND_MAX;

							if (!is_present(scheletro,temp)){flag=true;
										(*configurazione).x[id]=temp.x;
										(*configurazione).y[id]=temp.y;
											}
 									}					 									 															
					
							}




void alloco_vertici(geometria* scheletro){  int i; 
					(*scheletro).punti=(vertici*)malloc((N+1)*sizeof(vertici));
					for(i=0;i<N;i++){
					printf("inserisci ascissa e ordinata del %d° vertice\n",i);
					scanf("%f",&(*scheletro).punti[i].x);scanf("%f",&(*scheletro).punti[i].y);
							}
					(*scheletro).punti[N].x=(*scheletro).punti[0].x;
					(*scheletro).punti[N].y=(*scheletro).punti[0].y;//sto creando una geometria con un vertice in più che coincide con il primo
					}

main(){
geometria scheletro;
//numero di vertici dell'oggetto
alloco_vertici(&scheletro);
geometria scheletro_dev;
cudaMalloc((vertici**)&scheletro_dev.punti,(N+1)*sizeof(vertici));
cudaMemcpy(scheletro_dev.punti,scheletro.punti,(N+1)*sizeof(vertici),cudaMemcpyHostToDevice);

//for(int i=0;i<N;i++) {scheletro.punti[i].x=0;scheletro.punti[i].y=0;}

//cudaMemcpy(scheletro.punti,scheletro_dev.punti,(N+1)*sizeof(vertici),cudaMemcpyDeviceToHost);
srand(10);
//for(int i=0;i<N;i++) {printf("%f	%f \n",scheletro.punti[i].x,scheletro.punti[i].y);}


printf("finito di allocare\n");
stato configurazione;

cudaMalloc((float**)&configurazione.x,particelle*sizeof(float));

cudaMalloc((float**)&configurazione.y,particelle*sizeof(float));
//configurazione.y=(float*)malloc(particelle*sizeof(float));

inizializza<<<blocks,threads>>>(&configurazione,scheletro);
/*
vertici a,b,c;
a.x=0.;
a.y=0.;	b.x=1;b.y=1.;c.x=5.5;c.y=5.7;
if(isLeft(a,b,c)>0) printf("è a sinistra\n");
else printf("è a destra\n");*/
stampa(configurazione);
}


