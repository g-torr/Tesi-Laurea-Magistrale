/*
 * This program uses the host CURAND API to generate 100 
 * pseudorandom floats.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>

void stampa(float* x,int n){int i; for(i = 0; i < n; i++) {
        printf("%1.4f \n", x[i]);
    }
    printf("\n");}
void crea_hist(float*x, int n,float *hist,float h,float k,float a){int i,j,l; 
						

                                                

                                                

 

                                           for(i=0;i<h;i++){//fa aumentare la classe di freq

                                           l=0;*(hist+i)=0;

                                           

                                                       for (j=0;j<n;j++){ // fa girare la x

                                                   

                                                                             if((*(x+j)<(a+((i+1)*k)))&&(*(x+j)>(a+(i*k))))

                                                                             {l=l+1;

										*(hist+i)=l; 

                                                                                                                     } } 

                                                                                                                   }
									stampa(hist,h);}           

                            
int main(int argc, char *argv[])
{
    size_t n = 2000;
    size_t i;
    curandGenerator_t gen;
    float *devData, *hostData;

    /* Allocate n floats on host */
    hostData = (float *)calloc(n, sizeof(float));

    /* Allocate n floats on device */
   cudaMalloc((void **)&devData, n*sizeof(float));

    /* Create pseudo-random number generator */
   curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);//Ho letto però che nelle Monte Carlo in più dimensioni è conveniente usare quasirandom
//    curandCreateGenerator(&gen,CURAND_RNG_QUASI_SCRAMBLED_SOBOL32); // è un quasirandon famiglia SOBOL 32 scrabled
    /* Set seed */
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

    /* Generate n floats on device */
  curandGenerateUniform(gen, devData, n);
/*	float mean=0; float stddev=1;
    curandGenerateNormal(gen,devData,n,mean,stddev);
*/	

    /* Copy device memory to host */
    cudaMemcpy(hostData, devData, n * sizeof(float), cudaMemcpyDeviceToHost);

    /* Show result */
  stampa(hostData,n);

/*	float* hist; int h=20;hist=(float*)malloc(h*sizeof(float));
	float k=(float)6/h;
	crea_hist(hostData,n,hist,h,k,-3);*/
	
    /* Cleanup */
    curandDestroyGenerator(gen);
    cudaFree(devData);
    free(hostData);    
    return EXIT_SUCCESS;
}





