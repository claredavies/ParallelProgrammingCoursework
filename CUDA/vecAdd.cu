#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int    NoofReal;
int    NoofRand;

int totaldegrees = 360;
int binsperdegree = 4;
float *h_real_rasc, *h_real_decl;
float *h_rand_rasc, *h_rand_decl;
__device__ int get_index(float rasc_1, float decl_1, float rasc_2, float decl_2);

// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(float *real_rasc, float *real_decl, float *rand_rasc, float *rand_decl, int *histogram_DD
    ,int *histogram_DR, int *histogram_RR, int n)
{
    // Get our global thread ID
    int idx = blockIdx.x*blockDim.x+threadIdx.x;

    // Make sure we do not go out of bounds
    if (idx < n) {
        printf("ID:   %d \n",idx);
        for(int j =0; j < n; ++j) {
                atomicAdd(&(histogram_DR[0]),(int)1);
        }
        //    printf("j:   %d \n",j);
        //    atomicAdd(&(histogram_DR[get_index(real_rasc[idx], real_decl[idx], rand_rasc[j], rand_decl[j])]),(int)1);
        //
     }
}

__device__ int get_index(float rasc_1, float decl_1, float rasc_2, float decl_2)
{
      float  pif;
      pif = acosf(-1.0f);
      float calcAngle = 0;
    	calcAngle = ((sin(decl_1)*sin(decl_2))+(cos(decl_1)*cos(decl_2))*cos(rasc_1-rasc_2));
   		if(calcAngle > 1.0)
          calcAngle = 1.0;
      else if(calcAngle < -1.0)
        calcAngle = -1.0;
      return (int) 4*acos(calcAngle)*180.0/pif;
}
 
int main( int argc, char* argv[] )
{
    int readdata(char *argv1, char *argv2);
    int getDevice(void);

    printf("here 1");
    if ( getDevice() != 0 ) return(-1);

    // Size of vectors
    int n = 10000;

    //Host output vector
    int *h_histogram_DR, *h_histogram_DD, *h_histogram_RR;

    FILE *outfil;
    if ( argc != 4 ) {printf("Usage: a.out real_data random_data output_data\n");return(-1);}

    // Device input vectors
    float *d_real_rasc; float *d_real_decl; float *d_rand_rasc; float *d_rand_decl;

    //Device output vector
    int *d_histogram_DR, *d_histogram_DD, *d_histogram_RR;
 
    // Size, in bytes, of each vector
    size_t bytes = n*sizeof(double);
    //size_t size_long_int = n*sizeof(long int);
    size_t size_int = n*sizeof(int);
 
    // Allocate memory for each vector on host
    h_real_rasc = (float*)malloc(bytes);
    h_real_decl = (float*)malloc(bytes);
    h_rand_rasc = (float*)malloc(bytes);
    h_rand_decl = (float*)malloc(bytes);

    h_histogram_DD = (int*)calloc((totaldegrees*binsperdegree+1L),size_int);
    h_histogram_DR = (int*)calloc((totaldegrees*binsperdegree+1L),size_int);
    h_histogram_RR = (int*)calloc((totaldegrees*binsperdegree+1L),size_int);

    // Allocate memory for each vector on GPU
    cudaMalloc(&d_real_rasc, bytes);
    cudaMalloc(&d_real_decl, bytes);
    cudaMalloc(&d_rand_rasc, bytes);
    cudaMalloc(&d_rand_decl, bytes);
    cudaMalloc(&d_histogram_DD, size_int );
    cudaMalloc(&d_histogram_DR, size_int );
    cudaMalloc(&d_histogram_RR, size_int );

    printf("here 2");

    if ( readdata(argv[1], argv[2]) != 0 ) return(-1);

    // Copy host vectors to device
    cudaMemcpy(d_real_rasc, h_real_rasc, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_real_decl, h_real_decl, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_rasc, h_rand_rasc, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_decl, h_rand_decl, bytes, cudaMemcpyHostToDevice);

    int blockSize, gridSize;
 
    // Number of threads in each thread block
    blockSize = 1024;
 
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)n/blockSize);
 
    // Execute the kernel
    vecAdd<<<gridSize, blockSize>>>(d_real_rasc, d_real_decl, d_rand_rasc, d_rand_decl, d_histogram_DD, d_histogram_DR
     , d_histogram_RR, n);
 
    // Copy array back to host
    cudaMemcpy( h_histogram_DD, d_histogram_DD, size_int, cudaMemcpyDeviceToHost );
    cudaMemcpy( h_histogram_DR, d_histogram_DR, size_int, cudaMemcpyDeviceToHost );
    cudaMemcpy( h_histogram_RR, d_histogram_RR, size_int, cudaMemcpyDeviceToHost );

    int histogramDRsum, histogramDDsum, histogramRRsum;

    histogramDRsum = 0;
    for ( int i = 0; i < binsperdegree*totaldegrees;++i ) {
         histogramDRsum += h_histogram_DR[i];
    }
    printf("   DR histogram sum = %ld\n",histogramDRsum);
    if ( histogramDRsum != 10000000000L )
     {   printf("   Incorrect histogram sum, exiting..\n");
         return(0);
     }

    histogramDDsum = 0L;
    for ( int i = 0; i < binsperdegree*totaldegrees;++i )
        histogramDDsum += h_histogram_DD[i];
    printf("   DD histogram sum = %ld\n",histogramDDsum);
    if ( histogramDDsum != 10000000000L ) {printf("   Incorrect histogram sum, exiting..\n");return(0);}

    histogramRRsum = 0L;
    for ( int i = 0; i < binsperdegree*totaldegrees;++i )
        histogramRRsum += h_histogram_RR[i];
    printf("   RR histogram sum = %ld\n",histogramRRsum);

    if ( histogramRRsum != 10000000000L ) {printf("   Incorrect histogram sum, exiting..\n");return(0);}


   printf("   Omega values:");

   outfil = fopen(argv[3],"w");
   if ( outfil == NULL ) {printf("Cannot open output file %s\n",argv[3]);return(-1);}

   fprintf(outfil,"bin start\tomega\t        hist_DD\t        hist_DR\t        hist_RR\n");
   for ( int i = 0; i < binsperdegree*totaldegrees; ++i )
       {
       if (h_histogram_RR[i] > 0 )
          {
          double omega =  (h_histogram_DD[i]-2*h_histogram_DR[i]+h_histogram_RR[i])/((double)(h_histogram_RR[i]));

          fprintf(outfil,"%6.3f\t%15lf\t%15ld\t%15ld\t%15ld\n",((float)i)/binsperdegree, omega,
             h_histogram_DD[i], h_histogram_DR[i], h_histogram_RR[i]);
          if ( i < 5 ) printf("   %6.4lf",omega);
          }
       else
          if ( i < 5 ) printf("         ");
       }

       printf("\n");

       fclose(outfil);

       printf("   Results written to file %s\n",argv[3]);

     // Release device memory
        cudaFree(d_real_rasc);
        cudaFree(d_real_decl);
        cudaFree(d_rand_rasc);
        cudaFree(d_rand_decl);
        cudaFree(d_histogram_DD);
        cudaFree(d_histogram_DR);
        cudaFree(d_histogram_RR);

        // Release host memory
        free(h_real_rasc);
        free(h_real_decl);
        free(h_rand_rasc);
        free(h_rand_decl);
        free(h_histogram_DD);
        free(h_histogram_DR);
        free(h_histogram_RR);

    return 0;
}

int readdata(char *argv1, char *argv2)
{
  int    i,linecount;
  char   inbuf[80];
  double ra, dec, dpi;
  FILE  *infil;

  printf("   Assuming data is in arc minutes!\n");
                          // phi   = ra/60.0 * dpi/180.0;
                          // theta = (90.0-dec/60.0)*dpi/180.0;
                          // otherwise use
                          // phi   = ra * dpi/180.0;
                          // theta = (90.0-dec)*dpi/180.0;

  dpi = acos(-1.0);
  infil = fopen(argv1,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv1);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv1, linecount-1);

  NoofReal = linecount-1;

  if ( NoofReal != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  h_real_rasc = (float *)calloc(NoofReal,sizeof(float));
  h_real_decl = (float *)calloc(NoofReal,sizeof(float));

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  i = 0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 )
         {
         printf("   Cannot read line %d in %s\n",i+1,argv1);
         fclose(infil);
         return(-1);
         }
      h_real_rasc[i] = (float)( ra/60.0*dpi/180.0);
      h_real_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal )
      {
      printf("   Cannot read %s correctly\n",argv1);
      return(-1);
      }

  infil = fopen(argv2,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv2);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv2, linecount-1);

  NoofRand = linecount-1;
  if ( NoofRand != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  h_rand_rasc = (float *)calloc(NoofRand,sizeof(float));
  h_rand_decl = (float *)calloc(NoofRand,sizeof(float));

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  i =0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 )
         {
         printf("   Cannot read line %d in %s\n",i+1,argv2);
         fclose(infil);
         return(-1);
         }
      h_rand_rasc[i] = (float)( ra/60.0*dpi/180.0);
      h_rand_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal )
      {
      printf("   Cannot read %s correctly\n",argv2);
      return(-1);
      }

  return(0);
}

int getDevice(void)
{

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  printf("   Found %d CUDA devices\n",deviceCount);
  if ( deviceCount < 0 || deviceCount > 128 ) return(-1);
  int device;
  cudaSetDevice(0);
  cudaGetDevice(&device);
  if ( device != 0 ) printf("   Unable to set device 0, using %d instead",device);
  else printf("   Using CUDA device %d\n\n", device);

  return(0);
}