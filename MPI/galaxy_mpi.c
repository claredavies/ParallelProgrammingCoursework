/*
 *  Prerequisties:
 *     This code runs using an MPI library, either OpenMPI or MPICH2.
 *     These libraries can be installed in either a cluster of computers
 *     or a multicore machine.
 *     
 *  How to compile:
 *     mpicc -o vec-add VA-MPI-simple.c
 *
 *  How to execute:
 *     mpirun -np 2 ./vec-add
 *
 *     Note that this executes the code on 2 processes, using the -np command line flag.
 *     See ideas for further exploration of MPI using this code at the end of this file.
 */


#include "mpi.h"      // must have a system with an MPI library
#include <stdio.h>    //printf
#include <stdlib.h>   //malloc
#include <math.h>
#include <sys/time.h>
/*
 * Definitions
 */
#define MASTER 0         //One process will take care of initialization
#define ARRAY_SIZE 10     //Size of arrays that will be added together.

/*
 *  In MPI programs, the main function for the program is run on every
 *  process that gets initialized when you start up this code using mpirun.
 */
float *real_rasc, *real_decl, *rand_rasc, *rand_decl;
float *real_rasc_d, *real_decl_d, *rand_rasc_d, *rand_decl_d;
float  pif;
long int MemoryAllocatedCPU = 0L;

int main (int argc, char *argv[]) 
{
	    // declaring functions
    int parseargs_readinput(int argc, char *argv[]);
    int get_index(float rasc_1, float decl_1, float rasc_2, float decl_2);

    pif = acosf(-1.0f);

	long int histogram_DD[360] = {0L};
    long int histogram_DR[360] = {0L};
    long int histogram_RR[360] = {0L};

	long int histogram_DD_total[360] = {0L};
	long int histogram_DR_total[360] = {0L};
	long int histogram_RR_total[360] = {0L};

    real_rasc        = (float *)calloc(100000L, sizeof(float));
    real_decl        = (float *)calloc(100000L, sizeof(float));
    rand_rasc        = (float *)calloc(100000L, sizeof(float));
    rand_decl        = (float *)calloc(100000L, sizeof(float));

    MemoryAllocatedCPU += 10L*100000L*sizeof(float);

	
	int total_proc;	 // total nuber of processes	
	int rank;        // rank of each process
	int n_per_proc;	// elements per process	
	int n = ARRAY_SIZE;   // number of array elements
	int i,j;       // loop index
	int limit = 10;
	int count = 0;

	MPI_Status status;   // not used in this arguably poor example
	                     // that is devoid of error checking.

	// 1. Initialization of MPI environment
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &total_proc);
	// 2. Now you know the total number of processes running in parallel
	MPI_Comm_rank (MPI_COMM_WORLD,&rank);
	// 3. Now you know the rank of the current process
	
	// 4. We choose process rank 0 to be the root, or master,
	// which will be used to  initialize the full arrays.
	if (rank == MASTER)  {
		real_rasc        = (float *)calloc(100000L, sizeof(float));
	    real_decl        = (float *)calloc(100000L, sizeof(float));
	    rand_rasc        = (float *)calloc(100000L, sizeof(float));
	    rand_decl        = (float *)calloc(100000L, sizeof(float));

	    if ( parseargs_readinput(argc, argv) != 0 ) {
	    	printf("   Program stopped.\n");return(0);
	    }
    	printf("Input data read, now calculating histograms\n");
	}
	n_per_proc = n/total_proc;
	
	// 5. Initialize my smaller subsections of the larger array
	real_rasc_d = (float *) malloc(sizeof(float)*n_per_proc);
	real_decl_d = (float *) malloc(sizeof(float)*n_per_proc);
	rand_rasc_d = (float *) malloc(sizeof(float)*n_per_proc);
	rand_decl_d = (float *) malloc(sizeof(float)*n_per_proc);
	
	// 6.
	//scattering array a from MASTER node out to the other nodes
	MPI_Scatter(real_rasc, n_per_proc, MPI_INT, real_rasc_d, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD); 
	MPI_Scatter(real_decl, n_per_proc, MPI_INT, real_decl_d, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD); 
	MPI_Scatter(rand_rasc, n_per_proc, MPI_INT, rand_rasc_d, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD); 
	MPI_Scatter(rand_decl, n_per_proc, MPI_INT, rand_decl_d, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD); 

	// for DR
    for(i = 0; i < n_per_proc; ++i) 
    {
      for(j =0; j < n_per_proc; ++j) 
      {
        ++histogram_DR[get_index(real_rasc[j], real_decl[j], rand_rasc[j], rand_decl[j])];
      }
    }

    // // for RR
    for(i = 0; i < n_per_proc; ++i) {
      for(j = i+1; j < n_per_proc; ++j) {    
        histogram_RR[get_index(rand_rasc[j], rand_decl[j], rand_rasc[j], rand_decl[j])] += 2;
      }
    }
    histogram_RR[0] += limit;

    // // for DD
    for(i = 0; i < n_per_proc; ++i) {
      for(j = i+1; j < n_per_proc; ++j) {  
      	count++;
        histogram_DD[get_index(real_rasc[j],real_decl[j],real_rasc[j],real_decl[j])] += 2;
      }
    }
    histogram_DD[0] += limit;
		
	// 8. MASTER node gathering array c from the workers
	MPI_Gather(histogram_DD_total, n_per_proc, MPI_INT, histogram_DD, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Gather(histogram_DR_total, n_per_proc, MPI_INT, histogram_DR, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Gather(histogram_RR_total, n_per_proc, MPI_INT, histogram_RR, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);

/////////////////////// all concurrent processes are finished once they all communicate
/////////////////////// data back to the master via the gather function.

	// Master process gets to here only when it has been able to gather from all processes
	if (rank == MASTER)  {			
		// sanity check the result  (a test we would eventually leave out)
		// int good = 1;
		// for(i=0;i<n;i++) {
		// 	//printf ("%d ", c[i]);
		// 	if (c[i] != a[i] + b[i]) {
		// 		printf("problem at index %lld\n", i);
		// 		good = 0;
		// 		break;
		// 	}
		// }
		// if (good) {
		// 	printf ("Values correct!\n");
		// }
		
	}

	// clean up memory
	if (rank == MASTER)  {
		free(real_rasc);  free(real_decl); free(rand_rasc); free(rand_decl);
	}
	free(real_rasc_d);  free(real_decl_d); free(rand_rasc_d); free(rand_decl_d);
	
	// 9. Terminate MPI Environment and Processes
	MPI_Finalize();  


	    return(0);
}

	int get_index(float rasc_1, float decl_1, float rasc_2, float decl_2)
    {
        float calcAngle = 0;
    	calcAngle = ((sin(decl_1)*sin(decl_2))+(cos(decl_1)*cos(decl_2))*cos(rasc_1-rasc_2));
   		if(calcAngle > 1.0)
          calcAngle = 1.0;
        else if(calcAngle < -1.0)
          calcAngle = -1.0;
        return (int) 4*acos(calcAngle)*180.0/pif;
    }


int parseargs_readinput(int argc, char *argv[])
    {
    FILE *real_data_file, *rand_data_file, *out_file;
    float arcmin2rad = 1.0f/60.0f/180.0f*pif;
    int Number_of_Galaxies;
  
    if ( argc != 4 ) 
       {
       printf("   Usage: galaxy real_data random_data output_file\n   All MPI processes will be killed\n");
       return(1);
       }
    if ( argc == 4 )
       {
       printf("   Running galaxy_openmp %s %s %s\n",argv[1], argv[2], argv[3]);

       real_data_file = fopen(argv[1],"r");
       if ( real_data_file == NULL ) 
          {
          printf("   Usage: galaxy  real_data  random_data  output_file\n");
          printf("   ERROR: Cannot open real data file %s\n",argv[1]);
          return(1);
          }
       else
    	{
          fscanf(real_data_file,"%d",&Number_of_Galaxies);
          for ( int i = 0; i < 100000; ++i ) 
              {
              float rasc, decl;
        if ( fscanf(real_data_file,"%f %f", &rasc, &decl ) != 2 )
           {
                 printf("   ERROR: Cannot read line %d in real data file %s\n",i+1,argv[1]);
                 fclose(real_data_file);
           return(1);
           }
        real_rasc[i] = rasc*arcmin2rad;
        real_decl[i] = decl*arcmin2rad;

        }
           fclose(real_data_file);
     printf("   Successfully read 100000 lines from %s\n",argv[1]);
     }

       rand_data_file = fopen(argv[2],"r");
       if ( rand_data_file == NULL ) 
          {
          printf("   Usage: galaxy  real_data  random_data  output_file\n");
          printf("   ERROR: Cannot open random data file %s\n",argv[2]);
          return(1);
          }
       else 
    {
          fscanf(rand_data_file,"%d",&Number_of_Galaxies);
          for ( int i = 0; i < 100000; ++i ) 
              {
              float rasc, decl;
        if ( fscanf(rand_data_file,"%f %f", &rasc, &decl ) != 2 )
           {
                 printf("   ERROR: Cannot read line %d in real data file %s\n",i+1,argv[2]);
                 fclose(rand_data_file);
           return(1);
           }
        rand_rasc[i] = rasc*arcmin2rad;
        rand_decl[i] = decl*arcmin2rad;
        }
          fclose(rand_data_file);
    printf("   Successfully read 100000 lines from %s\n",argv[2]);
    }
       out_file = fopen(argv[3],"w");
       if ( out_file == NULL ) 
          {
          printf("   Usage: galaxy  real_data  random_data  output_file\n");
          printf("   ERROR: Cannot open output file %s\n",argv[3]);
          return(1);
          }
       else fclose(out_file);
       }
	
	return 0;
}