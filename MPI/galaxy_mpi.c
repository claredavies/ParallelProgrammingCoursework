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
/*
 *  In MPI programs, the main function for the program is run on every
 *  process that gets initialized when you start up this code using mpirun.
 */
float *real_rasc, *real_decl, *rand_rasc, *rand_decl;
float  pif;
long int MemoryAllocatedCPU = 0L;
double start,time,stop;

int main (int argc, char *argv[]) 
{
	// declaring functions
    int parseargs_readinput(int argc, char *argv[]);
    int get_index(float rasc_1, float decl_1, float rasc_2, float decl_2);

    start = MPI_Wtime();
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
  	int i,j;       // loop index
  	int limit = 5;
    int result_wanted = limit*limit;
    printf("   Result wanted %d \n", result_wanted);

  	MPI_Status status;   

  	MPI_Init (&argc, &argv);   // Initialization of MPI environment
  	MPI_Comm_size (MPI_COMM_WORLD, &total_proc); //  total number of processes running in parallel
  	MPI_Comm_rank (MPI_COMM_WORLD,&rank); //rank of the current process
  	
  	if (rank == MASTER)  {
  	   if ( parseargs_readinput(argc, argv) != 0 ) {
        printf("   Program stopped.\n");return(0);
      }
      printf("Input data read, now calculating histograms\n");
    }

    MPI_Bcast(real_rasc, limit, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(real_decl, limit, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(rand_rasc, limit, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(rand_decl, limit, MPI_INT, MASTER, MPI_COMM_WORLD);

    int range = limit/total_proc;
    int remainder = limit - (range*total_proc);
    int even_division_remainder = rank%2;

    int end,start;
   
    if(rank == MASTER) {
      start = rank*range;
      end = start + range;
    }

    else if(rank == (total_proc-1)){
      start = (rank*range)+1;
      end = (start-1) + range + remainder;
      limit = limit + 2;
    }
    
    else {
      start = (rank*range)+1;
      end = (start-1) + range;
    }

    printf("rank %d: start = %d, end = %d\n",rank,start,end);

  // for DR
    for(i = start; i <= end; ++i) 
    {
      for(j =0; j < (limit-1); ++j) {
        if( j == start && i == (start+1)) {

        }
        else {
               printf("    i currently:   %d  j currently:  %d \n",i,j);
            ++histogram_DR[get_index(real_rasc[i], real_decl[i], rand_rasc[j], rand_decl[j])];
        }
        }
    }

    // // // for RR
    // for(i = start; i <= end; ++i) {
    //   for(j =0; j < (limit-1); ++j) {   
    //        ++histogram_RR[get_index(real_rasc[i], real_decl[i], real_rasc[j], real_decl[j])];
    //   }
    // }

    // // for DD
    // for(i = start; i <= end; ++i) {
    //   for(j =0; j < (limit-1); ++j)  {  
    //     ++histogram_DD[get_index(rand_rasc[i], rand_decl[i], rand_rasc[j], rand_decl[j])];
    //   }
    // }
    
	MPI_Reduce(histogram_DD,histogram_DD_total,360,MPI_LONG,MPI_SUM,MASTER,MPI_COMM_WORLD);
	MPI_Reduce(histogram_DR,histogram_DR_total,360,MPI_LONG,MPI_SUM,MASTER,MPI_COMM_WORLD);
    MPI_Reduce(histogram_RR,histogram_RR_total,360,MPI_LONG,MPI_SUM,MASTER,MPI_COMM_WORLD);

    MPI_Finalize();  


	// Master process gets to here only when it has been able to gather from all processes
	if (rank == MASTER)  {	
     // check point: the sum of all historgram entries should be 10 000 000 000
    long int histsum = 0L;
    int correct_value=1;

    for ( int i = 0; i < 360; ++i ) histsum += histogram_DD_total[i];
    printf("   Histogram DD : sum = %ld\n",histsum);
    if ( histsum != result_wanted ) correct_value = 0;

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) 
      histsum += histogram_DR_total[i];
    printf("   Histogram DR : sum = %ld\n",histsum);
    if ( histsum != result_wanted ) correct_value = 0;

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_RR_total[i];
    printf("   Histogram RR : sum = %ld\n",histsum);
    if ( histsum != result_wanted ) correct_value = 0;

    if ( correct_value != 1 ) 
       {printf("   Histogram sums should be %d  Ending program prematurely\n",result_wanted);return(0);}

    printf("Omega values for the histograms:\n");
    float omega[360];
    for ( int i = 0; i < 360; ++i ) 
    {
        if ( histogram_RR_total[i] != 0L )
        {
           omega[i] = (histogram_DD_total[i] - 2L*histogram_DR_total[i] + histogram_RR_total[i])/((float)(histogram_RR_total[i]));
           if ( i < 10 ) 
            printf("   angle %.2f deg. -> %.2f deg. : %.3f\n", i*0.25, (i+1)*0.25, omega[i]);
        }
    }

    FILE *out_file = fopen(argv[3],"w");
    if ( out_file == NULL ) printf("   ERROR: Cannot open output file %s\n",argv[3]);
    else
       {
       for ( int i = 0; i < 360; ++i ) 
       {
           if ( histogram_RR_total[i] != 0L )
              fprintf(out_file,"%.2f  : %.3f\n", i*0.25, omega[i] ); 
       }
          
       fclose(out_file);
       printf("   Omega values written to file %s\n",argv[3]);
       }
       
    stop = MPI_Wtime();
    time = stop-start;	
	}

	// clean up memory
	if (rank == MASTER)  {
		free(real_rasc);  free(real_decl); free(rand_rasc); free(rand_decl);
	}
	
	// 9. Terminate MPI Environment and Processes
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
  
    if ( argc != 4 ) {
       printf("   Usage: galaxy real_data random_data output_file\n   All MPI processes will be killed\n");
       return(1);
    }
    if ( argc == 4 ){
       printf("   Running galaxy_openmp %s %s %s\n",argv[1], argv[2], argv[3]);

       real_data_file = fopen(argv[1],"r");
       if ( real_data_file == NULL )  {
          printf("   Usage: galaxy  real_data  random_data  output_file\n");
          printf("   ERROR: Cannot open real data file %s\n",argv[1]);
          return(1);
      }
       else
    	{
        fscanf(real_data_file,"%d",&Number_of_Galaxies);
        for ( int i = 0; i < 100000; ++i ) {
          float rasc, decl;
          if ( fscanf(real_data_file,"%f %f", &rasc, &decl ) != 2 ){
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
       if ( rand_data_file == NULL ) {
          printf("   Usage: galaxy  real_data  random_data  output_file\n");
          printf("   ERROR: Cannot open random data file %s\n",argv[2]);
          return(1);
        }
       else {
          fscanf(rand_data_file,"%d",&Number_of_Galaxies);
          for ( int i = 0; i < 100000; ++i ) {
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
      if ( out_file == NULL ) {
        printf("   Usage: galaxy  real_data  random_data  output_file\n");
        printf("   ERROR: Cannot open output file %s\n",argv[3]);
        return(1);
        }
      else fclose(out_file);
    }
	
	  return 0;
}