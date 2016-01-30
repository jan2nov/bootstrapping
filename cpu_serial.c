#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#include "colors.h"

#define ALIGNSIZE 64
#define NR_ITER 10

int main( int argc, char** argv) 
{
 printf("\n  #######################\n   Blocked bootstrapping \n  #######################\n");


	//Define basic variables
	int nr_elms = 10240000;
	int nr_boots = 16384;	
	int elms_per_bin = 128;
	int nr_bins = nr_elms/elms_per_bin;

	//time variables
	double time_bins_start, time_bins_end;

 	float *data_array;
	float *bin_array;
	int *rand_array;
	float *boots_array;
	float boots_mean;
	float data_mean;
	float boots_sd;
	float partial;

	// Allocate
	data_array = (float*) _mm_malloc (nr_elms*sizeof(float),ALIGNSIZE);
	bin_array = (float*) _mm_malloc (nr_bins*sizeof(float),ALIGNSIZE);
	rand_array = (int*) _mm_malloc (nr_bins*nr_boots*sizeof(int),ALIGNSIZE);
	boots_array= (float*) _mm_malloc (nr_boots*sizeof(float),ALIGNSIZE);

	bin_array[0:nr_bins] = 0.0f;
	boots_array[0:nr_boots] = 0.0f;

	//random generator
	srand (time(NULL));	

	// Array notation set number to 0.0f
//	data_array[0:nr_elms] = 1.0f; //__sec_implicit_index(0);
	printf("\n\tCreate data array.\n");
	time_bins_start = omp_get_wtime();
	for (int i=0; i<nr_elms;i++)
		data_array[i] = floorf(1000*rand()/(float)RAND_MAX);	
	//	data_array[i] = (rand()/(float)RAND_MAX);	
	time_bins_end = omp_get_wtime() - time_bins_start;
	
	printf("\n\t\tThe time to create %i numbers: " COLOR_RED "%lf\n" COLOR_RESET, nr_elms,time_bins_end);

	//Random indices for boots - second time from CPU COSIT
	printf("\n\tGenerate random numbers to indicies.\n");
	time_bins_start = omp_get_wtime();
	for (int i=0; i<nr_bins*nr_boots;i++)
		rand_array[i] = (nr_bins-1)*(rand()/(float)RAND_MAX);
	time_bins_end = omp_get_wtime() - time_bins_start;

	printf("\n\t\tThe time to create %i random numbers for indices: " COLOR_RED "%lf" COLOR_RESET, nr_bins*nr_boots,time_bins_end);
	printf("\n\t\tFirst three elements in random indeces: %i %i %i",rand_array[0],rand_array[1],rand_array[2]);
	printf("\n\t\tThe max index is: %i. The min index is: %i\n",__sec_reduce_max(rand_array[0:nr_bins*nr_boots]),
								  __sec_reduce_min(rand_array[0:nr_bins*nr_boots]));

	//Create bins -- third time of CPU from COSIT
	printf("\n\tReduce the lenght to bins.\n");
	time_bins_end = 0;
	for(int j=0;j<NR_ITER;j++){
		time_bins_start = omp_get_wtime();
		for(int i = 0; i<nr_bins; i++)
			bin_array[i] = __sec_reduce_add(data_array[i*elms_per_bin:elms_per_bin]);
		time_bins_end += omp_get_wtime() - time_bins_start;
	}

	printf(COLOR_GREEN "\n\t\tResults in bin_array[%5i]: %f",0,bin_array[0]);
	printf("\n\t\tResults in bin_array[%5i]: %f",1,bin_array[1]);
	printf("\n\t\tResults in bin_array[%5i]: %f\n" COLOR_RESET,nr_bins-1,bin_array[nr_bins-1]);


	// print time of reduce to bins
	printf("\n\t\tThe time to create %i sums to bin_array: " COLOR_RED "%lf" COLOR_RESET, nr_bins,time_bins_end/NR_ITER);
	printf("\n\t\tNumber of operation %i. GFlops: " COLOR_RED "%4.3lf\n" COLOR_RESET, nr_elms, (nr_elms*NR_ITER)/(time_bins_end*1e9));

	// bootstrap_mean - 4th CPU COSIT
	printf("\n\tGenerate the #boots resample and sums.\n");
	time_bins_start = omp_get_wtime();
	for(int i = 0; i < nr_boots; i++)
		boots_array[i] = __sec_reduce_add(bin_array[rand_array[i*nr_bins:nr_bins]]);
	boots_array[0:nr_boots]	 /= nr_elms;
	time_bins_end = omp_get_wtime() - time_bins_start;

	printf(COLOR_GREEN "\n\t\tResults in boots_array[%5i]: %f", 0, boots_array[0]);
	printf("\n\t\tResults in boots_array[%5i]: %f", 1, boots_array[1]);
	printf("\n\t\tResults in boots_array[%5i]: %f\n" COLOR_RESET, nr_boots-1, boots_array[nr_boots-1]);

	// print time of bootstrap
	printf("\n\t\tThe time to create %i boots and reduce_add: " COLOR_RED "%lf" COLOR_RESET, nr_boots, time_bins_end);
	printf("\n\t\tNumber of operation %i. GFlops: " COLOR_RED "%4.3lf\n" COLOR_RESET, nr_bins*nr_boots, (nr_bins*nr_boots)/(time_bins_end*1e9));

	//Compute the mean of bootstraps and statistics
	printf("\n\tCompute mean of data and boots.\n");
	time_bins_start = omp_get_wtime();
		boots_mean = __sec_reduce_add(boots_array[0:nr_boots]);
		boots_mean /= nr_boots;
	time_bins_end = omp_get_wtime() - time_bins_start;
	printf("\n\t\tThe time to create boots mean: " COLOR_RED "%lf" COLOR_RESET, time_bins_end);

	time_bins_start = omp_get_wtime();
		data_mean = __sec_reduce_add(data_array[0:nr_elms]);
		data_mean /= nr_elms;
	time_bins_end = omp_get_wtime() - time_bins_start;
	printf("\n\t\tThe time to create data_mean: " COLOR_RED "%lf" COLOR_RESET, time_bins_end);

	time_bins_start = omp_get_wtime();
		for(int i = 0; i < nr_boots; i++)	
			partial += (boots_array[i] - boots_mean)*(boots_array[i] - boots_mean);
		boots_sd = sqrt(partial/(nr_boots - 1));
	time_bins_end = omp_get_wtime() - time_bins_start;
	printf("\n\t\tThe time to create standard deviation: " COLOR_RED "%lf\n" COLOR_RESET, time_bins_end);

	printf("\n\t\tMean of the boots: " COLOR_RED "%lf" COLOR_RESET, boots_mean);
	printf("\n\t\tDeviation of  boots: " COLOR_RED "%lf" COLOR_RESET, boots_sd);
	printf("\n\t\tMean of the data: " COLOR_RED "%lf\n" COLOR_RESET, data_mean);


/*	time_bins_start = omp_get_wtime();
	for (int i = 0; i < nr_boots; i++)
		for(int j = 0; j < nr_bins; j++)	
			boots_array[i] += bin_array[rand_array[i*nr_bins + j]];	
	time_bins_end = omp_get_wtime() - time_bins_start;
	printf("\n\t\tThe time to create %i boots and reduce_add: " COLOR_RED "%lf" COLOR_RESET, nr_boots, time_bins_end);
*/

	//Clean-up
	_mm_free(data_array);
	_mm_free(bin_array);
	_mm_free(rand_array);
	_mm_free(boots_array);

 
 printf("\n  ################### \n   That's All folks!\n  ###################\n");
}
