#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>


int main(int argc, char ** argv) {

	int i, j, k, m, n;
	double *A, *X;
	double normA, normR, tmp;

	srand(0);

    	m = 20;
    	n = 10;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	A = (double *) malloc( m * m * sizeof(double));
 	for(j = 0; j < m; j++){
 		for(i = 0; i < m; i++){
			A[i+j*m] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
		}
	}

	X = (double *) malloc( m * n * sizeof(double));
 	for(k = 0; k < n; k++){
 		for(i = 0; i < m; i++){
			X[i+k*m] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
		}
	}


	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	for ( j = 0 ; j < m ; j++ ){
 		for(k = 0; k < n; k++) X[j+k*m] = X[j+k*m] / A[j+j*m];
		for ( i = j+1 ; i < m ; i++ ){
 			for(k = 0; k < n; k++) X[i+k*m] -= A[i+j*n] * X[j+k*m];
	}
	}

	clock_gettime(CLOCK_MONOTONIC, &end);

/*
//      
//      I need to write a check here
//
	normR = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			tmp = B[i][j];
			for (k = 0; k <= j; k++) {
				tmp -= A[i][k]*A[j][k];
			}
			normR += ( ( i == j ) ? 1.0e+00 : 2.00e+00 ) * tmp * tmp;
		}
	}
	normR = sqrt( normR );
*/

	double time_taken;
	time_taken = end.tv_sec - start.tv_sec;
	time_taken = time_taken + (end.tv_nsec - start.tv_nsec) * 1e-9;

	double perf = ((double) m)*((double) m)*((double) n)/time_taken*1e-9;

	printf("[ TRSM ] m = %5d, n = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n", m, n, time_taken, perf);

	free( X );
	free( A );

	return 0;
}
