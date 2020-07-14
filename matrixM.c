#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <papi.h>

#define N 144
#define NUM_EVENTS 2
#define BLOCK_SIZE 8
#define THREADS 1

int main(int argc, char* argv){

	int i, j, k, i2, j2, k2, temp;
	int tmp[8];

	// Memory allocation for the matrices, depending on the given size
	int** A ;
	int** B ;
	int** C ;

	A = (int**) malloc(sizeof(int*) * N);
	B = (int**) malloc(sizeof(int*) * N);
	C = (int**) malloc(sizeof(int*) * N);

	for(int n = 0; n < N; n++){
		A[n] = (int*) malloc(sizeof(int) * N);
		B[n] = (int*) malloc(sizeof(int) * N);
		C[n] = (int*) malloc(sizeof(int) * N);
	}

	//srand (time(NULL)); 
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			A[i][j] = rand() % 9;
			B[i][j] = 1;
			C[i][j] = 0;
		}
	}

	// Get actual time
	double start_time = omp_get_wtime();

	// PAPI
	int retval;
	int events[NUM_EVENTS] = {PAPI_L2_TCM, PAPI_L2_TCA};
	int eventSet;
	long long values[NUM_EVENTS];

	retval = PAPI_library_init(PAPI_VER_CURRENT);
	if(retval != PAPI_VER_CURRENT)
		printf("Error 1!\n");

	int num_counters = PAPI_num_counters();
	printf("%d counters!\n");

	retval = PAPI_create_eventset(&eventSet);
	if(retval != PAPI_OK)
		printf("Error 2!\n");

	retval = PAPI_add_events(eventSet, events, NUM_EVENTS);
	if(retval != PAPI_OK)
		printf("Error 3!\n");

	retval = PAPI_start(eventSet);

	// matrixM_ijk(A, B, C);
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			for(k = 0; k < N; k++){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	// matrixM_ikj(A, B, C);
	/*for(i = 0; i < N; i++){
		int* rowA = A[i];
		int* rowC = C[i];
		for(k = 0; k < N; k++){
			int* rowB = B[k];
			int elemA = rowA[k];
			for(j = 0; j < N; j++){
				rowC[j] += elemA * rowB[j];
			}
		}
	}*/

	// matrixM_jki(A, B, C);
	/*for(j = 0; j < N; j++){
		for(k = 0; k < N; k++){
			int elemB = B[k][j];
			for(i = 0; i < N; i++){
				C[i][j] += A[i][k] * elemB;
			}
		}
	}*/

	// matrixM_jkiTranspose(A, B, C);
	/*for(j = 0; j < N; j++){
		for(k = 0; k < N; k++){
			int elemB = B[j][k];
			for(i = 0; i < N; i++){
				C[i][j] += A[k][i] * elemB;
			}
		}
	}*/

	// Block Optimization Code
	/*for(i2 = 0; i2 < N; i2 = i2 + BLOCK_SIZE){
		for(j2 = 0; j2 < N; j2 = j2 + BLOCK_SIZE){
			for(k2 = 0; k2 < N; k2 = k2 + BLOCK_SIZE){
				for(i = i2; i < i2 + BLOCK_SIZE; i++){
					for(j = j2; j < j2 + BLOCK_SIZE; j++){
						temp = 0;
						for(k = k2; k < k2 + BLOCK_SIZE; k++){
							temp += A[i][k] * B[j][k];
						}
						C[i][j] += temp;
					}
				}
			}
		}
	}*/

	// Vectorization Code
	/*for(i2 = 0; i2 < N; i2 = i2 + BLOCK_SIZE){
		for(j2 = 0; j2 < N; j2 = j2 + BLOCK_SIZE){
			for(k2 = 0; k2 < N; k2 = k2 + BLOCK_SIZE){
				for(i = i2; i < i2 + BLOCK_SIZE; i++){
					for(j = j2; j < j2 + BLOCK_SIZE; j++){
						for(k = 0; k < BLOCK_SIZE; k++){
							tmp[k] = A[i][k2 + k] * B[j][k2 + k];
						}
						C[i][j] += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
					}
				}
			}
		}
	}*/

	// Paralelization Code
	/*#pragma omp parallel for private (i2, j2, k2, i, j, k, tmp) num_threads(THREADS)
	for(i2 = 0; i2 < N; i2 = i2 + BLOCK_SIZE){
		for(j2 = 0; j2 < N; j2 = j2 + BLOCK_SIZE){
			for(k2 = 0; k2 < N; k2 = k2 + BLOCK_SIZE){
				for(i = i2; i < i2 + BLOCK_SIZE; i++){
					for(j = j2; j < j2 + BLOCK_SIZE; j++){
						for(k = 0; k < BLOCK_SIZE; k++){
							tmp[k] = A[i][k2 + k] * B[j][k2 + k];
						}
						C[i][j] += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
					}
				}
			}
		}
	}*/

	// PAPI
	retval = PAPI_stop(eventSet, values);	
	printf("PAPI_L2_TCM = %lld \n", values[0]);
	printf("PAPI_L2_TCA = %lld \n", values[1]); 

	// Get execution time
	double time = omp_get_wtime() - start_time;
	printf("Execution time: %f\n", time);

	// Result Matrix
	/*for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			printf("%d ", C[i][j]);
		}
		printf("\n");
	}*/

	return 1;
    
}
