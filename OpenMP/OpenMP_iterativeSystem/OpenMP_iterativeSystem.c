#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#define N 15000 // Size

/* Function to validate the number of arguments */
void argumentsOK(int argc){
     
     if(argc != 5){ // If argc is diferent of 3 argument --> error
        printf("./ite file iterations numThreads data.txt\n");
        exit(-1);
    }

}

double **createMatrix(){ // CREATE DYNAMIC MATRIX
    
    double **matrix;
    matrix = (double**)malloc(N * sizeof(double*)); // Reserve memory for rows
    for(int i = 0; i < N; i ++){
        matrix[i] = (double*)malloc(N * sizeof(double)); // Reserve memory for columns
    }
    return matrix;

}

double *createArray(){ // CREATE DYNAMIC ARRAY
    
    double *array;

    array = (double*)malloc(N * sizeof(double)); // Reserve memory for array
   
    return array;

}

void AddValueArrayIdentity(double *array){ // FILL ARRAY VALUE 1
 
    for(int i = 0; i < N; i++){
        array[i] = (double) 1; // All array = 1.0
    }

}

/*  In this function we open the file that we indicate by argument.
    Read the file by lines
    Store it in an matrix */
void readFile (double **matrix, char *file_entry){ // READ MATRIX BY ROWS
 
    FILE *f = fopen(file_entry, "rb");
    if(f) { 
        for (int i = 0; i < N; i++) {
            int bytes_read = fread(matrix[i], sizeof(double), N, f);  
        }
        fclose(f);
    }

}

/* In this function we open a new file and write inside by rows the values of a matrix */
void writeFile(double **matrix, char *file_entry){
    FILE *f = fopen(file_entry, "wb");
    for (int i = 0; i < N; i++) {
        size_t r1 = fwrite(matrix[i], sizeof(double), N, f);
    }
    fclose(f);
}

/* This function is used to show us a few data from the matrix of the image that we read */
void showMatrix(double **matrix){  // SHOW MATRIX 10 ELEMENTS

    printf("\n\nORIGINAL MATRIX\n\n"); 
    for(int i = 0; i < 10; i++){  
        for(int j = 0; j < 10;j++){   
            printf("%.2f   ",matrix[i][j]);  
        }   
        printf("\n"); } 
        
}

/* This function is used to show us a few data from the array of array */
void showArray(double *array){  // SHOW ARRAY 10 ELEMENTS

    printf("\nARRAY\n"); 
    for(int i = 0; i < 10; i++){  
        printf("%f   ",array[i]);  
        printf("\n"); } 
        
}

/* Assign random values to the new matrix */
void valueMatrix(double **matrix){ // GIVE VALUES TO THE MATRIX
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i == j){ // DIAGONAL VALUES = 1
                matrix[i][j] = 1.0;
            }else if(i > j){ // POSITIVE VALUE (BOTTOM) 
                matrix[i][j] = (double)50*(i+1)*(j+1)/((double)N*N*10000);;
            }else{ // NEGATIVE VALUE (TOP) 
                matrix[i][j] = (double)-50*(i+1)*(j+1)/((double)N*N*10000);
            }
        }
    }
}

int main(int argc, char *argv[]){ // gcc -fopenmp -o it OpenMP_iterativeSystem.c && ./it matrix.bin 5 3
    
    argumentsOK(argc); // Funtion control arguments

    int numThreads = atoi(argv[3]);
    int m = atoi(argv[2]); // Fetch the value and cast it to int
    int iam, np; // OpenMP variables
    int i, j, k, absHigher = 0; //  Declaration of auxiliary variables that are used to perform different operations
    int rowsProcessed;
    double accumulate = 0; 
    double absCompare = 0; // Is the value of maximum representable finite floating-point (double) number
    double absReal;

    double absMax = 0;
    double absMin = 0;

    double *x1 = createArray(); // Create a array
    double *x2 = createArray(); // Create a new array
    double **matrix = createMatrix(); // Create new matrix
    double startTime, finishTime; // Store in seconds the moment of start an end of the execution

    FILE *f = fopen(argv[1], "rb"); // open file by argument

    if(! f){ // if the file.bin no exists ...
        printf("Create new file...\n"); 
        valueMatrix(matrix); // Value random for a new matrix
        writeFile(matrix, argv[1]); // Write the matrix
    }
     
    readFile(matrix, argv[1]); // If file exists, read this file
    
    AddValueArrayIdentity(x2); // New array = 1

    FILE *f2 = fopen(argv[4], "w");
    fprintf(f2, "\nDATA -- ITERATIVE SYSTEM \n");
    fprintf(f2, "Matrix File Name: %s \n", argv[1]);
    fprintf(f2, "Output File Name: %s \n", argv[4]);
    fprintf(f2, "Number of processes: %i\n", numThreads);

    // ITERATIONS

    startTime = omp_get_wtime(); // Run time

    #pragma omp parallel num_threads(numThreads) shared(matrix, x1, x2, m, absMax, absMin, f2) private(iam, np, i, j, k, rowsProcessed, accumulate, absCompare) default(none)
    {

        np = omp_get_num_threads(); // Get the number os threads
        iam = omp_get_thread_num(); // Get the thread id

        rowsProcessed = N / np; // number of rows to process
        
        #pragma omp for schedule (static, rowsProcessed) // Distribute work statically using a fixed number os elements
        for(i = 0; i < N; i++){ // First iteration
            accumulate = 0;
            for(j = 0; j < N; j++){
                accumulate += matrix[i][j] * x2[i]; // The result vector is the matrix * the unit vector
            }
            x1[i] = accumulate;
            
        }

        #pragma omp barrier // Synchronize the threads
        
        for(k = 1; k < m; k ++){ // ALL THREADS DO ALL ITERATIONS
            
            absCompare = 0; // ALL threads must initialize to 0 at each iteration because it is private memory
            #pragma omp for schedule (static, rowsProcessed) reduction(max : absMax) reduction(min : absMin) // Get maximum and minimum
            for(i = 0; i < N; i ++){ // Access each element of the matrix to make the product for each element of the previous vector
                accumulate = 0;
                for(j = 0; j < N; j ++){
                    accumulate += matrix[i][j]*x1[j]; // The result vector is the matrix * the vector
                }
                x2[i] = accumulate;
                
                if(fabs(absCompare) < fabs(accumulate)){ // Calculate the maximum and minimum of each iteration
                   // printf("accumulate: %1.1f\n", accumulate);
                    absCompare = accumulate;

                    if(absMax < absCompare){
                        absMax = absCompare;
                    }

                    if(absMin > absCompare){
                        absMin = absCompare;
                    }
                }
            
            }
           
            #pragma omp master // Only thread 0 works
            {   
                if(fabs(absMax) < fabs(absMin)){ // if the |minimum| is greater than |maximum| is our new maximum
                    absMax = absMin;
                } 
                printf("Max in iteration %i: %1.1f\n", k, absMax);  // Print the results
                fprintf(f2, "   Max in iteration %i: %1.1f\n", k, absMax);
               // showArray(x2);
            }
           
            #pragma omp barrier  // Synchronize the threads
            //printf("Data %i: %1.1f (%i) ABSREAL: %1.1f\n", k, absCompare, absHigher, absReal);  // Print the results
            // printf("AbsMax: %1.1f -> fabs(absMax): %1.1f\n", absMax, fabs(absMax));
            
            #pragma omp for schedule (static, rowsProcessed) // Distribute work statically using a fixed number os elements
            for(i = 0; i < N; i ++){   
                
                if(fabs(absMax) > 25){ // Divide each position of the vector by the element with the largest absolute value.
                
                    x2[i] = x2[i] / absMax; // Fill the vector x1 with the results to use it in the next iteration           
                }

                x1[i] = x2[i]; // Copy it to the vector with which the iteration begins
                 
            }
           
           #pragma omp single // A single thread initializes variables to 0 and also synchronizes all other threads
           {
                absMax = 0;
                absMin = 0;
               
           }
           
        }

    }


    finishTime = omp_get_wtime(); // Ends time
    printf("Time: %f \n", finishTime - startTime);

  
    fprintf(f2, "Time: %f \n", finishTime - startTime);

    free(x1); // Free memory array
    free(x2); // Free memory array
   
    free(matrix);  // Free memory matrix

    return 0;
}