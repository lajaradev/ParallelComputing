#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#define N 15000 // Size

double **createMatrix(){ // CREATE DYNAMIC MATRIX WITH DEFINE
    
    double **matrix;

    matrix = (double**)malloc(N * sizeof(double*)); // Reserve memory for rows

    for(int i = 0; i < N; i ++){

        matrix[i] = (double*)malloc(N * sizeof(double)); // Reserve memory for columns

    }

    return matrix;

}

double **createMatrix2(int ROWS){ // CREATE DYNAMIC MATRIX WITH INT
    
    double **matrix;

    matrix = (double**)malloc(ROWS * sizeof(double*)); // Reserve memory for rows

    for(int i = 0; i < ROWS; i ++){

        matrix[i] = (double*)malloc(ROWS * sizeof(double)); // Reserve memory for columns

    }

    return matrix;

}

double *createArray(){ // CREATE DYNAMIC ARRAY WITH DEFINE
    
    double *array;

    array = (double*)malloc(N * sizeof(double)); // Reserve memory for array

    return array;

}

double *createArray2(int ROWS){ // CREATE DYNAMIC ARRAY WITH INT
    
    double *array;

    array = (double*)malloc(ROWS * sizeof(double)); // Reserve memory for array

    return array;

}

/* Assign random values to the new matrix */
void valueMatrix(double **matrix){ // GIVE VALUES TO THE MATRIX
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i == j){ // DIAGONAL VALUES = 1
                matrix[i][j] = 1.0;
            }else if(i > j){ // POSITIVE VALUE (BOTTOM) 
                 matrix[i][j] = (double)50*(i+1)*(j+1)/((double)N*N*10000); 
            }else{ // NEGATIVE VALUE (TOP) 
                matrix[i][j] = (double)-50*(i+1)*(j+1)/((double)N*N*10000);
            }
        }
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

void AddValueArrayIdentity(double *arrayIdentity){ // FILL ARRAY VALUE 1
 
    for(int i = 0; i < N; i++){
        arrayIdentity[i] = (double) 1; // All array = 1.0
    }

}

/* This function is used to show us a few data from the array of array */
void showArray(double *array){  // SHOW ARRAY 10 ELEMENTS

    printf("\nARRAY\n"); 
    for(int i = 0; i < 10; i++){  
        printf("%f   ",array[i]);  
        printf("\n"); } 
        
}

/*  Function that calculates, according to the number of processes, 
    the number of rows to process and the number of rows to send.*/
void rowsForEachProcess(int nproces, int *rowsSend, int *rowStart){

    int i;

    for(i = 0; i < nproces; i ++){ // NUMBER OS ELEMENTS
        
        rowsSend[i] = N / nproces; // Share the work equally
        
    }

    if(N % nproces != 0){ // With the module obtain the rest of the division in case lines have been left without taking into account, and they are added to the last process
        
        rowsSend[nproces - 1] += N % nproces;
       
    }

    rowStart[0] = 0; // firts process, rowstart is 0

    for (i = 1; i < nproces; i ++) {  // Calculation of the position from the start of sending to each process
        
        rowStart[i] = rowStart[i - 1] + rowsSend[i - 1];

    }
    
    /*for(i = 0; i < nproces; i ++){
        printf("Process [%i]\n", i);  
        printf("ELEMENTS SEND %i \n", rowsSend[i]);
        printf("ROW START %i \n", rowStart[i]);
        printf("\n");
    }*/

    
}

int main(int argc,char *argv[]){

    int nproces, myrank;
    int i, j, m, k;
    int absHigher;

    double accumulate; 
    double absCompare;
    double absCompareMax;
    double absCompareMin;
    double absReal;    
    double **matrix; // Create new matrix
    double *x1;
    double *x2;
    double starttime, finishtime; // Store in seconds the moment of start an end of the execution

    int *rowsSend; // Array that stores the number of rows to send
    int *rowStart;

    FILE *f2 = fopen(argv[3], "w");

    MPI_Status status;

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&nproces);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
    rowsSend = (int*) malloc(nproces * sizeof(int));
    rowStart = (int*) malloc(nproces * sizeof(int));

    rowsForEachProcess(nproces, rowsSend, rowStart); // calculate rows for each process

    if(myrank == 0){

        x1 = createArray(); // Create a array
        x2 = createArray();
        matrix = createMatrix();

    }
    else{ // Fair memory allocation for each process

        x2 = createArray(); // Create a array
        
        x1 = (double*)malloc(rowsSend[myrank] * sizeof(double));
        matrix = (double**)malloc(rowsSend[myrank] * sizeof(double*)); // Reserve memory for rows
   
        for (int i = 0; i < rowsSend[myrank]; i++) {
            matrix[i] = (double*)malloc(N * sizeof(double)); // Reserve memory for columns
        }

    }

    if(myrank == 0){
        
        FILE *f = fopen(argv[1], "rb"); // open file by argument
       
        if(! f){ // if the file.bin no exists ...
            printf("Create new file...\n"); 
            valueMatrix(matrix); // Value random for a new matrix
            writeFile(matrix, argv[1]); // Write the matrix
        }
        
        readFile(matrix, argv[1]); // If file exists, read this file
        //showMatrix(matrix);
    }
    
    
    if(myrank == 0){ // SEND DATA TO ALL PROCESS

        int rowStart = 0; // Firts row start at row 0
        int rowEnd = rowsSend[0] - 1 ; // Last row is all send - 1

        //printf("Process [%i]: rowsUse %i \n", myrank, rowsSend[0]);
        //printf("Row Start %i -> Row End %i\n\n", rowStart, rowEnd);
    
        for(i = 1; i < nproces; i ++){ // Minus to process 0
            
        rowStart += rowsSend[i - 1]; // Row start became always add row processed - 1
        rowEnd = rowStart + (rowsSend[i] - 1); 

            for(j = rowStart; j <= rowEnd; j++){ //MPI SEND

                // Send (Data, columns, type, where, tag)
                MPI_Send(matrix[j], N, MPI_DOUBLE, i, 5, MPI_COMM_WORLD); // Send the rows to process

            } 
                    
        //printf("\nPROCESS [%i] SEND MATRIX ORIGINAL TO PROCESS [%i] \n", myrank, i);
        //printf("Rows Send: %i \n", rowsSend[i]);
        //printf("SEND: Row Start %i -> Row End %i\n\n", rowStart, rowEnd);
        
        } 
    }
    else{

         for(i = 0; i < rowsSend[myrank]; i ++){

            // Recive (Data, columns, type, where, tag)
            MPI_Recv(matrix[i], N, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status); // Receive the rows to process

         }

        //printf("\nPROCESS [%i] RECIVE MATRIX ORIGINAL TO PROCESS [%i]\n", myrank, 0);
        //printf("Rows Receives: %i \n", rowsSend[myrank]);
        
        }
        
    // ALL PROCESS DOES ITERATIONS

    m = atoi(argv[2]); // Fetch the value and cast it to int 
  
    AddValueArrayIdentity(x2); // New array = 1
    
    starttime = MPI_Wtime(); // Moment of start   

    for(i = 0; i < rowsSend[myrank]; i++){ // First iteration
        accumulate = 0;
        for(j = 0; j < N; j++){
            accumulate += matrix[i][j] * x2[i]; // The result vector is the matrix * the unit vector
        }
        x1[i] = accumulate;
    }
     
    //showArray(x1);
    //printf("P[%i] -> rowsSend: %i, rowStart: %i\n", myrank, rowsSend[myrank], rowStart[myrank] );
    
    // MPI_Allgatherv(void *SendData, int SendDataNum, MPI_Datatype, void *ReceiveData, int ReceiveDataNum, int *Displacements, MPI_Datatype, MPI_Comm)
    // Collect the data from all x1 to enter them in x2, and we place them according to the displacement that is given by the row where they begin 
    
    MPI_Allgatherv(x1, rowsSend[myrank], MPI_DOUBLE, x2, rowsSend, rowStart ,MPI_DOUBLE, MPI_COMM_WORLD);
    
    //showArray(x2);
    
    for(k = 1; k < m; k ++){ // Loop iterations

        for(i = 0; i < rowsSend[myrank]; i ++){ // Access each element of the matrix to make the product for each element of the previous vector
            accumulate = 0;
            for(j = 0; j < N; j ++){
                accumulate += matrix[i][j]*x2[j]; // The result vector is the matrix * the vector
            }
            x1[i] = accumulate;
        }

        //printf("I'm P [%i] -> Iteration: %i\n", myrank, k);
        //showArray(x1);

        absCompare = 0;
        for(i = 0; i < rowsSend[myrank]; i ++){ // Get the absolute value and compare it with the auxiliary variable that is updated
            if(absCompare < fabs(x1[i])){
                absCompare = fabs(x1[i]);
                absReal = x1[i];
            }
        }
       
        //printf("P[%i] -> ABS: %1.1f ->  -> REAL: %1.1f\n",myrank, absCompare, absReal);
        
        //MPI_Allreduce(void *operando, void *result, int NumData, MPI_Datatype, MPI_Op, MPI_Coom)
        // Send the absolute maximums calculated in all processes
       
        MPI_Allreduce(&absReal, &absCompareMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&absReal, &absCompareMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        if(fabs(absCompareMax) > (fabs(absCompareMin))){
            absReal = absCompareMax;
        }
        else{
            absReal = absCompareMin;
        }
        
        //printf("P[%i]. absCompare: %1.1f - Max: %1.1f. - Min: %1.1f. - Real: %1.1f.\n",myrank, absCompare, absCompareMax,absCompareMin, absReal);
        //printf("P[%i]. absCompare: %.2f\n",myrank, absCompare);
        //printf("P[%i]. absComparemax: %.2f\n",myrank, absCompareMax);
        
        if(absCompare > 25.0){
            for(i = 0; i < rowsSend[myrank]; i ++){
                x1[i] = x1[i] / absReal; // Fill the vector x1 with the results to use it in the next iteration  
            }
        }
        
        //showArray(x1);

        //printf("P[%i] -> SendData: %i\n", myrank, rowsSend[myrank] );
        
        // MPI_Allgatherv(void *SendData, int SendDataNum, MPI_Datatype, void *ReceiveData, int ReceiveDataNum, int *Displacements, MPI_Datatype, MPI_Comm)
        // Collect the data from all x1 to enter them in x2, and we place them according to the displacement that is given by the row where they begin 

        MPI_Allgatherv(x1, rowsSend[myrank], MPI_DOUBLE, x2, rowsSend,rowStart, MPI_DOUBLE, MPI_COMM_WORLD);
        
        //showArray(x2);
        
        if(myrank == 0){ // Look for the position of the maximum
            fprintf(f2, "Max in iteration %i: %1.1f\n", k, absReal);  // Print the results
            
            // Show index of the absolute value
            /*absHigher = 0;
            for(i = 0; i < N; i ++){ // Get the absolute value and compare it with the auxiliary variable that is updated
                
                if(absReal == fabs(x2[i])){
                    absHigher = i; // Tells us the position of the maximum value 
                }
                else if(1.00 == fabs(x2[i])){
                    absHigher = i;
                }  
            }*/
            //printf("Max in iteration %i: %1.1f (%i)\n", k, absReal, absHigher);  // Print the results
        }


    }

    // END LOOP ITERATIONS

    if(myrank == 0){

        finishtime = MPI_Wtime(); // moment of end
        printf("Time: %f \n", finishtime - starttime);
        
        fprintf(f2, "\nDATA -- ITERATIVE SYSTEM \n");
        fprintf(f2, "Matrix File Name: %s \n", argv[1]);
        fprintf(f2, "Output File Name: %s \n", argv[3]);
        fprintf(f2, "Number of processes: %i\n", nproces);
        fprintf(f2, "Time: %f \n", finishtime - starttime);

       
    }








    free(matrix);  // Free memory matrix
    free(x1); // Free memory array
    free(x2); // Free memory array

    MPI_Finalize();
}
