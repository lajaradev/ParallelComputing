#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

/* Function to create the two matrix dynamically */
unsigned char** createMatrix(int ROWS, int COLUMNS){ 
    
    unsigned char **matrix;
    
    matrix = (unsigned char**)malloc(ROWS * sizeof(unsigned char*)); // Reserve memory for rows
   
    for (int i = 0; i < ROWS; i++) {
        matrix[i] = (unsigned char*)malloc(COLUMNS * sizeof(unsigned char)); // Reserve memory for columns
    }

    return matrix;
}

/*  In this function we open the file that we indicate by argument.
    Read the file by lines
    Store it in an array */
void readFile (unsigned char **matrix, char *file_entry, int *ROWS, int *COLUMNS){ // READ MATRIX BY ROWS
    FILE *f = fopen(file_entry, "rb");
    if(f) { 
        for (int i = 0; i < *ROWS; i++) {
            int bytes_read = fread(matrix[i], sizeof(unsigned char), *COLUMNS, f);  
        }
        fclose(f);
    }
}

/* In this function we open a new file and write inside by rows the values of a matrix */
void writeFile(unsigned char **matrix, int ROWS, int COLUMNS, char *output){
    FILE *f = fopen(output, "wb");
    for (int i = 0; i < ROWS; i++) {
        size_t r1 = fwrite(matrix[i], sizeof(unsigned char), COLUMNS, f);
    }
    fclose(f);
}

void writeTXTMatrix(unsigned char **matrix, int ROWS, int COLUMNS){

    FILE *f2 = fopen("matrix.txt", "w");
    for(int i = 0; i < ROWS; i ++){
        for(int j = 0; j < COLUMNS; j ++){
            fprintf(f2, "%i   ", matrix[i][j]);  
        }
        fprintf(f2, "\n");
    }

}

/* This function is used to show us a few data from the matrix of the image that we read */
void showMatrixOrigianl(unsigned char **matrixOriginal){
    
    /*SHOW ORIGINAL MATRIX 10 ELEMENTS*/ 

    printf("\nORIGINAL MATRIX\n\n"); 
    for(int i = 0; i < 10; i++){  
        for(int j = 0; j < 10;j++){   
            printf("%i   ",matrixOriginal[i][j]);  
        }   
        printf("\n"); 
    } 
        
}

/* This function is used to show us a few data of the matrix that we are going to write */
void showMatrixFiltered(unsigned char **matrixFiltered){
    
   /*SHOW FILTERED MATRIX 10 ELEMENTS*/
       
    printf("\n\nFILTERED MATRIX \n\n"); 
    for(int i = 0; i < 10; i++){  
        for(int j = 0; j < 10;j++){   
            printf("%i   ",matrixFiltered[i][j]);  
        } 
        printf("\n");   
    }   
    
}

/* Function to swap values in a vector. */
void exchangeQuicksort(int *V, unsigned int l, unsigned int r){ 
    int aux;
    aux = V[l];
    V[l] = V[r];
    V[r] = aux;
}

/* Function to sort values in a vector recursively. */
void quicksort(int *V, unsigned int left, unsigned int right){

    unsigned int l, r, p;
    int pivot;

    p = (left + right) / 2;;

    if(p > 0){
        pivot = V[p];
        l = left;
        r = right;

        while(l <= r){
            while(V[l] < pivot)
                l++;

            while(V[r] > pivot)
                r--;

            if(l <= r){
                exchangeQuicksort(V, l, r);

                l++;
                r--;
            }
        }

        if(left < r)
            quicksort(V, left, r);
        if(l < right)
            quicksort(V, l, right);
    }
}

/*  Function that calculates, according to the number of processes, 
    the number of rows to process and the number of rows to send.*/
void rowsForEachProcess(int nproces, int row, int *rowsSend, int *rowsProcessed){

    int i, aux = 0;
    row = row - 2; // Subtract the edges

    for(i = 0; i < nproces; i ++){ // NUMBER OS ELEMENTS
        
        rowsProcessed[i] = row / nproces; // Share the work equally
        rowsSend[i] = rowsProcessed[i] + 2; // Need two more for calculates
        
    }

    if(row % nproces != 0){ // With the module obtain the rest of the division in case lines have been left without taking into account, and they are added to the last process
        
        rowsProcessed[nproces - 1] += row % nproces;
        rowsSend[nproces - 1] = rowsProcessed[nproces - 1] + 2;
    }
 
    /*for(i = 0; i < nproces; i ++){
        printf("Process [%i]\n", i);  
        printf("ELEMENTS SEND %i \n", rowsSend[i]);
        printf("ELEMENTS PROCESS %i \n", rowsProcessed[i]);  
        printf("\n");
    }*/

    
}

int main(int argc,char *argv[]){
    
    int row = atoi(argv[2]); // Fetch the value and cast it to int
    int col = atoi(argv[3]);

    int nproces, myrank, i, j; // Store the number of processes and the identifier of each process
    int *rowsSend; // Array that stores the number of rows to send
    int *rowsProcessed; // Array that stores the number of rows to process
    unsigned char **matrixOriginal; 
    unsigned char **matrixFiltered;
    int topLeft, topCenter, topRight, centerLeft, centerCenter, centerRight, downLeft, downCenter, downRight, average, median; 
    double starttime, finishtime; // Store in seconds the moment of start an end of the execution

    MPI_Status status;

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&nproces);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    rowsSend = (int*) malloc(nproces * sizeof(int));
	rowsProcessed = (int*) malloc(nproces*sizeof(int));

    if(strcmp(argv[4], "sobel") == 0){ // Send 2 rows more
        rowsForEachProcess(nproces, row + 2, rowsSend, rowsProcessed); // calculate rows for each process
    }
    else{
        rowsForEachProcess(nproces, row, rowsSend, rowsProcessed); // calculate rows for each process
    }
     
    if(myrank == 0){ // PROCESS 0

        matrixOriginal = createMatrix(row, col); // reserve memory for the original matrix

        FILE *f = fopen(argv[1], "rb"); // open file by argument
        // printf("Image %s \n", argv[1]);

        if(strcmp(argv[4], "sobel") == 0){
            
            matrixFiltered = createMatrix(row + 2, col + 2); // reserve memory for the filtered matrix
            
            if(f) { 
                for (int i = 1; i < row + 1; i++) { // Read file from row [1][1] to write in position [1][1] of the filtered matrix
                    for(int j = 1; j < col + 1; j ++){
                        int bytes_read = fread(&matrixFiltered[i][j], sizeof(unsigned char), 1, f);  
                    }
                }
                fclose(f);

                starttime = MPI_Wtime(); // moment of start       

                for(i = 1; i < row + 1; i ++){ // In all rows (symmetry)
                    matrixFiltered[i][0] = matrixFiltered[i][2]; // New column equal a position [][] original file
                    matrixFiltered[i][col + 1] = matrixFiltered[i][col - 1]; // The same but int last column
                }

                for(i = 1; i < col + 1; i ++){ // In all column (symmetry)
                    matrixFiltered[0][i] = matrixFiltered[2][i]; // New row equal a position [][] original file
                    matrixFiltered[row + 1][i] = matrixFiltered[row - 1][i]; // The same but int last column
                }

                matrixFiltered[0][0] = matrixFiltered[2][2]; // Give values to the corners
                matrixFiltered[0][col + 1] = matrixFiltered[2][col - 1];
                matrixFiltered[row + 1][0] = matrixFiltered[row - 1][2];
                matrixFiltered[row + 1][col + 1] = matrixFiltered[row - 1][col - 1];
            }
        }
        else{
            
            matrixFiltered = createMatrix(row, col); // reserve memory for the filtered matrix        

            if(f){ // if the file.raw exists ...
        
                readFile(matrixOriginal, argv[1], &row, &col);
                // showMatrixOrigianl(matrixOriginal);

                starttime = MPI_Wtime(); // moment of start  

                for(int i = 0; i < row; i++){ // COPY BORDERS
                    memcpy(matrixFiltered[i], matrixOriginal[i], row * sizeof(char) );}

                for(int j = 0; j < col; j++){ // COPY BORDERS
                    memcpy(matrixFiltered[j], matrixOriginal[j], col * sizeof(char) );}
            }       

       }
   
    }
    else{ // rest of process

        matrixOriginal = createMatrix(rowsSend[myrank], col);   // allocate memory for element blocks
          
        if(strcmp(argv[4], "sobel") == 0){
            matrixFiltered = createMatrix(rowsSend[myrank], col + 2);   
        }else{
            matrixFiltered = createMatrix(rowsSend[myrank], col);   
        }
    }

    if(myrank == 0){ // SEND DATA TO ALL PROCESS
    
        int rowStart = 0; // Firts row start at row 0
        int rowEnd = rowsSend[0] - 1 ; // Last row is all send - 1

    
        //printf("Process [%i]: elementsNeed %i, rowsProcessed %i \n", myrank, rowsSend[0], rowsProcessed[0]);
        //printf("Row Start %i -> Row End %i\n\n", rowStart, rowEnd);
    
        for(i = 1; i < nproces; i ++){ // Minus to process 0
            
        rowStart += rowsProcessed[i - 1]; // Row start became always add row processed - 1
        rowEnd = rowStart + (rowsSend[i] - 1); 

            for(j = rowStart; j <= rowEnd; j++){ //MPI SEND
                // Send (Data, columns, type, where, tag)
                if(strcmp(argv[4], "sobel") == 0){
                    MPI_Send(matrixFiltered[j], col + 2, MPI_UNSIGNED_CHAR, i, 5, MPI_COMM_WORLD); 
                }else{
                    MPI_Send(matrixOriginal[j], col, MPI_UNSIGNED_CHAR, i, 5, MPI_COMM_WORLD); 
                }                       
            }   

                       
        //printf("\nPROCESS [%i] SEND MATRIX ORIGINAL TO PROCESS [%i] \n", myrank, i);
        //printf("Rows Send: %i, Rows Processed %i \n", rowsSend[i], rowsProcessed[i]);
        //printf("SEND: Row Start %i -> Row End %i\n\n", rowStart, rowEnd);
        //showMatrixOrigianl(matrixOriginal);
        

        }
    }
    else{ // DATA RECIVE
      
        for(i = 0; i < rowsSend[myrank]; i ++){
            // Recive (Data, columns, type, where, tag)
            if(strcmp(argv[4], "sobel") == 0){
                MPI_Recv(matrixFiltered[i], col + 2, MPI_UNSIGNED_CHAR, 0, 5, MPI_COMM_WORLD, &status);
            }
            else{
                MPI_Recv(matrixOriginal[i], col, MPI_UNSIGNED_CHAR, 0, 5, MPI_COMM_WORLD, &status);
            }
            
        }
    
        // printf("\nPROCESS [%i] RECIVE MATRIX ORIGINAL TO PROCESS [%i]\n", myrank, 0);
        // printf("Rows Receives: %i \n", rowsSend[myrank]);
        // showMatrixOrigianl(matrixOriginal);
    
    }

	// ALL PROCESS START CALCULATIONS 

        //  printf("P [%i] -> Rows: %i\n", myrank, rowsProcessed[myrank]);
        //  printf("COLUMS: %i\n", col);
        
        if(strcmp(argv[4], "average") == 0){ // compare two strings, if they are the same then = 0
        
            for(i = 1; i <= rowsProcessed[myrank]; i ++){
                
                for(j = 1; j < col - 1; j ++){

                    // COPY THE VALUE THAT IS IN THE POSITIONS
                    topLeft = (int) matrixOriginal[i-1][j-1];
                    topCenter = (int) matrixOriginal[i-1][j];
                    topRight = (int) matrixOriginal[i-1][j+1];
                    centerLeft = (int) matrixOriginal[i][j-1];
                    centerCenter = (int) matrixOriginal[i][j];
                    centerRight = (int) matrixOriginal[i][j+1];
                    downLeft  = (int) matrixOriginal[i+1][j-1];
                    downCenter  = (int) matrixOriginal[i+1][j];
                    downRight  = (int) matrixOriginal[i+1][j+1];

                    average = (topLeft + topCenter + topRight + centerLeft + centerCenter + centerRight + downLeft + downCenter + downRight) / 9;         
                    matrixFiltered[i][j] = average;  // IN THE POSITION I COPY THE VALUE OF THE CALCULATED AVERAGE
                      
                }
                 
            }
            
        }

        else if(strcmp(argv[4], "median") == 0){ // compare two strings, if they are the same then = 0

            for(i = 1; i <= rowsProcessed[myrank]; i ++){
                
                for(j = 1; j < col - 1; j ++){

                    // COPY THE VALUE THAT IS IN THE POSITIONS
                    topLeft = (int) matrixOriginal[i-1][j-1];
                    topCenter = (int) matrixOriginal[i-1][j];
                    topRight = (int) matrixOriginal[i-1][j+1];
                    centerLeft = (int) matrixOriginal[i][j-1];
                    centerCenter = (int) matrixOriginal[i][j];
                    centerRight = (int) matrixOriginal[i][j+1];
                    downLeft  = (int) matrixOriginal[i+1][j-1];
                    downCenter  = (int) matrixOriginal[i+1][j];
                    downRight  = (int) matrixOriginal[i+1][j+1];

                    int V[9] = {topLeft, topCenter, topRight, centerLeft, centerCenter, centerRight, downLeft, downCenter, downRight};

                    quicksort(V,0,9); // ORDER ARRAY POSITION [0 - 8]
                                
                    median = V[4]; // TAKE CENTER VALUE

                    matrixFiltered[i][j] = median; // IN THE POSITION I COPY THE VALUE OF THE CALCULATED MEDIAN 
                }
                 
            }


        }

        else if(strcmp(argv[4], "sobel") == 0){  // compare two strings, if they are the same then = 0
             
            for(i = 1; i <= rowsProcessed[myrank]; i ++){
                
                for(j = 1; j < col + 1; j ++){

                    // value of c
                    int C = 
                    (matrixFiltered[i-1][j-1] * -1) +
                    (matrixFiltered[i-1][j]   *  0) +
                    (matrixFiltered[i-1][j+1] *  1) +
                    (matrixFiltered[i][j-1]   * -2) +
                    (matrixFiltered[i][j]     *  0) +
                    (matrixFiltered[i][j+1]   *  2) +
                    (matrixFiltered[i+1][j-1] * -1) +
                    (matrixFiltered[i+1][j]   *  0) +
                    (matrixFiltered[i+1][j+1] *  1);

                    // value of f
                    int F =
                    (matrixFiltered[i-1][j-1] * -1) +
                    (matrixFiltered[i-1][j]   * -2) +
                    (matrixFiltered[i-1][j+1] * -1) +
                    (matrixFiltered[i][j-1]   *  0) +
                    (matrixFiltered[i][j]     *  0) +
                    (matrixFiltered[i][j+1]   *  0) +
                    (matrixFiltered[i+1][j-1] *  1) +
                    (matrixFiltered[i+1][j]   *  2) +
                    (matrixFiltered[i+1][j+1] *  1);

                    int sobel = sqrt(pow(C,2) + pow(F,2)); // formula
                    matrixOriginal[i - 1][j - 1] = sobel; // IN THE POSITION I COPY THE VALUE OF THE CALCULATED SOBEL 
                     
                }
                 
            }
           
        }
        else{
            printf("UNDEFINED FUNCTION");
        }

        //showMatrixFiltered(matrixFiltered);

        if(myrank == 0){  // Process 0 receives data
            
            int rowStart = 0; // Firts row start at row 0
            for(i = 1; i < nproces; i ++){
                
                rowStart += rowsProcessed[i - 1]; // Row start became always add row processed - 1
                //printf("Rows Start Recv: %i\n", rowStart);
            
                for(j = 0; j < rowsProcessed[i]; j ++){ // Receives rows
                    if(strcmp(argv[4], "sobel") == 0){ 
                        MPI_Recv(matrixOriginal[rowStart + j] , col, MPI_UNSIGNED_CHAR, i, 7, MPI_COMM_WORLD, &status);
                    }
                    else{
                        MPI_Recv(matrixFiltered[rowStart + j] , col, MPI_UNSIGNED_CHAR, i, 7, MPI_COMM_WORLD, &status);
                    }
                   
                }
                
                //printf("P[%i], Send to me %i rows\n", i, rowsProcessed[i]);
                //showMatrixFiltered(matrixFiltered);
            }

        }
        else{  // All process send to process 0 the processed data                            
            
            if(strcmp(argv[4], "sobel") == 0){

                for(j = 0; j < rowsProcessed[myrank]; j ++){
                    
                    MPI_Send(&matrixOriginal[j][0], col, MPI_UNSIGNED_CHAR, 0, 7, MPI_COMM_WORLD);         
                            
                }
            }
            else{

                for(j = 1; j <= rowsProcessed[myrank]; j ++){
                   
                    MPI_Send(&matrixFiltered[j][1], col - 2, MPI_UNSIGNED_CHAR, 0, 7, MPI_COMM_WORLD);
                    
                }
                
                // Send from column 1 to the size to send, for this we use &[][]

                // printf(" I'm Process [%i] - Send to [%i] -> rowsProcessed %i\n", myrank, 0, rowsProcessed[myrank]);
            }
            
                         
        }

        if(myrank == 0){

            finishtime = MPI_Wtime(); // moment of end
            printf("Time: %f \n", finishtime - starttime);
            
            if(strcmp(argv[4], "sobel") == 0){ 
                writeFile(matrixOriginal, row, col, argv[5]);
                //writeTXTMatrix(matrixOriginal, row, col);
            }
            else{
                writeFile(matrixFiltered, row, col, argv[5]);
                writeTXTMatrix(matrixFiltered, row, col);
            }     
            
            FILE *f2 = fopen("data.txt", "w");
            
            fprintf(f2, "DATA -- PROCESSED IMAGE \n");
            fprintf(f2, "Input File Name: %s \n", argv[1]);
            fprintf(f2, "Image Size: Rows %i, Columns %i \n", row, col);
            fprintf(f2, "Output File Name: %s \n", argv[5]);
            fprintf(f2, "Number of processes: %i\n", nproces);
            fprintf(f2, "Time: %f \n", finishtime - starttime);
            
        }

    free(matrixOriginal);
    free(matrixFiltered);
    
    MPI_Finalize();
}
