//Farhana Zaman Glory
//7856215
//COMP 7850
//Advances in parallel Computing
//Project Topic: Hybrid (MPI+OpenMP)Implementation of Yang's Algorithm for Longest Common Subsequence problem


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

// these are the block decomposition macros of MPI, these are needed for balanced loads of blocks in the processors
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p)) // finds low index
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)  //find high index
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n)) //find owner for given index

#define ALL_CHARS_SET 26

//function prototypes 

// this is the function which always returns 1, this is needed for increment scores 1 in score table if matching is found
short score(int x);

//for determining size of the block
int block_width(int strlen);

//for determining number of blocks per processor
int blocks_per_proc(int strlen, int block_width);

//for determining the last block size
int last_block_size(int strlen, int block_width);

// these two functions are needed for computing score table
char C(int index);
int letter_index(char letter);


//main
int main(int argc, char *argv) {
    int rank, nprocs;
    
    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;
    FILE *file;
    int input1_size = 0;
    char *input1, *input2;
    int input2_size = 0;
	
	  double start_time; // use these for timing
    double stop_time;

    file = fopen("Input4.txt", "r");
    if (file != NULL) {
        fscanf(file, "%d %d", &input1_size, &input2_size);

        input1 = (char *) calloc(input1_size + 1, sizeof(char));
        input2 = (char *) calloc(input2_size + 1, sizeof(char));

        fscanf(file, "%s %s", input1, input2);

    } else {
        printf("%s\n", "Unable to open file");
        exit(-1);
    }

    fclose(file);

    /*char C[ALL_CHARS_SET] = { 'A', 'C', 'G', 'T' };*/
    int *Pivot;
    
    //memory allocation for pivot table
    Pivot = (int *) calloc(ALL_CHARS_SET * (input2_size + 1), sizeof(int));


    int i, j, k, b;
    int width = input2_size + 1;
    int block_low = BLOCK_LOW(rank, nprocs, ALL_CHARS_SET);
    int block_high = BLOCK_HIGH(rank, nprocs, ALL_CHARS_SET);
    int block_size = BLOCK_SIZE(rank, nprocs, ALL_CHARS_SET);
    
    //pivot is making a table from input2 and with the C sequence, so we are taking block sizes according to input2
    int input2_width = block_width(input2_size + 1);
    int input2_per_proc = blocks_per_proc(input2_size + 1, input2_width);
    int last_input2_blck_size = last_block_size(input2_size + 1, input2_width);

    //the function forms a barrier, and no processes in the communicator can pass the barrier until all of them call the function.
    MPI_Barrier(MPI_COMM_WORLD);

    /* Compute Pivot table */
    #pragma omp for
    for (i = 0; i < ALL_CHARS_SET; i++) {
        int index = i * width;
        char c = C(i);
        for (j = 1; j < input2_size + 1; j++) {
            char b = input2[j - 1];
            if (b == c) {
                Pivot[index + j] = j;
            } else {
                Pivot[index + j] = Pivot[index + j - 1];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    block_low = BLOCK_LOW(rank, nprocs, input1_size + 1);
    block_high = BLOCK_HIGH(rank, nprocs, input1_size + 1);
    block_size = BLOCK_SIZE(rank, nprocs, input1_size + 1);
    int first = BLOCK_OWNER(0, nprocs, input1_size + 1); //index
    int last = BLOCK_OWNER(input1_size, nprocs, input1_size + 1); //index
    int row_size;
	  start_time = MPI_Wtime(); // start time

    if (rank == first) {
        row_size = block_size;
    } else {
        row_size = block_size + 1;
    }

    //memory allocation for score table
    int *Score_table = (int*) calloc(row_size * (input2_size + 1), sizeof(int));

    if (Score_table == NULL) {
        printf("You are out of memory. Please try again\n");
        exit(-1);
    }
    
    
    // computing the final score table
    #pragma omp for
    for(k = 0; k < input2_per_proc; k++) {
        int columns = input2_width;  // total column numbers in the score table
        int current_column = k * columns;  //buffer

        if (k + 1 == input2_per_proc) {
            columns = last_input2_blck_size; //holds the last block 
        }
        if (rank != first && block_size > 0) {
            MPI_Recv(Score_table + current_column, columns, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        }

        b = block_low == 0? 1 : block_low;
        //here i count row wise and j counts column wise
        for (i = 1; b < block_high + 1 && i < input1_size + 1; i++, b++) {
            int index = i * width;

            
            for (j = 0; j < columns; j++) {
                if (k != 0 || j != 0) {
                    char xi = input1[b - 1];
                    int column = current_column + j;
                    char yj = input2[column - 1];
                    int l_val = letter_index(xi);
                    int p_val = Pivot[l_val * width + column]; // here calculating p[c,j]

                    int t = (0 - p_val) < 0? 1 : 0;
                    int s;
                    // yang's algorithm of computing score table
                    if (p_val == 0) {
                        s = (0 - (Score_table[index - width + column] - t * Score_table[index - width])) < 0? 1 : 0;
                    } else {
                        s = (0 - (Score_table[index - width + column] - t * Score_table[index - width + p_val - 1])) < 0? 1 : 0;
                    }

                    if (xi == yj) score(i);
                    Score_table[index + column] = Score_table[index - width + column] + (t * (s ^ 1)); /* if value matches , then here it increases the value of the prev row's same indexed j's value with 1*/
                }
            }
        }

        if (block_size > 0 && block_high < input1_size) {
            int next = BLOCK_OWNER(block_high + 1, nprocs, input1_size + 1);  // rank of the destination process
            if (rank == first) {
                MPI_Send(&Score_table[block_high * width + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD); // root sending the blocks of input2 to the processors 
            } else {
                MPI_Send(&Score_table[block_size * width + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD);
            }
        }
    }

    //printing the LCS 
    i = row_size - 1;
    j = input2_size;
    char xi, yj;
    int last_cell;
    int lcs_len;
    int *lcs_reconstruct; //LCS reconstructing array

    // receiving LCS length and lcs from the processors
    if (rank != last) {
        MPI_Recv(&lcs_len, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        lcs_reconstruct = (int*) calloc(lcs_len + 1, sizeof(int));
        MPI_Recv(lcs_reconstruct, lcs_len + 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        j = lcs_reconstruct[0];
        last_cell = Score_table[i * width + j];
    } else {
        last_cell = Score_table[i * width + j];
        lcs_len = last_cell;
        lcs_reconstruct = (int*) calloc(lcs_len + 1, sizeof(int));
    }

    while(i > 0 && last_cell > 0) {
        int index = i * width + j;
        if (rank == first) {
            xi = input1[(block_low + i) - 1];
        } else {
            xi = input1[(block_low + i) - 2];
        }
        yj = input2[j - 1];
        if (xi == yj) {
            lcs_reconstruct[last_cell] = xi;
            last_cell--;
            i--;
            j--;
        } else if (Score_table[index - width] > Score_table[index - 1]) {
            i--; // take the character of prev row
        } else {
            j--; // take the character of prev column
        }
    }


    // sending lcs and length of LCS to the root
    if (rank != first)
    {
        lcs_reconstruct[0] = j;
        int prev = BLOCK_OWNER(block_low - 1, nprocs, input1_size + 1);
        MPI_Send(&lcs_len, 1, MPI_INT, prev, 0, MPI_COMM_WORLD);
        MPI_Send(lcs_reconstruct, lcs_len + 1, MPI_INT, prev, 0, MPI_COMM_WORLD);
    } 
    else 
    {
	    	stop_time = MPI_Wtime(); // stopping timer
        printf("%d\n", lcs_len);
        for (i = 1; i < lcs_len + 1; i++) 
        {
            printf("%c", lcs_reconstruct[i]);  // root is printing the whole LCS string
        }

        printf("\n");
	    	printf("Total execution time (sec) needed to solve LCS problem is: %f\n", stop_time - start_time);
    }

    free(Score_table);
    free(Pivot);
    free(input1);
    free(input2);
    free(lcs_reconstruct);

    MPI_Finalize();

    return 0;
}

// this function returns always score 1
short score(int x) {
    double dscore = 0;
	  dscore = dscore + pow(sin((double) x), 2) + pow(cos((double) x), 2);
    return (short) (dscore);
}

//for determining size of the block
int block_width(int strlen) {
    if (strlen < 100) {
        return strlen / 4;
    } else {
        return strlen / 50;
    }
}

//for determining number of blocks per processor
int blocks_per_proc(int strlen, int block_width) {
    float x = strlen / (block_width * 1.0);
    int result = x;

    if (x - result > 0.5) {
        result++;
    }
    return result;
}

//for determining the last block size
int last_block_size(int strlen, int block_width) {
    float x = strlen / (block_width * 1.0);
    int result = x;

    if (x - result > 0.5) {
        return strlen - (block_width * result);
    } else {
        return block_width + (strlen - (block_width * result));
    }
}

// these two functions are needed for computing score table
char C(int index) {
    return 'A' + index;
}

int letter_index(char letter) {
    return letter - 'A';
}
