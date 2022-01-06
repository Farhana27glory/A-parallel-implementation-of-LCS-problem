//Farhana Zaman Glory
//7856215
//COMP 7850
//Advances in parallel Computing
//Project Topic: Serial Implementation of Yang's Algorithm for Longest Common Subsequence problem


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALL_CHARS_SET 26

// function prototypes

// this is the function which always returns 1, this is needed for increment scores 1 in score table if matching is found
short score(int x);

// these two functions are needed for computing score table
char C(int index);
int letter_index(char letter);


int main(int argc, char const *argv) {
    FILE *file;

    int input1_size = 0;
    char *input1, *input2;
    int input2_size = 0;

    double start_time=0; // use these for timing
    double stop_time=0;


	//reading inputs from file
    file = fopen("Input0.txt", "r");
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

	//memory allocation for pivot table
    int *Pivot = (int *) calloc(ALL_CHARS_SET * (input2_size + 1), sizeof(int));

	//memory allocation for score table
    int *Score_table = (int*) calloc((input1_size + 1) * (input2_size + 1), sizeof(int));
    if (Score_table == NULL) {
        printf("You are out of memory, Please try again later\n");
        exit(-1);
    }

    int i, j;
    const int width = input2_size + 1;



    /* Compute Pivot table */

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

           // printf(" %d", Pivot[index + j]);
        }
    }

    printf("\n");

	// computing main score table

    for (i = 1; i < input1_size + 1; i++) {
        int index = i * width;
        for (j = 1; j < input2_size + 1; j++) {
            char xi = input1[i - 1];
            char yj = input2[j - 1];
            int l_val = letter_index(xi);
            int p_val = Pivot[l_val * width + j];

            int t = (0 - p_val) < 0? 1 : 0;
            int s = (0 - (Score_table[index - width + j] - t * Score_table[index - width + p_val - 1])) < 0? 1 : 0;

            if (xi == yj) score(i);
            {
                Score_table[index + j] = Score_table[index - width + j] + t * (s ^ 1);
                printf(" %d", Score_table[index + j]);

            }

        }
    }

    i = input1_size;
    j = input2_size;
    char xi, yj;
    int last_cell = Score_table[i * width + j];
    int lcs_len = last_cell;
    char lcs_reconstruct[lcs_len + 1];
    lcs_reconstruct[lcs_len] = '\0';


	//reconstructing LCS from the score table
    while(last_cell > 0) {
        int index = i * width + j;
        xi = input1[i - 1];
        yj = input2[j - 1];
        if (xi == yj) {
            lcs_reconstruct[last_cell - 1] = xi;
            last_cell--;
            i--;
            j--;
        } else if (Score_table[index - width] > Score_table[index - 1]) {
            i--;
        } else {
            j--;
        }
    }


    printf("\n%d\n%s\n", lcs_len, lcs_reconstruct);

    free(Score_table);
    free(Pivot);
    free(input1);
    free(input2);

    return 0;
}

// this is the function which always returns 1, this is needed for increment scores 1 in score table if matching is found
short score(int x) {
    double dscore = 0;
	  dscore = dscore + pow(sin((double) x), 2) + pow(cos((double) x), 2);
    return (short) (dscore);
}

// these two functions are needed for computing score table
char C(int index) {
    return 'A' + index;
}

int letter_index(char letter) {
    return letter - 'A';
}
