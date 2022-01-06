//Farhana Zaman Glory
//7856215
//COMP 7850
//Advances in Parallel Computing
//Serial Algorithm Implementation of Longest Common Subsequence Problem


#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <time.h>

//macros
#define max(x,y) ((x)>(y)?(x):(y))


//global variables
char *string_A;
char *string_B;
char *unique_chars_C; //unique alphabets
int c_len;
short **DP_Results;


//function prototypes
void print_matrix(int **x, int row, int col);
short lcs(short **DP, char *A, char *B, int m, int n);
char* lcs_reconstruct(short **DP_Results, char *string_A, char *string_B);



void print_matrix(int **x, int row, int col)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            printf("%d ",x[i][j]);
        }
        printf("\n");
    }
}


short lcs(short **DP, char *A, char *B, int m, int n)
{

    for(int i=1;i<(m+1);i++)
    {
        for(int j=1;j<(n+1);j++)
        {
            if(A[i-1] == B[j-1])
            {
                DP[i][j] = DP[i-1][j-1] + 1;
                
            }
            else
            {
                DP[i][j] = max(DP[i-1][j],DP[i][j-1]);
                
            }
        }


    }
    
    return DP[m][n];
}

//reconstructing LCS from score table
char* lcs_reconstruct(short **DP_Results, char *string_A, char *string_B)  {
  int row = strlen(string_A);
  int col = strlen(string_B);
  
  //number of elements in longest subsequence
  int len = DP_Results[row][col];
  char *LCS = (char *)malloc(len);
  
  //initialize values in array
  memset(LCS, ' ', len);
  
  //algorithm to reconstruct the lcs, uses the method from hw6
  while(row > 0 && col > 0) {
    if (string_A[row] == string_B[col]) {
      LCS[len] = string_A[row];
      --len;
      --row;
      --col;
    }
    else if (DP_Results[row - 1][col] >= DP_Results[row][col - 1]) {
      --row;
    }
    else {
      --col;
    }
  }
    if(string_A[row] == string_B[col]) {
     LCS[len] = string_A[row];
   }
  return LCS;
}

int main(int argc, char *argv)
{

    FILE *fp;
    int len_a,len_b;

    fp = fopen("Input2.txt", "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    printf("Sequence lengths : %d %d\n", len_a, len_b);

    string_A = (char *) calloc((len_a+1), sizeof(char *));
    string_B = (char *) calloc((len_b+1), sizeof(char *));
    unique_chars_C = (char *) calloc((c_len+1), sizeof(char *));

    fscanf(fp, "%s %s %s", string_A,string_B,unique_chars_C);


    //allocate memory for DP Results
    DP_Results = (short **) calloc((len_a+1), sizeof(short *));
    for(int k=0;k<len_a+1;k++)
    {
        DP_Results[k] = (short *)calloc((len_b+1), sizeof(short));
        

    }


    printf("\nLength of LCS is: %d\n",lcs(DP_Results,string_A,string_B,len_a,len_b));
    char * LCS = lcs_reconstruct(DP_Results, string_A,string_B);
    printf("Longest subsequence is %s\n", LCS);

    free(DP_Results);
    return 0;
}
