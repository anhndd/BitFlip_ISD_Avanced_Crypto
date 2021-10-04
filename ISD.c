#include <stdio.h>
#include <stdlib.h>
#include <time.h> // To generate random number
#include <unistd.h>
#include <string.h>

int checkTime(clock_t begin){
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(time_spent > 300){
        printf("Time limt\n");
        return -1;
    }
    return 0;
}

int valueinarray(int val, int arr[], int size)
{
    int i;
    for(i = 0; i < size; i++)
    {
        if(arr[i] == val)
            return 1;
    }
    return 0;
}

void swapRowMatrix(int** M, int size, int row1, int row2){
    int i;
    for(i=0;i<size;i++){
        int temp = M[row1][i];
        M[row1][i] = M[row2][i];
        M[row2][i] = temp;
    }
}

void cloneMatrix(int**rs, int **M, int sizeRow, int sizeCol){
    int i,j;
    for(i=0; i<sizeRow; i++){
        for(j=0; j<sizeCol; j++){
            rs[i][j] = M[i][j];
        }
    }
}

void randomMatrix(int**Mat, int sizeRow, int sizeCol){
    int i,j;
    for(i=0; i<sizeRow; i++){
        for(j=0; j<sizeCol; j++){
            Mat[i][j] = rand() % 2;
        }
    }
}

void generateInitMatrix(int**Mat, int size){
    int i,j;
    for(i=0; i<size; i++){
        for(j=0; j<size; j++){
            if(i==j){
                Mat[i][j] = 1;
            } else{
                Mat[i][j] = 0;
            }
        }
    }
}

void printMatrix(int** M, int sizeRow, int sizeCol){
    int i,j;
    for(i=0; i<sizeRow; i++){
        for(j=0; j<sizeCol; j++){
            printf("%d ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int pivotGauss(int** M,int** M_inverse, int size){
    int i,j,k;

    for(i=0;i< size;i++){
        if(M[i][i] == 0){
            for(j=i+1;j< size;j++){
                if(M[j][i] == 1){
                    for(k=0;k< size;k++){
                        M[i][k] = (M[i][k] + M[j][k]) % 2;
                        M_inverse[i][k] = (M_inverse[i][k] + M_inverse[j][k]) % 2;
                    }    
                }
            }
        }

        for(j=0;j< size;j++){
            if (i!=j && M[j][i] == 1){
                for(k=0;k< size;k++){
                    M[j][k] = (M[j][k] + M[i][k]) % 2;
                    M_inverse[j][k] = (M_inverse[j][k] + M_inverse[i][k]) % 2;
                }
                
            }
        }

        if (M[i][i] == 0) {
            return -1;
        }
    }
    return 0;
}

int** initMatrix(int **M, int sizeRow, int sizeCol){
    int i;
    M = (int**) malloc(sizeRow* sizeof(int*));
    for(i=0; i<sizeRow; i++) {
        M[i] = (int *) malloc(sizeCol * sizeof(int));
    }
    return M;
}

int** multiMatrix(int** rs,int **M1, int** M2, int sizeRow, int sizeCol, int sizeMiddle){
    int i, j;

    for(i = 0;i < sizeRow; i++){
        for(j = 0;j < sizeCol; j++){
            rs[i][j] = 0;
        }
    }
    for(i = 0;i < sizeRow; i++){
        for(j = 0;j < sizeCol; j++){
            int k;
            for(k = 0;k < sizeMiddle; k++){
                rs[i][j] += M1[i][k] * M2[k][j];
                rs[i][j] %= 2;
            }
        }
    }
}

int **pickRandomColMatrix(int**rs, int**M, int sizeRow, int sizeCol){
    
    int i,k;
    int j = 0;
    int chosenCol[sizeRow];

    for(i=0;i<sizeRow;i++){
        int randomCol;
        do {
            randomCol = rand() % sizeCol;
        } while(valueinarray(randomCol, chosenCol, sizeRow) == 1);
        chosenCol[j] = randomCol;
        for(k=0;k<sizeRow;k++) {
            rs[k][j] = M[k][randomCol];
        }
        j++;
    }
    return rs;
}

int hammingweight(int **M, int sizeRow, int sizeCol) {
    int i, j, sum = 0;
    for(i=0;i<sizeRow; i++){
        for(j=0;j<sizeCol; j++){
            sum += M[i][j];
        }
    }
    return sum;
}

void main()
{
    int i, j;
    int **h, **x, **s, **a, **a_inverse, **xb, **a_clone;
    int n; // Length of h0 and h1
    int w; // Weight of h0 and h1
    int k;
    int nk;
    int flag_found = 0;

//     // Param
    n = 1000; // Length of h0 and h1
    w = 10; // Weight of h0 and h1
    k = 500;
    nk = n - k;

    srand(time(NULL));
    h = initMatrix(h, nk, n);
    x = initMatrix(x, n, 1);
    s = initMatrix(s, nk, 1);
    a = initMatrix(a, nk, nk);
    a_clone = initMatrix(a_clone, nk, nk);
    a_inverse = initMatrix(a_inverse, nk, nk);
    xb = initMatrix(xb, n, 1);

    randomMatrix(h, nk, n);
    randomMatrix(x, n, 1);
    multiMatrix(s, h, x, nk, 1, n);

    printf("H:\n");
    printMatrix(h, nk, n);

    printf("X:\n");
    printMatrix(x, n, 1);

    printf("S:\n");
    printMatrix(s, nk, 1);
    

    clock_t begin = clock();
    int fail = 0;
    /* here, do your time-consuming job */



    do {
        do {
            generateInitMatrix(a_inverse, nk);
            pickRandomColMatrix(a, h, nk, n);
            // cloneMatrix(a_clone,a,nk,nk);
            // printf("checking invertible\n");
            if(checkTime(begin) != 0){
                fail = 1;
                break;
            }
        } while(pivotGauss(a, a_inverse, nk) != 0);
        if(fail != 0){
            break;
        }

        printf("A:\n");
        printMatrix(a, nk, nk);

        printf("A inverse:\n");
        printMatrix(a_inverse, nk, nk);

        multiMatrix(xb, a_inverse, s, nk, 1, nk);

        printf("X':\n");
        printMatrix(xb, nk, 1);

        int hammingX = hammingweight(xb, nk, 1);

        printf("Check %d\n", hammingX);
        if (hammingX <= w){
            flag_found = 1;
            clock_t end = clock();
            double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            printf("Found - time: %lf\n", time_spent);
        } else{
            printf("Not satisfiable\n");
        }

        
    } while(flag_found == 0);

    free(h);
    free(s);
    free(x);
    free(a);
    free(a_inverse);
    free(xb);
    free(a_clone);
}
