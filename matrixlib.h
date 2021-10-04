#include <stdio.h>
#include <stdlib.h>

int valueinarray(int val, int arr[], int size) {
  int i;
  for (i = 0; i < size; i++) {
    if (arr[i] == val)
      return 1;
  }
  return 0;
}

int **initMatrix(int **M, int sizeRow, int sizeCol) {
  int i;
  M = (int **)malloc(sizeRow * sizeof(int *));
  for (i = 0; i < sizeRow; i++) {
    M[i] = (int *)malloc(sizeCol * sizeof(int));
  }
  return M;
}

void swapRowMatrix(int **M, int size, int row1, int row2) {
  int i;
  for (i = 0; i < size; i++) {
    int temp = M[row1][i];
    M[row1][i] = M[row2][i];
    M[row2][i] = temp;
  }
}

void addMatrix(int **rs, int **M1, int **M2, int sizeRow, int sizeCol) {
  int i, j;
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = M1[i][j] + M2[i][j];
      rs[i][j] %= 2;
    }
  }
}

void addMatrixTrans(int **rs, int **M1, int **M2, int sizeRow, int sizeCol) {
  int i, j;

  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = M1[i][j] + M2[j][i];
      rs[i][j] %= 2;
    }
  }
}

void subtractMatrixNoModulo(int **rs, int **M1, int **M2, int sizeRow,
                            int sizeCol) {
  int i, j;
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = M1[i][j] - M2[i][j];
    }
  }
}

void multiMatrix(int **rs, int **M1, int **M2, int sizeRow, int sizeCol,
                 int sizeMiddle) {
  int i, j;

  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = 0;
    }
  }
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      int k;
      for (k = 0; k < sizeMiddle; k++) {
        rs[i][j] += M1[i][k] * M2[k][j];
        rs[i][j] %= 2;
      }
    }
  }
}

void multiMatrixNoModulo(int **rs, int **M1, int **M2, int sizeRow, int sizeCol,
                         int sizeMiddle) {
  int i, j;

  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = 0;
    }
  }
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      int k;
      for (k = 0; k < sizeMiddle; k++) {
        rs[i][j] += M1[i][k] * M2[k][j];
      }
    }
  }
}

void cloneMatrix(int **rs, int **M, int sizeRow, int sizeCol) {
  int i, j;
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      rs[i][j] = M[i][j];
    }
  }
}

void generateInitMatrix(int **Mat, int size) {
  int i, j;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (i == j) {
        Mat[i][j] = 1;
      } else {
        Mat[i][j] = 0;
      }
    }
  }
}

void randomMatrix(int **Mat, int sizeRow, int sizeCol) {
  int i, j;
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      Mat[i][j] = rand() % 2;
    }
  }
}

void randomMatrixWithHammingWeight(int **Mat, int sizeRow, int sizeCol,
                                   int weight) {
  int i, j;

  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      Mat[i][j] = 0;
    }
  }

  while (weight > 0) {
    int row = rand() % sizeRow;
    int col = rand() % sizeCol;
    if (Mat[row][col] != 1) {
      Mat[row][col] = 1;
      weight--;
    }
  }
}

int **pickRandomColMatrix(int **rs, int **M, int sizeRow, int sizeCol) {

  int i, k;
  int j = 0;
  int chosenCol[sizeRow];
  // for(i=0;i<sizeRow;i++){
  //     chosenCol[i] = -1;
  // }

  for (i = 0; i < sizeRow; i++) {
    int randomCol;
    do {
      randomCol = rand() % sizeCol;
    } while (valueinarray(randomCol, chosenCol, sizeRow) == 1);
    chosenCol[j] = randomCol;
    // printf("\n%d %d\n",j, randomCol);
    for (k = 0; k < sizeRow; k++) {
      rs[k][j] = M[k][randomCol];
    }
    j++;
  }
  return rs;
}

void rot(int **mat_rot, int **vec, int size) {
  int i, j;

  // First row
  for (j = 0; j < size; j++)
    mat_rot[0][j] = vec[0][j];

  // Rotate the next rows
  for (i = 1; i < size; i++) {
    mat_rot[i][0] = mat_rot[i - 1][size - 1];
    for (j = 1; j < size; j++)
      mat_rot[i][j] = mat_rot[i - 1][j - 1];
  }
}

void transposeMatrix(int **rs, int **M, int sizeRowRs, int sizeColRs) {
  int i, j;
  for (i = 0; i < sizeRowRs; i++) {
    for (j = 0; j < sizeColRs; j++) {
      rs[i][j] = M[j][i];
    }
  }
}

void printMatrix(int **M, int sizeRow, int sizeCol) {
  // int i, j;
  // for (i = 0; i < sizeRow; i++) {
  //   for (j = 0; j < sizeCol; j++) {
  //     printf("%d ", M[i][j]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");
}
