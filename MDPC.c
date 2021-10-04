#include "formulalib.h"
#include "matrixlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // To generate random number

int checkTime(clock_t begin) {
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  if (time_spent > 300) {
    printf("Time limt\n");
    return -1;
  }
  return 0;
}

int pivotGauss(int **M, int **M_reverse, int size) {
  int i, j, k;

  for (i = 0; i < size; i++) {
    if (M[i][i] == 0) {
      for (j = i + 1; j < size; j++) {
        if (M[j][i] == 1) {
          for (k = 0; k < size; k++) {
            M[i][k] = (M[i][k] + M[j][k]) % 2;
            M_reverse[i][k] = (M_reverse[i][k] + M_reverse[j][k]) % 2;
          }
        }
      }
    }

    for (j = 0; j < size; j++) {
      if (i != j && M[j][i] == 1) {
        for (k = 0; k < size; k++) {
          M[j][k] = (M[j][k] + M[i][k]) % 2;
          M_reverse[j][k] = (M_reverse[j][k] + M_reverse[i][k]) % 2;
        }
      }
    }

    if (M[i][i] == 0) {
      return -1;
    }
  }
  return 0;
}

int hammingweight(int **M, int sizeRow, int sizeCol) {
  int i, j, sum = 0;
  for (i = 0; i < sizeRow; i++) {
    for (j = 0; j < sizeCol; j++) {
      sum += M[i][j];
    }
  }
  return sum;
}

void printArray(int array[], int size) {
  int i;
  printf("\n");
  for (i = 0; i < size; i++) {
    printf("%d ", array[i]);
  }
  printf("\n");
}

void BitFlipping(int **h0, int **h1, int **s, int T, int t, int size) {
  int i, j, **u, **v, **h0_rot, **h1_rot, **h0_transpose, **h1_transpose, **H,
      **flipped_positions, **uv, **syndrome, **sum,
      **flipped_positions_transpose, **H_flipped_transpose, **uv_transpose,
      **H_uv_transpose;
  int size_H = 2 * size;
  u = initMatrix(u, 1, size);
  v = initMatrix(v, 1, size);
  for (i = 0; i < size; i++) {
    u[0][i] = 0;
    v[0][i] = 0;
  }

  h0_rot = initMatrix(h0_rot, size, size);
  h1_rot = initMatrix(h1_rot, size, size);

  rot(h0_rot, h0, size);
  rot(h1_rot, h1, size);

  h0_transpose = initMatrix(h0_transpose, size, size);
  h1_transpose = initMatrix(h1_transpose, size, size);

  transposeMatrix(h0_transpose, h0_rot, size, size);
  transposeMatrix(h1_transpose, h0_rot, size, size);

  printf("h0_rot:\n");
  printMatrix(h0_rot, size, size);
  printf("h0_transpose:\n");
  printMatrix(h0_transpose, size, size);

  printf("h1_rot:\n");
  printMatrix(h1_rot, size, size);
  printf("h1_transpose:\n");
  printMatrix(h1_transpose, size, size);

  H = initMatrix(H, size, size_H);
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      H[i][j] = h0_transpose[i][j];
      H[i][j + size] = h1_transpose[i][j];
    }
  }

  printf("H:\n");
  printMatrix(H, size, size_H);

  syndrome = initMatrix(syndrome, 1, size);
  cloneMatrix(syndrome, s, 1, size);

  sum = initMatrix(sum, 1, size_H);
  flipped_positions = initMatrix(flipped_positions, 1, size_H);
  flipped_positions_transpose =
      initMatrix(flipped_positions_transpose, size_H, 1);
  H_flipped_transpose = initMatrix(H_flipped_transpose, size, 1);

  int count = 100;
  while (
      ((hammingweight(u, 1, size) != t) || (hammingweight(v, 1, size) != t)) &&
      hammingweight(syndrome, 1, size) != 0) {
    multiMatrixNoModulo(sum, syndrome, H, 1, size_H, size);
    for (i = 0; i < size; i++) {
      flipped_positions[0][i] = 0;
      if (sum[0][i] >= T) {
        flipped_positions[0][i] = (flipped_positions[0][i] + 1) % 2;
      }
      u[0][i] = (flipped_positions[0][i] + u[0][i]) % 2;

      int j = i + size;
      flipped_positions[0][j] = 0;
      if (sum[0][j] >= T) {
        flipped_positions[0][j] = (flipped_positions[0][j] + 1) % 2;
      }
      v[0][i] = (flipped_positions[0][j] + v[0][i]) % 2;
    }

    transposeMatrix(flipped_positions_transpose, flipped_positions, size_H, 1);
    multiMatrix(H_flipped_transpose, H, flipped_positions_transpose,
                        size, 1, size_H);

    // subtractMatrixNoModulo(syndrome, syndrome, H_flipped_transpose, size, 1);
    addMatrixTrans(syndrome, syndrome, H_flipped_transpose, 1, size);

    count--;
    if (count == 0) {
      printf("Fail 100 try times\n");
      break;
    }
  }

  printf("Check\n");
  uv = initMatrix(uv, 1, size_H);
  uv_transpose = initMatrix(uv_transpose, size_H, 1);

  for (i = 0; i < size; i++) {
    uv[0][i] = u[0][i];
    uv[0][i + size] = v[0][i];
  }
  printf("Check 1\n");
  transposeMatrix(uv_transpose, uv, size_H, 1);

  printf("Check 2\n");
  H_uv_transpose = initMatrix(H_uv_transpose, size, 1);
  multiMatrixNoModulo(H_uv_transpose, H, uv_transpose, size, 1, size_H);
  addMatrixTrans(s, s, H_uv_transpose, 1, size);

  printf("Check 3\n");
  if (hammingweight(s, 1, size) != 0) {
    printf("Fail to find e0, e1\n");
  } else {
    printf("u:\n");
    printMatrix(u, 1, size);

    printf("v:\n");
    printMatrix(v, 1, size);
  }

  // subtractMatrixNoModulo(s, s, H_uv_transpose, 1, size);

  free(u);
  free(v);
  free(h0_rot);
  free(h1_rot);
  free(h0_transpose);
  free(h1_transpose);
  free(H);
  free(syndrome);
  free(sum);
  free(flipped_positions);
  free(flipped_positions_transpose);
  free(H_flipped_transpose);
  free(uv);
  free(uv_transpose);
  free(H_uv_transpose);
}

void retrive_e_From_c(int **e0, int **e1, int **c, int **h0, int **h1, int T,
                      int t, int size) {
  int **s, **rot_c;
  s = initMatrix(s, 1, size);
  rot_c = initMatrix(rot_c, size, size);
  rot(rot_c, c, size);
  multiMatrix(s, h0, rot_c, 1, size, size);
  printf("s:\n");
  printMatrix(s, 1, size);

  BitFlipping(h0, h1, s, T, t, size);

  free(s);
  free(rot_c);
}

int main() {
  int i, j;
  int **h, **h0, **h1, **e0, **e1, **h0_reverse, **h0_reverse_rot, **e1_rot,
      **c, **he1;
  int n; // Length of h0 and h1
  int w; // Weight of h0 and h1
  int k;
  int t;
  int T;
  int flag_found = 0;

  // Param
  n = 4813; // Length of h0 and h1
  t = 39;
  w = 39;
  T = 26;

  h0 = initMatrix(h0, 1, n);
  h0_reverse = initMatrix(h0_reverse, 1, n);

  srand(time(NULL));
  int inversible = -1;
  while (inversible == -1) {
    randomMatrixWithHammingWeight(h0, 1, n, w);

    inversible = invertFormula(h0, h0_reverse, n, w);
    if (inversible == -1) {
      printf("Irreversible\n");
    } else {
      printf("h0:\n");
      printMatrix(h0, 1, n);
      printf("h_reverse:\n");
      printMatrix(h0_reverse, 1, n);
    }
  }

  h1 = initMatrix(h1, 1, n);
  randomMatrixWithHammingWeight(h1, 1, n, w);

  printf("h1:\n");
  printMatrix(h1, 1, n);

  h = initMatrix(h, 1, n);
  h0_reverse_rot = initMatrix(h0_reverse_rot, n, n);
  rot(h0_reverse_rot, h0_reverse, n);
  multiMatrix(h, h1, h0_reverse_rot, 1, n, n);
  printf("H:\n");
  printMatrix(h, 1, n);

  e0 = initMatrix(e0, 1, n);
  e1 = initMatrix(e1, 1, n);

  randomMatrixWithHammingWeight(e0, 1, n, t);
  randomMatrixWithHammingWeight(e1, 1, n, t);

  c = initMatrix(c, 1, n);
  he1 = initMatrix(he1, 1, n);
  e1_rot = initMatrix(e1_rot, n, n);
  rot(e1_rot, e1, n);
  multiMatrix(he1, h, e1_rot, 1, n, n);
  addMatrix(c, e0, he1, 1, n);

  printf("e0:\n");
  printMatrix(e0, 1, n);

  printf("e1:\n");
  printMatrix(e1, 1, n);

  printf("c:\n");
  printMatrix(c, 1, n);

  // receive c
  retrive_e_From_c(e0, e1, c, h0, h1, T, t, n);

  // clock_t begin = clock();
  // int fail = 0;
  /* here, do your time-consuming job */

  free(h);
  free(h0);
  free(h1);
  free(e0);
  free(e1);
  free(h0_reverse);
  free(h0_reverse_rot);
  free(e1_rot);
  free(c);
  free(he1);

  return 0;
}