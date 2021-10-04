void modFormula(int dividend[], int divisor[], int quotient[], int size) {
  int i, j = size - 1;

  while (divisor[j] == 0) {
    j--;
  }

  for (i = size - 1; i > -1; i--) {
    if (dividend[i] != 0) {
      // printf("%d, %d\n", i , j);
      if (i >= j) {
        int rs = i - j;
        quotient[rs] = 1;

        int k;
        for (k = 0; k < size; k++) {
          if (divisor[k] != 0) {
            int temp = k + rs;
            dividend[temp] = (dividend[temp] + divisor[k]) % 2;
          }
        }
      } else {
        break;
      }
    }
  }
}

void swapFormula(int f1[], int f2[], int size) {
  int i, j;
  for (i = 0; i < size; i++) {
    int temp = f1[i];
    f1[i] = f2[i];
    f2[i] = temp;
  }
}

int endMod(int f[], int size) {
  int i;
  for (i = 1; i < size; i++) {
    if (f[i] != 0) {
      return 0;
    }
  }

  if (f[0] == 0) {
    return -1;
  }
  return 1;
}

int invertFormula(int **M, int **M_reverse, int size, int weight) {
  int dividend[size], i, j;
  int divisor[size], v0[size], v1[size], quotient[size], v_temp[size];
  for (i = 0; i < size; i++) {
    dividend[i] = 0;
    divisor[i] = M[0][i];
    v0[i] = 0;
    v1[i] = 0;
    v_temp[i] = 0;
    quotient[i] = 0;
  }
  dividend[0] = 1;
  dividend[weight] = 1;
  v1[0] = 1;

  modFormula(divisor, dividend, quotient, size);

  while (1) {
    for (i = 0; i < size; i++) {
      quotient[i] = 0;
    }
    // printArray(divisor, size);
    // printArray(dividend, size);
    modFormula(dividend, divisor, quotient, size);
    // printArray(divisor, size);
    // printArray(dividend, size);

    int end = endMod(dividend, size);
    if (end == -1) {
      return -1;
    }

    for (i = 0; i < size; i++) {
      if (v1[i] != 0) {
        for (j = 0; j < size; j++) {
          if (quotient[j] != 0) {
            int temp = i + j;
            v_temp[temp] =
                (v_temp[temp] + (v1[i] * quotient[j]) + v0[temp]) % 2;
          }
        }
      }
    }
    swapFormula(v0, v1, size);
    swapFormula(v1, v_temp, size);
    // printArray(v0, size);
    // printArray(v1, size);

    if (end == 1) {
      break;
    }

    swapFormula(dividend, divisor, size);
  }

  for (i = 0; i < size; i++) {
    M_reverse[0][i] = v1[i];
  }

  return 0;
}
