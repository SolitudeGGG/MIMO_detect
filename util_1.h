#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sys_config_1.h"
#include "MyComplex_1.h"
#include "time.h"

#pragma warning(disable:4996)

extern FILE* f;
void c_eye_generate(MyComplex* Mat, int n, float val);
void c_ones_generate(MyComplex* Mat, int row, int col);
void Inverse_LU(MyComplex* A, int row, int col);

void my_complex_add(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]);
void c_matmultiple(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res);

void generateUniformRandoms_int(int Nt, int* x_init, int mu);
void generateUniformRandoms_float(int Nt, like_float p_uni[9]);
void my_complex_scal(const int N, const like_float alpha, MyComplex* X, const int incX);
void read_gaussian_data(const char* filename, MyComplex* array, int n, int offset);
int argmin(float* array, int n);
int unequal_times(int* array1, int* array2, int n);
void map(int mu, int Nt, like_float& dqam, MyComplex* x, MyComplex* x_hat);

void my_complex_sub(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]);
void my_complex_copy(const int N, const MyComplex* X, const int incX, MyComplex* Y, const int incY);