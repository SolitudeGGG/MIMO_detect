#pragma once
#include "MyComplex.h"
#include "hls_math.h"

typedef struct _IO_FILE FILE;
#define RAND_MAX 0x7fffffff


#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/*算法参数*/
static const int Ntr_1 = 8;/*发射&接收天线数*/
static const int Ntr_2 = 64;/*NtorNr^2*/
static const int iter_1 = 32;/*采样器的采样数*/
static const int mu_1 = 4;/*调制阶数*/
static const int mmse_init_1 = 0;/*是否使用MMSE检测的结果作为MCMC采样的初始值*/
static const int lr_approx_1 = 0;
static const int max_iter_1 = 100;/*希望仿真的最大轮数*/

void read_gaussian_data(const char* filename, MyComplex* array, int n, int offset);
void QAM_Demodulation(MyComplex* x_hat, int Nt, int mu, int* bits_demod);
void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
int argmin(float* array, int n);
int unequal_times(int* array1, int* array2, int n);
MyComplex complex_conjugate(MyComplex a);
void c_eye_generate(MyComplex* Mat, int n, float val);
void c_matmultiple(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res);
void my_complex_add(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]);
void initMatrix(MyComplex* A, int row, int col);
MyComplex complex_divide(MyComplex a, MyComplex b);
MyComplex complex_multiply(MyComplex a, MyComplex b);
MyComplex complex_add(MyComplex a, MyComplex b);
MyComplex complex_subtract(MyComplex a, MyComplex b);
void my_complex_sub(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]);
void MulMatrix(const MyComplex* A, const MyComplex* B, MyComplex* C, int row, int AB_rc, int col);
void Inverse_LU(MyComplex* A, int row, int col);
void my_complex_scal(const int N, const float alpha, MyComplex* X, const int incX);
void map(int mu, int Nt, float dqam, MyComplex* x, MyComplex* x_hat);
void generateUniformRandoms_int(int Nt, int* x_init, int mu);
void my_complex_copy(const int N, const MyComplex* X, const int incX, MyComplex* Y, const int incY);
void generateUniformRandoms_float(int Nt, float* p_uni);
void out_hw(const MyComplex* X, const int incX, float* Y_real, float* Y_imag, int incY);
void MHGD_detect_accel_sw(
    float* x_hat_real, float* x_hat_imag, 
    float* H_real, float* H_imag, 
    float* y_real, float* y_imag, 
    float sigma2,
    float* v_tb_real, float* v_tb_imag
);
