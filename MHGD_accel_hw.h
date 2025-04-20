#pragma once
#include "MyComplex_1.h"
#include "hls_math.h"

typedef struct _IO_FILE FILE;

#define LCG_A 1664525
#define LCG_C 1013904223
#define LIMIT_MAX 0x7fffffff
#define LIMIT_MAX_F 2147483647.0f
#define M_PI 3.14159265358979323846
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

void read_gaussian_data_hw(const char* filename, MyComplex* array, int n, int offset);
void QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int mu, int* bits_demod);
void QPSK_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod);
void _16QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod);void _64QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod);
int argmin_hw(like_float* array, int n);
int unequal_times_hw(int* array1, int* array2, int n);

unsigned int lcg_rand_hw();
like_float lcg_rand_1_hw();
void get_dqam_hw(like_float dqam);
void c_eye_generate_hw(MyComplex* Mat, float val);
void data_local(MyComplex* H, MyComplex* y, MyComplex* v_tb, Myreal* H_real, Myimage* H_imag, Myreal* y_real, Myimage* y_imag, Myreal* v_tb_real, Myimage* v_tb_imag);
void c_matmultiple_hw(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res);
void my_complex_add_hw(const MyComplex a[], const MyComplex b[], MyComplex r[]);
void initMatrix_hw(MyComplex* A);
MyComplex complex_divide_hw(MyComplex a, MyComplex b);
MyComplex complex_multiply_hw(MyComplex a, MyComplex b);
MyComplex complex_add_hw(MyComplex a, MyComplex b);
MyComplex complex_subtract_hw(MyComplex a, MyComplex b);
void my_complex_sub_hw(const MyComplex a[], const MyComplex b[], MyComplex r[]);
void MulMatrix_hw(const MyComplex* A, const MyComplex* B, MyComplex* C);
void Inverse_LU_hw(MyComplex* A);
void my_complex_scal_hw(const like_float alpha, MyComplex* X, const int incX);
void map_hw(like_float dqam, MyComplex* x, MyComplex* x_hat);
void generateUniformRandoms_int_hw(int* x_init);
void my_complex_copy_hw(const MyComplex* X, const int incX, MyComplex* Y, const int incY);
void generateUniformRandoms_float_hw(like_float* p_uni);
void grad_preconditioner_updater_hw(MyComplex* H, int transA, int transB, MyComplex* HH_H, MyComplex* sigma2eye, MyComplex* grad_preconditioner, int Nt, int Nr);
void learning_rate_line_search_hw(int lr_approx, MyComplex* H, int transA, int transB, MyComplex* grad_preconditioner, int Nr, int Nt, MyComplex* temp_NtNr, MyComplex* pmat);
void x_initialize_hw(int mmse_init, MyComplex* sigma2eye, int Nt, int Nr, float sigma2, MyComplex* HH_H, MyComplex* temp_NtNt, MyComplex* H, int transA, int transB, MyComplex* y, MyComplex* temp_Nt, MyComplex* x_mmse, like_float dqam, MyComplex* x_hat, int* x_init, MyComplex* constellation_norm);
void r_hw(MyComplex* H, int transB, int transA, MyComplex* x_hat, int Nr, int Nt, MyComplex* r, MyComplex* y);
void r_cal_hw(MyComplex* r, int transA, int transB, int Nr, int Nt, MyComplex* temp_1, MyComplex* x_hat, MyComplex* x_survivor, like_float r_norm, like_float r_norm_survivor);
void lr_hw(int lr_approx, MyComplex* pmat, int transB, int transA, MyComplex* r, int Nr, int Nt, MyComplex* pr_prev, MyComplex* temp_1, MyComplex* _temp_1, like_float lr);
void z_grad_hw(MyComplex* H, int transA, int transB, MyComplex* temp_Nt, MyComplex* grad_preconditioner, MyComplex* z_grad, like_float lr, MyComplex* x_hat);
void gauss_add_hw(MyComplex* v, MyComplex* v_tb, int offset, like_float step_size, MyComplex* z_grad, MyComplex* z_prop);
void r_newnorm_hw(MyComplex* H, int transB, MyComplex* x_prop, int transA, MyComplex* temp_Nr, MyComplex* y, MyComplex* r_prop, MyComplex* temp_1, like_float r_norm_prop);
void survivor_hw(like_float r_norm_survivor, like_float r_norm_prop, MyComplex* x_prop, MyComplex* x_survivor);
void acceptance_hw(int transB, int transA, like_float r_norm_prop, like_float r_norm, like_float log_pacc, like_float p_acc, like_float* p_uni, MyComplex* x_prop , MyComplex* x_hat, MyComplex* r_prop, MyComplex* r, MyComplex* pmat, MyComplex* pr_prev, MyComplex* temp_1, MyComplex* _temp_1, like_float lr, like_float step_size, like_float dqam, like_float alpha);
void out_hw(const MyComplex* X, const int incX, Myreal* Y_real, Myimage* Y_imag, int incY);

void Inverse_LU(MyComplex_f* A);
void initMatrix(MyComplex_f* A);
void MulMatrix(const MyComplex_f* A, const MyComplex_f* B, MyComplex_f* C);
MyComplex_f complex_divide(MyComplex_f a, MyComplex_f b);
MyComplex_f complex_multiply(MyComplex_f a, MyComplex_f b);
MyComplex_f complex_add(MyComplex_f a, MyComplex_f b);
MyComplex_f complex_subtract(MyComplex_f a, MyComplex_f b);
MyComplex_f complex_divide(MyComplex_f a, MyComplex_f b);

void MHGD_detect_accel_hw(
    Myreal* x_hat_real, Myimage* x_hat_imag, 
    Myreal* H_real, Myimage* H_imag, 
    Myreal* y_real, Myimage* y_imag, 
    float sigma2,
    Myreal* v_tb_real, Myimage* v_tb_imag
);
