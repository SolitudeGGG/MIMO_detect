#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "util_1.h"
#include "math.h"
#include "time.h"
#include "MyComplex_1.h"

///// vitis hls 相关库文件
///#include "ap_fixed.h"	//ap_fixed<18,6,AP_TRN_ZERO,AP_SAT>		<W,I,Q,O,N>
///#include "ap_int.h"	//ap_int<N> or ap_uint<N>, 1<=N<=1024
///#include "hls_math.h"	//data_t s = hls::sinf(angle);
///#include "hls_stream.h"



#define time_size (x_copy + 1)

// for循环范围限定参数定义区
static const int Nr_min = 1;
static const int Nr_max = 10;
static const int Nr_2_min = 1;
static const int Nr_2_max = 100;

static const int Nt_min = 1;
static const int Nt_max = 10;
static const int Nt_2_min = 1;
static const int Nt_2_max = 100;

static const int iter_min = 5000;
static const int iter_max = 15000;



extern MyComplex QPSK_Constellation[4];
extern MyComplex _16QAM_Constellation[16];
extern MyComplex _64QAM_Constellation[64];

extern float* MHGD_start;
extern float* MHGD_end;
extern float* MHGD_dur;

typedef enum
{
	gradpre = 0,
	covar_cal,
	lr_cal,
	mmse_and_r_init,
	zgrad,
	randomwalk,
	map_op,
	new_r,
	accept,
	x_copy
}time_seperate_t;

void MHGD_Init_accel(int Nt, int Nr, int mu, int iter, MyComplex* constellation);
void MHGD_detect_accel(Myreal* x_hat_real, Myimage* x_hat_imag,int Nt, int Nr, int mu, MyComplex* H, MyComplex* y, float sigma2, int mmse_init, int lr_approx, int iter, MyComplex* v_tb);
void MHGD_free_accel(int Nt, int Nr, int mu);
void Hx_accel(int Nt, int Nr, int mu, float dqam, MyComplex* x, MyComplex* x_result);

void QAM_Demodulation(MyComplex* x_hat, int Nt, int mu, int* bits_demod);
void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);

