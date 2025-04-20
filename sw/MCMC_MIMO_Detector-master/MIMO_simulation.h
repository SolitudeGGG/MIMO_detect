#pragma once

#include "util.h"
#include "math.h"
#include "time.h"
#include "MHGD_accel.h"
#include "MMSE.h"
#include "MyComplex.h"
#include <stdio.h>
#include <stdlib.h>



extern float BER;
extern float SER;
extern char *target_file;
extern float* Throughput;

typedef enum
{
	MMSE = 0,
	MHGD,
	NAG_MCMC,
	EP,
	ML

}detect_type_t;

typedef enum
{
	rayleigh = 0,
	corr = 1

}channel_type_t;

typedef struct
{
	int mu;
	int Nt;
	int Nr;
	int samples;
	int sampler;
	int EP_iter;
	detect_type_t detect_type;
	channel_type_t channel_type;
	int mmse_init;
	int lr_approx;
	int SNR_start;
	int SNR_end;
	int SNR_step;
	int target_err_bits;
	int max_iter;
}MIMO_sys_t;

void Parameter_init(int Nt, int Nr, int mu, int iter, int samples, detect_type_t detect_type,int max_iter);
float MIMO_detect(MIMO_sys_t MIMO_sys, MyComplex* H, MyComplex* y, MyComplex* x_hat, float SNR);
float MMSE_detect(int Nt, int Nr, MyComplex* y, MyComplex* H, MyComplex* x_hat, float sigma2);
void QAM_Demodulation(MyComplex* x_hat, int Nt, int mu, int* bits_demod);
void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod);
void MIMO_free(MIMO_sys_t MIMO_sys);
float detection(MIMO_sys_t MIMO_sys, float SNR);
