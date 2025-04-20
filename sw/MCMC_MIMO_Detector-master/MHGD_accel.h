#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "math.h"
#include "time.h"
#include "MyComplex.h"




#define time_size (x_copy + 1)

extern MyComplex QPSK_Constellation[4];
extern MyComplex _16QAM_Constellation[16];

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
float MHGD_detect_accel(MyComplex* x_hat, int Nt, int Nr, int mu, MyComplex* H, MyComplex* y, float sigma2, int mmse_init, int lr_approx, int iter);
void MHGD_free_accel(int Nt, int Nr, int mu);
void Hx_accel(int Nt, int Nr, int mu, float dqam, MyComplex* x, MyComplex* x_result);