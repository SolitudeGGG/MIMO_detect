#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "MIMO_simulation.h"
#include "MyComplex.h"

MIMO_sys_t MIMO_sys = {
.mu = 4,/*调制阶数*/
.Nt = 8,/*发射天线数*/
.Nr = 8,/*接收天线数*/
.sampler = 1,/*在实现C语言并行多线程之前，MHGD只能运行单采样器*/
.samples = 32,/*每个采样器的采样数*/
.detect_type = MHGD,/*检测类型，本工程只实现了MMSE和MHGD*/
.channel_type = rayleigh,/*仿真的信道类型*/
.mmse_init = 0,/*是否使用MMSE检测的结果作为MCMC采样的初始值*/
.lr_approx = 0,
.SNR_start = 5,/*仿真的起始SNR*/
.SNR_end = 25,/*仿真的终止SNR*/
.SNR_step = 5,/*仿真的SNR步长*/
.target_err_bits = 5000,
.max_iter = 100,/*希望仿真的最大轮数*/
};

#define SNR_len ((MIMO_sys.SNR_end - MIMO_sys.SNR_start) / MIMO_sys.SNR_step + 1)/*SNR长度*/ 
#pragma warning(disable:4996)

int main()
{
#ifdef run_detection/*在MIMO_sys定义的不同SNR下进行MIMO检测*/
	/*Initialization and memory allocation*/
	float* SNR = (float*)malloc(SNR_len * sizeof(float));
	float* BER = (float*)malloc(SNR_len * sizeof(float)); 
	Throughput = (float*)malloc(SNR_len * sizeof(float)); 
	int i = 0;
	for (i = 0; i < SNR_len; i++)
		SNR[i] = MIMO_sys.SNR_start + i * MIMO_sys.SNR_step;
	/*init all the parameters*/
	Parameter_init(MIMO_sys.Nt, MIMO_sys.Nr, MIMO_sys.mu, MIMO_sys.samples, MIMO_sys.samples, MIMO_sys.detect_type, MIMO_sys.max_iter);
	/*MIMO_simulation*/
	for (i = 0; i < SNR_len; i++)
		BER[i] = detection(MIMO_sys, SNR[i]);
	/*print the result*/
	for (i = 0; i < SNR_len; i++)
		printf("%.6f ", BER[i]);
	printf("\n");
	for (i = 0; i < SNR_len; i++)
		printf("%.6f ", Throughput[i]);/*输出所有SNR下的吞吐量*/
	/*free the memory*/
	MIMO_free(MIMO_sys);
#endif

	return 0;
}
