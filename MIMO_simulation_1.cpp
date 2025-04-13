#include "MIMO_simulation_1.h"
#include "MyComplex_1.h"
#include <string.h>
#include "MMSE_1.h"
#include "MHGD_accel_1.h"
#include <stdio.h>

typedef unsigned char uint8_t; // 定义uint8_t

float detect_time = 0;
float start = 0;
float end = 0;
MyComplex QPSK_Constellation[4];
MyComplex _16QAM_Constellation[16];
MyComplex _64QAM_Constellation[64];
char *target_file;
float* Throughput;

int* bits;
MyComplex* x;
MyComplex* H;
MyComplex* y;
MyComplex* noise;
// MyComplex* x_hat;
Myreal* x_hat_real;
Myimage* x_hat_imag;
int* bits_demod;
MyComplex* H_real;
MyComplex* H_image;
MyComplex* noise_real;
MyComplex* noise_image;
MyComplex* input_H;//信道H
MyComplex* input_y;//接收矢量y
int* origin_bits;

/*参数初始化与指针内存空间分配*/
void Parameter_init(int Nt, int Nr, int mu, int iter, int samples, detect_type_t detect_type, int max_iter)
{
	int i = 0;	float real = 0; float imag = 0;
	
	/*根据调制阶数初始化星座点*/
	switch (mu)
	{
	case 2:
		QPSK_Constellation[0].real = -1; QPSK_Constellation[0].imag = -1;
		QPSK_Constellation[1].real = -1; QPSK_Constellation[1].imag = 1;
		QPSK_Constellation[2].real = 1; QPSK_Constellation[2].imag = -1;
		QPSK_Constellation[3].real = 1; QPSK_Constellation[3].imag = 1;
		break;
	case 4:
		f = fopen("/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt", "r");	
		while (!feof(f))
		{
			fscanf(f, "%f %f\n", &real, &imag);
			_16QAM_Constellation[i].real = real; _16QAM_Constellation[i].imag = imag;
			i++;
		}
		fclose(f);
		break;
	case 6:
		f = fopen("/home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt", "r");
		while (!feof(f))
		{
			fscanf(f, "%f %f\n", &real, &imag);
			_64QAM_Constellation[i].real = real; _64QAM_Constellation[i].imag = imag;
			i++;
		}
		fclose(f);
		break;
	default:
		break;
	}

	/*针对不同的检测方式(以及不同调制阶数）提前分配相应的指针空间*/
	switch (detect_type)
	{
	case MMSE:
		//MMSE_Init(Nt, Nr);
		//printf("not mmse");
		break;
	case MHGD:
		switch (mu)
		{
		case 2:
			MHGD_Init_accel(Nt, Nr, mu, iter, QPSK_Constellation);
			break;
		case 4:
			MHGD_Init_accel(Nt, Nr, mu, iter, _16QAM_Constellation);
			break;
		case 6:
			MHGD_Init_accel(Nt, Nr, mu, iter, _64QAM_Constellation);
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
	/*提前分配指针空间*/
	bits = (int*)malloc(sizeof(int) * Nt * mu);/*transmitted bits*/
	x = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	H = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	y = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	noise = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	// x_hat = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	x_hat_real = (Myreal*)malloc(Nt * sizeof(Myreal));
	x_hat_imag = (Myreal*)malloc(Nt * sizeof(Myimage));
	bits_demod = (int*)malloc(sizeof(int) * Nt * mu);//received bits after demodulation
	H_real = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	H_image = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	noise_real = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	noise_image = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	/*集中检测用到的矩阵*/
	input_H = (MyComplex*)malloc(max_iter * Nt * Nr * sizeof(MyComplex));
	input_y = (MyComplex*)malloc(max_iter * Nr * sizeof(MyComplex));
	origin_bits = (int*)malloc(max_iter * Nt * mu * sizeof(int));
	
	/*记录当前时间和相关的实验参数*/
	f = fopen("/home/ggg_wufuqi/hls/MHGD/C_time_seperate.csv", "a");
	time_t timep;
	struct tm* p;
	time(&timep);
	p = gmtime(&timep);
	fprintf(f, "%d-\t%d-\t%d-\t%d:\t%d:\t%d,", 1900 + p->tm_year, 1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec);
	fprintf(f, "Nt:%d,Nr:%d,mu:%d,samples:%d\n", Nt, Nr, mu, samples);
	fclose(f);
}

/*分配的指针空间在仿真结束后需要统一释放*/
void MIMO_free(MIMO_sys_t MIMO_sys)
{
	free(bits);
	free(x);
	free(H);
	free(y);
	free(noise);
	// free(x_hat);
	free(x_hat_real);
	free(x_hat_imag);
	free(H_real);
	free(H_image);
	free(noise_real);
	free(noise_image);
	free(bits_demod);

	switch (MIMO_sys.detect_type)
	{
	case MMSE:
	    //printf("not mmse");
		break;
	case MHGD:
	MHGD_free_accel(MIMO_sys.Nt,MIMO_sys.Nr,MIMO_sys.mu);
	break;
	default:
		break;
	}
}


/*
 * @brief  根据接收信号y 信道H检测原始发送信号x，解调成比特后打印输出到文件中，并返回本次SNR下检测的误比特率BER
 * @note
 * @param  SNR 信噪比
 * @retval 本次检测的误比特率BER
 */
float detection(MIMO_sys_t MIMO_sys, float SNR)
{
	int Nt = MIMO_sys.Nt; int Nr = MIMO_sys.Nr; int mu = MIMO_sys.mu;
	int total_error_bits = 0; int total_bits = 0; int count = 0;
	float BER = 0;
	int max_iter = MIMO_sys.max_iter;
	detect_time = 0;
	/*字符串拼接，根据信噪比不同写入不同的文本文件*/
	char bits_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/reference_file/bits_SNR=";
	char H_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=";
	char y_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=";
	char bits_output_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=";
	char txt[] = ".txt";
	MyComplex x_hat[8];
	sprintf(bits_file + strlen(bits_file), "%d", (int)SNR);
	strcat(bits_file, txt);
	sprintf(H_file + strlen(H_file), "%d", (int)SNR);
	strcat(H_file, txt);
	sprintf(y_file + strlen(y_file), "%d", (int)SNR);
	strcat(y_file, txt);
	sprintf(bits_output_file + strlen(bits_output_file), "%d", (int)SNR);
	strcat(bits_output_file, txt);

	float real_temp;
	float imag_temp;
	int b;
	int i = 0; int j = 0;
	/*读取输入测试文件数据（信道数据、接受信号数据、比特数据）*/
	printf("SNR=%f loading file...", SNR);
	f = fopen(H_file, "r");
	while (!feof(f))
	{
		if (i >= Nt * Nr * max_iter) { // 防止数组越界
			fprintf(stderr, "H_file数据超过预期大小\n");
			break;
		}
		fscanf(f, "%f %f\n", &real_temp, &imag_temp);
		input_H[i].real = real_temp; input_H[i].imag = imag_temp;
		i++;
	}
	fclose(f); i = 0;

	f = fopen(y_file, "r");
	while (!feof(f))
	{
		if (i >= Nr * max_iter) { // 防止数组越界
			fprintf(stderr, "y_file数据超过预期大小\n");
			break;
		}
		fscanf(f, "%f %f\n", &real_temp, &imag_temp);
		input_y[i].real = real_temp; input_y[i].imag = imag_temp;
		i++;
	}
	fclose(f); i = 0;

	f = fopen(bits_file, "r");
	while (!feof(f))
	{
		if (i >= Nt * mu * max_iter) { // 防止数组越界
			fprintf(stderr, "origin_bits数据超过预期大小\n");
			break;
		}
		fscanf(f, "%d\n", &b);
		origin_bits[i] = b;
		i++;
	}
	fclose(f); i = 0;

	f = fopen(bits_output_file, "w");
	fclose(f);
	printf("\rSNR=%f loading completed...\n", SNR);

	/*开始检测*/
	for (i = 0; i < max_iter; i++)
	{
		/*将读取的数据存入变量并处理*/
		for (j = 0; j < Nr * Nt; j++)
			H[j] = input_H[Nr * Nt * i + j];
		for (j = 0; j < Nr; j++)
			y[j] = input_y[Nr * i + j];
		for (j = 0; j < Nt * mu; j++)
			bits[j] = origin_bits[Nt * mu * i + j];
		/*MIMO检测，检测类型可在main.c的MIMO_sys中更改，可选MHGD与MMSE*/
		float MSE = MIMO_detect(MIMO_sys, H, y, x_hat, SNR);
		/*解调，检测的结果比特存储在bits_demod中*/
		QAM_Demodulation(x_hat, Nt, mu, bits_demod);
		/*将解调比特结果输出到相应的文件中*/
		f = fopen(bits_output_file, "a");
		for (j = 0; j < Nt * mu; j++)
			fprintf(f, "%d\n", bits_demod[j]);
		fclose(f);
		/*将结果输出到控制台*/
		int err_bits = 0;
		err_bits = unequal_times(bits_demod, bits, Nt * mu); // 计算误比特数
		total_error_bits += err_bits; total_bits += mu * Nt;
		count++;
		printf("-------------error bits:%d, total err_bits:%d, round:%d-----------\r", err_bits, total_error_bits, i + 1);
	}
	BER = (float)total_error_bits / (float)total_bits ;
	/*av_time是每一次检测的平均用时，单位为s*/
	printf("SNR = %.2f, BER = %.8f, detect_time = %.4es, av_time = %.4es\n", SNR, BER, detect_time, detect_time / (double)i);
	Throughput[((int)SNR - MIMO_sys.SNR_start) / MIMO_sys.SNR_step] = (double)count / detect_time;
	printf("\n");

	/*返回本次SNR下的误码率，方便查看*/
	return BER;
}

/*
 * @brief  MIMO检测函数
 * @note   检测结果存放在x_hat中,同时返回MSE值
 * @param
 * @retval MSE
 */
float MIMO_detect(MIMO_sys_t MIMO_sys, MyComplex* H, MyComplex* y, MyComplex x_hat[8], float SNR)
{
	MyComplex* v_tb;

	v_tb = (MyComplex*)malloc(MIMO_sys.Nt*MIMO_sys.samples*sizeof(MyComplex));
	read_gaussian_data("/home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt", v_tb, MIMO_sys.Nt*MIMO_sys.samples, 0);
	/*接收参数*/
	int mu = MIMO_sys.mu; int Nt = MIMO_sys.Nt; int Nr = MIMO_sys.Nr;
	int samples = MIMO_sys.samples; int sampler = MIMO_sys.sampler;
	detect_type_t detect_type = MIMO_sys.detect_type;
	int EP_iter = MIMO_sys.EP_iter;
	uint8_t mmse_init = MIMO_sys.mmse_init;
	uint8_t lr_approx = MIMO_sys.lr_approx;
	int i = 0;
	float MSE = 0;
	float signal_power = 0;
	float sigma2 = 0;
	// Myreal x_hat_real[8];
	// Myimage x_hat_imag[8];
	signal_power = (float)Nt / (float)Nr;
	sigma2 = signal_power * pow(10.0f, -SNR / 10.0f);

	switch (detect_type)
	{
	/*------------------MMSE detect----------------*/
	case MMSE:
		start = clock();
		//printf("not mmse");
		end = clock();
		detect_time += (end - start) / (float)CLOCKS_PER_SEC;
		break;
	/*------------------MHGD detect----------------*/
	case MHGD:
		start = clock();
#ifdef accel
		MHGD_detect_accel(x_hat_real, x_hat_imag, Nt, Nr, mu, H, y, sigma2, mmse_init, lr_approx, samples, v_tb);
		for(int i=0; i<8; i++){
			x_hat[i].real = *(x_hat_real+i);
			x_hat[i].imag = *(x_hat_imag+i);
		}
#else
		MSE = MHGD_detect(x_hat, Nt, Nr, mu, x, H, y, sigma2, 1, 0, samples);
#endif // accel
		end = clock();
		detect_time += (end - start) / (float)CLOCKS_PER_SEC;
		break;
	default:
		break;
	}

	return MSE;
}

/*
 * @brief  解调函数
 * @note   将检测结果x_hat映射到离散的星座点上，解调结果保存在bits_demod中
 * @param  检测结果x_hat，发射天线数Nt，调制阶数mu
 * @retval none
 */
void QAM_Demodulation(MyComplex* x_hat, int Nt, int mu, int* bits_demod)
{
	switch (mu)
	{
	case 2:
		QPSK_Demodulation(x_hat, Nt, bits_demod);
		break;
	case 4:
		_16QAM_Demodulation(x_hat, Nt, bits_demod);
		break;
	case 6:
		_64QAM_Demodulation(x_hat, Nt, bits_demod);
		break;
	default:
		break;
	}
}

/*QPSK解调*/
void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex* temp = (MyComplex*)malloc(4 * sizeof(MyComplex));
	float distance[4] = { 0 };
	/*for one in x_hat*/
	for (i = 0; i < Nt; i++)
	{
		for (j = 0; j < 4; j++)
		{
			temp[j].real = x_hat[i].real * (like_float)hls::sqrt(2.0) - QPSK_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)hls::sqrt(2.0) - QPSK_Constellation[j].imag;
		}
		for (j = 0; j < 4; j++)
		{
			distance[j] = std::sqrt((float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag));
		}
		best_id = argmin(distance, 4);
		switch (best_id)
		{
		case 0:
			bits_demod[2 * i] = 0; bits_demod[2 * i + 1] = 0;
			break;
		case 1:
			bits_demod[2 * i] = 0; bits_demod[2 * i + 1] = 1;
			break;
		case 2:
			bits_demod[2 * i] = 1; bits_demod[2 * i + 1] = 0;
			break;
		case 3:
			bits_demod[2 * i] = 1; bits_demod[2 * i + 1] = 1;
			break;
		default:
			break;
		}
	}

	free(temp);
}

/*16QAM解调*/
void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex* temp = (MyComplex*)malloc(16 * sizeof(MyComplex));
	float distance[16] = { 0 };
	/*for one in x_hat*/
	for (i = 0; i < Nt; i++)
	{
		for (j = 0; j < 16; j++)
		{
			temp[j].real = x_hat[i].real * (like_float)hls::sqrt(10.0) - _16QAM_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)hls::sqrt(10.0) - _16QAM_Constellation[j].imag;
		}
		/*find the closest one*/
		for (j = 0; j < 16; j++)
		{
			distance[j] = std::sqrt((float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag));
		}
		best_id = argmin(distance, 16);
		bits_demod[4 * i] = (best_id >> 3) & 0x01;
		bits_demod[4 * i + 1] = (best_id >> 2) & 0x01;
		bits_demod[4 * i + 2] = (best_id >> 1) & 0x01;
		bits_demod[4 * i + 3] = best_id & 0x01;
	}

	free(temp);
}

/*64QAM解调*/
void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex* temp = (MyComplex*)malloc(64 * sizeof(MyComplex));
	float distance[64] = { 0 };
	int int_real; int int_imag;
	/*for one in x_hat*/
	for (i = 0; i < Nt; i++)
	{
		for (j = 0; j < 64; j++)
		{
			temp[j].real = x_hat[i].real * (like_float)sqrt(42.0) - _64QAM_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)sqrt(42.0) - _64QAM_Constellation[j].imag;
			distance[j] = (float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
		}
		best_id = argmin(distance, 64);
		switch ((int)_64QAM_Constellation[best_id].real)
		{
		case -7:
			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 0;
			break;
		case -5:
			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 1;
			break;
		case -3:
			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 1;
			break;
		case -1:
			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 0;
			break;
		case 1:
			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 0;
			break;
		case 3:
			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 1;
			break;
		case 5:
			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 1;
			break;
		case 7:
			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 0;
			break;
		default:
			break;
		}
		switch ((int)_64QAM_Constellation[best_id].imag)
		{
		case -7:
			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 0;
			break;
		case -5:
			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 1;
			break;
		case -3:
			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 1;
			break;
		case -1:
			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 0;
			break;
		case 1:
			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 0;
			break;
		case 3:
			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 1;
			break;
		case 5:
			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 1;
			break;
		case 7:
			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 0;
			break;
		default:
			break;
		}
	}

	free(temp);
}