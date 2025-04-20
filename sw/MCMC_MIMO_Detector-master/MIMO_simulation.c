#include "MIMO_simulation.h"
#include "MyComplex.h"
#include "string.h"
#include <stdio.h>


typedef unsigned char uint8_t; // ????uint8_t

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
MyComplex* x_hat;
int* bits_demod;
// float* H_real;
// float* H_image;
// float* noise_real;
// float* noise_image;
MyComplex* input_H;//???H
MyComplex* input_y;//???????y
int* origin_bits;

/*??????????????????????*/
void Parameter_init(int Nt, int Nr, int mu, int iter, int samples, detect_type_t detect_type, int max_iter)
{
	int i = 0;	float real = 0; float imag = 0;
	
	/*?????????????????????*/
	switch (mu)
	{
	case 2:
		QPSK_Constellation[0].real = -1; QPSK_Constellation[0].imag = -1;
		QPSK_Constellation[1].real = -1; QPSK_Constellation[1].imag = 1;
		QPSK_Constellation[2].real = 1; QPSK_Constellation[2].imag = -1;
		QPSK_Constellation[3].real = 1; QPSK_Constellation[3].imag = 1;
		break;
	case 4:
		f = fopen("_16QAM_Constellation.txt", "r");	
		while (!feof(f))
		{
			fscanf(f, "%f %f\n", &real, &imag);
			_16QAM_Constellation[i].real = real; _16QAM_Constellation[i].imag = imag;
			i++;
		}
		fclose(f);
		break;
	case 6:
		f = fopen("_64QAM_Constellation.txt", "r");
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

	/*??????????(????????????????????????????????*/
	switch (detect_type)
	{
	case MMSE:
		MMSE_Init(Nt, Nr);
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
	/*????????????*/
	bits = (int*)malloc(sizeof(int) * Nt * mu);/*transmitted bits*/
	x = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	H = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	y = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	noise = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	x_hat = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	bits_demod = (int*)malloc(sizeof(int) * Nt * mu);//received bits after demodulation
	// H_real = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	// H_image = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	// noise_real = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	// noise_image = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	/*???м??????????*/
	input_H = (MyComplex*)malloc(max_iter * Nt * Nr * sizeof(MyComplex));
	input_y = (MyComplex*)malloc(max_iter * Nr * sizeof(MyComplex));
	origin_bits = (int*)malloc(max_iter * Nt * mu * sizeof(int));
	
	// /*????????????????????*/
	// f = fopen("./data/C_time_seperate.csv", "a");
	// time_t timep;
	// struct tm* p;
	// time(&timep);
	// p = gmtime(&timep);
	// fprintf(f, "%d-\t%d-\t%d-\t%d:\t%d:\t%d,", 1900 + p->tm_year, 1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec);
	// fprintf(f, "Nt:%d,Nr:%d,mu:%d,samples:%d\n", Nt, Nr, mu, samples);
	// fclose(f);
}

/*????????????????????????????*/
void MIMO_free(MIMO_sys_t MIMO_sys)
{
	free(bits);
	free(x);
	free(H);
	free(y);
	free(noise);
	free(x_hat);
	// free(H_real);
	// free(H_image);
	// free(noise_real);
	// free(noise_image);
	free(bits_demod);

	switch (MIMO_sys.detect_type)
	{
	case MMSE:
		MMSE_free();
		break;
	default:
		break;
	}
}


/*
 * @brief  ??????????y ???H????????????x??????????????????????У??????????SNR????????????BER
 * @note
 * @param  SNR ?????
 * @retval ???μ??????????BER
 */
float detection(MIMO_sys_t MIMO_sys, float SNR)
{
	int Nt = MIMO_sys.Nt; int Nr = MIMO_sys.Nr; int mu = MIMO_sys.mu;
	int total_error_bits = 0; int total_bits = 0; int count = 0;
	float BER = 0;
	int max_iter = MIMO_sys.max_iter;
	detect_time = 0;
	/*????????????????????д???????????*/
	char bits_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/reference_file/bits_SNR=";;
	char H_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=";
	char y_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=";
	char bits_output_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=";
	char txt[] = ".txt";
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
	/*???????????????????????????????????????????????*/
	printf("SNR=%f loading file...", SNR);
	f = fopen(H_file, "r");
	while (!feof(f))
	{
		if (i >= Nt * Nr * max_iter) { // 防止数组越界
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

	/*??????*/
	for (i = 0; i < max_iter; i++)
	{
		/*???????????????????????*/
		for (j = 0; j < Nr * Nt; j++)
			H[j] = input_H[Nr * Nt * i + j];
		for (j = 0; j < Nr; j++)
			y[j] = input_y[Nr * i + j];
		for (j = 0; j < Nt * mu; j++)
			bits[j] = origin_bits[Nt * mu * i + j];
		/*MIMO?????????????main.c??MIMO_sys?и???????MHGD??MMSE*/
		float MSE = MIMO_detect(MIMO_sys, H, y, x_hat, SNR);
		/*???????????????洢??bits_demod??*/
		QAM_Demodulation(x_hat, Nt, mu, bits_demod);
		/*??????????????????????????*/
		f = fopen(bits_output_file, "a");
		for (j = 0; j < Nt * mu; j++)
			fprintf(f, "%d\n", bits_demod[j]);
		fclose(f);
		/*???????????????*/
		int err_bits = 0;
		err_bits = unequal_times(bits_demod, bits, Nt * mu); // ???????????
		total_error_bits += err_bits; total_bits += mu * Nt;
		count++;
		printf("-------------error bits:%d, total err_bits:%d, round:%d-----------\r", err_bits, total_error_bits, i + 1);
	}
	BER = (float)total_error_bits / (float)total_bits;
	/*av_time?????μ?????????????λ?s*/
	printf("SNR = %.2f, BER = %.8f, detect_time = %.4es, av_time = %.4es\n", SNR, BER, detect_time, detect_time / (double)i);
	Throughput[((int)SNR - MIMO_sys.SNR_start) / MIMO_sys.SNR_step] = (double)count / detect_time;
	printf("\n");

	/*???????SNR???????????????*/
	return BER;
}

/*
 * @brief  MIMO?????
 * @note   ??????????x_hat??,??????MSE?
 * @param
 * @retval MSE
 */
float MIMO_detect(MIMO_sys_t MIMO_sys, MyComplex* H, MyComplex* y, MyComplex* x_hat, float SNR)
{
	/*???????*/
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
	signal_power = (float)Nt / (float)Nr;
	sigma2 = signal_power * pow(10.0f, -SNR / 10.0f);

	switch (detect_type)
	{
	/*------------------MMSE detect----------------*/
	case MMSE:
		start = clock();
		MSE = MMSE_detect(Nt, Nr, y, H, x_hat, sigma2);	
		end = clock();
		detect_time += (end - start) / (float)CLOCKS_PER_SEC;
		break;
	/*------------------MHGD detect----------------*/
	case MHGD:
		start = clock();
#ifdef accel
		MSE = MHGD_detect_accel(x_hat, Nt, Nr, mu, H, y, sigma2, mmse_init, lr_approx, samples);
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
 * @brief  ???????
 * @note   ???????x_hat?????????????????????????????bits_demod??
 * @param  ?????x_hat????????????Nt?????????mu
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

/*QPSK???*/
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
			temp[j].real = x_hat[i].real * sqrt(2.0f) - QPSK_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * sqrt(2.0f) - QPSK_Constellation[j].imag;
		}
		for (j = 0; j < 4; j++)
		{
			distance[j] = sqrt(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
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

/*16QAM???*/
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
			temp[j].real = x_hat[i].real * sqrt(10.0f) - _16QAM_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * sqrt(10.0f) - _16QAM_Constellation[j].imag;
		}
		/*find the closest one*/
		for (j = 0; j < 16; j++)
		{
			distance[j] = sqrt(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
		}
		best_id = argmin(distance, 16);
		bits_demod[4 * i] = (best_id >> 3) & 0x01;
		bits_demod[4 * i + 1] = (best_id >> 2) & 0x01;
		bits_demod[4 * i + 2] = (best_id >> 1) & 0x01;
		bits_demod[4 * i + 3] = best_id & 0x01;
	}

	free(temp);
}

/*64QAM???*/
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
			temp[j].real = x_hat[i].real * sqrt(42.0f) - _64QAM_Constellation[j].real;
			temp[j].imag = x_hat[i].imag * sqrt(42.0f) - _64QAM_Constellation[j].imag;
			distance[j] = temp[j].real * temp[j].real + temp[j].imag * temp[j].imag;
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