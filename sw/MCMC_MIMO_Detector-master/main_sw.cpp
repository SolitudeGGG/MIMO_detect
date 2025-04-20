#include "MHGD_accel_sw.h"
#include "MyComplex.h"
#include "hls_math.h"
#include <string.h>
#include <stdio.h>

int main()
{
    /*变量定义*/
    FILE* ff;
    int bits[Ntr_1 * mu_1];
    float SNR = 25;
    float detect_time = 0;
    int bits_demod[Ntr_1 * mu_1];
    int origin_bits[max_iter_1 * Ntr_1 * mu_1];
    int Nt = Ntr_1; int Nr = Ntr_1; int mu = mu_1;
	int total_error_bits = 0; int total_bits = 0; int count = 0;
    float BER = 0;
	int max_iter = max_iter_1;
    int i = 0; float real = 0; float imag = 0; int j = 0; int b;
    float real_temp;
	float imag_temp;
    MyComplex x[Ntr_1];
    MyComplex H[Ntr_1 * Ntr_1];
    MyComplex y[Ntr_1];
    MyComplex noise[Ntr_1];
    MyComplex x_hat[Ntr_1];
    float x_hat_real[Ntr_1];
    float x_hat_imag[Ntr_1];
    float H_real[Ntr_1 * Ntr_1];
    float H_imag[Ntr_1 * Ntr_1];
    float y_real[Ntr_1];
    float y_imag[Ntr_1];
    MyComplex v_tb[Ntr_1 * iter_1];
    float v_tb_real[Ntr_1 * iter_1];
    float v_tb_imag[Ntr_1 * iter_1];
    MyComplex noise_real[Ntr_1];
    MyComplex noise_image[Ntr_1];
    MyComplex input_H[max_iter_1 * Ntr_1 * Ntr_1];
    MyComplex input_y[max_iter_1 * Ntr_1];

    /*字符串拼接，根据信噪比不同写入不同的文本文件*/
	char bits_file[1024] = "/home/ggg_wufuqi/hls/MHGD/MHGD/reference_file/bits_SNR=";
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

    /*读取输入测试文件数据（信道数据、接受信号数据、比特数据）*/
	printf("SNR=%f loading file...", SNR);
	ff = fopen(H_file, "r");
	while (!feof(ff))
	{
		if (i >= Nt * Nr * max_iter) { // 防止数组越界
			break;
		}
		fscanf(ff, "%f %f\n", &real_temp, &imag_temp);
		input_H[i].real = real_temp; input_H[i].imag = imag_temp;
		i++;
	}
	fclose(ff); i = 0;

	ff = fopen(y_file, "r");
	while (!feof(ff))
	{
		if (i >= Nr * max_iter) { // 防止数组越界
			break;
		}
		fscanf(ff, "%f %f\n", &real_temp, &imag_temp);
		input_y[i].real = real_temp; input_y[i].imag = imag_temp;
		i++;
	}
	fclose(ff); i = 0;

	ff = fopen(bits_file, "r");
	while (!feof(ff))
	{
		if (i >= Nt * mu * max_iter) { // 防止数组越界
			break;
		}
		fscanf(ff, "%d\n", &b);
		origin_bits[i] = b;
		i++;
	}
	fclose(ff); i = 0;

	ff = fopen(bits_output_file, "w");
	fclose(ff);
	printf("\rSNR=%f loading completed...\n", SNR);

    /*开始检测*/
    for (i = 0; i < max_iter; i++)
    {
        int l;
        float MSE = 0;/*?*/
        /*将读取的数据存入变量并处理*/
		for (j = 0; j < Nr * Nt; j++)
            H[j] = input_H[Nr * Nt * i + j];
        for (j = 0; j < Nr; j++)
            y[j] = input_y[Nr * i + j];
        for (j = 0; j < Nt * mu; j++)
            bits[j] = origin_bits[Nt * mu * i + j];
        /*MIMO检测，检测类型可在main.c的MIMO_sys中更改，可选MHGD与MMSE*/
        float signal_power = 0;
	    float sigma2 = 0;
        signal_power = (float)Nt / (float)Nr;
        sigma2 = signal_power * pow(10.0f, -SNR / 10.0f);
        read_gaussian_data("/home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt", v_tb, Nt*iter_1, 0);
        for (l = 0; l < Nt*iter_1; l++){
            v_tb_real[l] = v_tb[l].real;
            v_tb_imag[l] = v_tb[l].imag;
        }
        for (l = 0; l < Nr * Nt; l++){
            H_real[l] = H[l].real;
            H_imag[l] = H[l].imag;
        }
        for (l = 0; l < Nr; l++){
            y_real[l] = y[l].real;
            y_imag[l] = y[l].imag;
        }
        MHGD_detect_accel_sw(x_hat_real, x_hat_imag, H_real, H_imag, y_real, y_imag, sigma2, v_tb_real, v_tb_imag);
        for(l = 0; l < Nt; l++){
			x_hat[l].real = x_hat_real[l];
			x_hat[l].imag = x_hat_imag[l];
		}
        /*解调，检测的结果比特存储在bits_demod中*/
        QAM_Demodulation(x_hat, Nt, mu, bits_demod);
        /*将解调比特结果输出到相应的文件中*/
        ff = fopen(bits_output_file, "a");
        for (l = 0; l < Nt * mu; l++)
            fprintf(ff, "%d\n", bits_demod[l]);
        fclose(ff);
        /*将结果输出到控制台*/
        int err_bits = 0;
        err_bits = unequal_times(bits_demod, bits, Nt * mu); // 计算误比特数
        total_error_bits += err_bits; total_bits += mu * Nt;
        count++;
        printf("-------------error bits:%d, total err_bits:%d, round:%d-----------\r", err_bits, total_error_bits, i + 1);
    }
    BER = (float)total_error_bits / (float)total_bits ;
	/*av_time是每一次检测的平均用时，单位为s*/
	printf("SNR = %.2f, BER = %.8f", SNR, BER);
	printf("\n");
    return 0;
}