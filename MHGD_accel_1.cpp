#include "MHGD_accel_1.h"
#include "MyComplex_1.h"
#include "hls_math.h"
#include <time.h>
#include "util_1.h"
//#include "hls_stream.h"



#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

//#define  __SYNTHESIS__

//MyComplex QPSK_Constellation[4];
//MyComplex _16QAM_Constellation[16];
//MyComplex _64QAM_Constellation[64];

MyComplex* HH_H;
MyComplex* H_inv;
MyComplex* grad_preconditioner;
MyComplex* constellation_norm;
MyComplex* x_mmse;
float* col_norm;
float* row_norm;
MyComplex* covar;
MyComplex* pmat;
MyComplex* x_survivor;
MyComplex* r;
MyComplex* pr_prev;
MyComplex* z_grad;
MyComplex* v;
MyComplex* z_prop;
MyComplex* x_prop;
MyComplex* r_prop;
float* temp_2mu;
MyComplex* temp_NtNt;
MyComplex* temp_NtNr;
MyComplex* temp_Nt2mu;
MyComplex* temp_1;
MyComplex* _temp_1;
MyComplex* temp_Nt;
MyComplex* temp_Nr;

MyComplex* one_Nt2mu;
MyComplex* sigma2eye;

MyComplex* v_temp;
float* p_uni_temp;
//float* p_uni;
MyComplex* v_real;
MyComplex* v_image;
int* x_init;
MyComplex** table_real;
MyComplex** table_imag;
MyComplex** temp_real;
MyComplex** temp_imag;

float* MHGD_start;
float* MHGD_end;
float* MHGD_dur;

time_seperate_t time_seperate;
unsigned int* seed;


/*
 * @brief  ��ʼ��MHGD�������Ҫ������
 * @note   
 * @param
 * @retval MSE
 */
void MHGD_Init_accel(int Nt, int Nr, int mu, int iter, MyComplex* constellation)
{
	int i;
	// ʹ�ñ�׼�� malloc ��� mkl_malloc�����Ƴ��� 64 �ֽڶ�����صĲ���
	HH_H = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	H_inv = (MyComplex*)malloc(Nr * Nt * sizeof(MyComplex));
	grad_preconditioner = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	constellation_norm = (MyComplex*)malloc(pow(2, mu) * sizeof(MyComplex));
	x_mmse = (MyComplex*)malloc(Nt * 1 * sizeof(MyComplex));
	col_norm = (float*)malloc(Nt * sizeof(float));
	row_norm = (float*)malloc(Nt * sizeof(float));
	covar = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	pmat = (MyComplex*)malloc(Nr * Nr * sizeof(MyComplex));
	x_survivor = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	r = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	pr_prev = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	z_grad = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	v = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	z_prop = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	x_prop = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	r_prop = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	temp_2mu = (float*)malloc(pow(2, mu) * sizeof(float));
	temp_NtNt = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	temp_NtNr = (MyComplex*)malloc(Nt * Nr * sizeof(MyComplex));
	temp_Nt2mu = (MyComplex*)malloc(Nt * pow(2, mu) * sizeof(MyComplex));
	temp_1 = (MyComplex*)malloc(sizeof(MyComplex));
	_temp_1 = (MyComplex*)malloc(sizeof(MyComplex));
	temp_Nt = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	temp_Nr = (MyComplex*)malloc(Nr * sizeof(MyComplex));

	one_Nt2mu = (MyComplex*)malloc(Nt * pow(2, mu) * sizeof(MyComplex));
	sigma2eye = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));

	v_temp = (MyComplex*)malloc(Nt * iter * sizeof(MyComplex));
	p_uni_temp = (float*)malloc(iter * sizeof(float));
	//p_uni = (float*)malloc(10 * sizeof(float));
	v_real = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	v_image = (MyComplex*)malloc(Nt * sizeof(MyComplex));

	x_init = (int*)malloc(Nt * sizeof(int));

	table_real = (MyComplex**)malloc(mu * sizeof(MyComplex*));
	for (i = 0; i < mu; i++)
		table_real[i] = (MyComplex*)malloc(Nr * Nt * sizeof(MyComplex));
	table_imag = (MyComplex**)malloc(mu * sizeof(MyComplex*));
	for (i = 0; i < mu; i++)
		table_imag[i] = (MyComplex*)malloc(Nr * Nt * sizeof(MyComplex));

	temp_real = (MyComplex**)malloc(Nt * sizeof(MyComplex*));
	for (i = 0; i < Nt; i++)
		temp_real[i] = (MyComplex*)malloc(Nr * sizeof(MyComplex));
	temp_imag = (MyComplex**)malloc(Nt * sizeof(MyComplex*));
	for (i = 0; i < Nt; i++)
		temp_imag[i] = (MyComplex*)malloc(Nr * sizeof(MyComplex));

	MHGD_start = (float*)malloc(time_size * sizeof(float));
	MHGD_end = (float*)malloc(time_size * sizeof(float));
	MHGD_dur = (float*)malloc(time_size * sizeof(float));

	float dqam = sqrt(1.5f / (pow(2.0f, (float)mu) - 1));
	my_complex_copy(pow(2, mu), constellation, 1, constellation_norm, 1);
	my_complex_scal(pow(2, mu), dqam, constellation_norm, 1);
	
#ifdef print_debug
	printf("constellation_norm: \n"); c_print_matrix(constellation_norm, pow(2, mu), 1);
#endif

	c_ones_generate(one_Nt2mu, Nt, pow(2, mu));
}

/*
 * @brief  �ͷ�ָ��ռ�
 * @note   
 * @param
 * @retval 
 */
void MHGD_free_accel(int Nt, int Nr, int mu)
{
	int i;
	free(HH_H);
	free(H_inv);
	free(grad_preconditioner);
	free(constellation_norm);
	free(x_mmse);
	free(col_norm);
	free(row_norm);
	free(covar);
	free(pmat);
	free(x_survivor);
	free(r);
	free(pr_prev);
	free(z_grad);
	free(v);
	free(z_prop);
	free(x_prop);
	free(r_prop);
	free(temp_2mu);
	free(temp_NtNt);
	free(temp_NtNr);
	free(temp_Nt2mu);
	free(temp_1);
	free(_temp_1);
	free(temp_Nt);
	free(temp_Nr);

	free(one_Nt2mu);
	free(sigma2eye);

	free(v_temp);
	free(p_uni_temp);
	//free(p_uni);

	free(v_real);
	free(v_image);

	free(x_init);

	for (i = 0; i < mu; i++)
		free(table_real[i]);
	for (i = 0; i < mu; i++)
		free(table_imag[i]);
	for (i = 0; i < Nt; i++)
		free(temp_real[i]);
	for (i = 0; i < Nt; i++)
		free(temp_imag[i]);

	free(table_real);
	free(table_imag);
	free(temp_real);
	free(temp_imag);

	free(seed);
}

/*
 * @brief  MIMO��⺯��
 * @note   ����������x_hat��,ͬʱ����MSEֵ 
 *		   ��ԭʼ��MHGD_detect��ȣ�����һ�����Ż���ʵ���ٶȸ��죬��˵���ҪMHGD�㷨���ʱͳһ���øú���
 * @param
 * @retval MSE
 */
void MHGD_detect_accel(Myreal* x_hat_real, Myimage* x_hat_imag,int Nt, int Nr, int mu, MyComplex* H, MyComplex* y, float sigma2, int mmse_init, int lr_approx, int iter, MyComplex* v_tb)
{
	#pragma HLS interface mode=m_axi port = x_hat_real bundle=X depth = 8
	#pragma HLS interface mode=m_axi port = x_hat_imag bundle=X depth = 8
	#pragma HLS interface mode=s_axilite port = Nt bundle=contrl
	#pragma HLS interface mode=s_axilite port = Nr bundle=contrl
	#pragma HLS interface mode=s_axilite port = mu bundle=contrl
	#pragma HLS interface mode=m_axi port = H bundle=H depth = 64
	#pragma HLS interface mode=m_axi port = y bundle=Y depth = 8
	#pragma HLS interface mode=s_axilite port = sigma2 bundle=contrl
	#pragma HLS interface mode=s_axilite port = mmse_init bundle=contrl
	#pragma HLS interface mode=s_axilite port = lr_approx bundle=contrl
	#pragma HLS interface mode=s_axilite port = iter bundle=contrl
	#pragma HLS interface mode=m_axi port = v_tb bundle=V_TB depth = 256

    #pragma HLS data_pack variable=x_hat_real  // 确保结构体连续存储
	#pragma HLS data_pack variable=x_hat_imag  // 确保结构体连续存储
    //  #pragma HLS DATA_PACK variable=x_hat_real
	//  #pragma HLS DATA_PACK variable=x_hat_imag
	//  #pragma HLS DATA_PACK variable=y
	//  #pragma HLS DATA_PACK variable=H
	//  #pragma HLS DATA_PACK variable=v_tb
	MyComplex x_hat_1[8];
	// for(int i= 0;i<8;i++){
		// #pragma HLS PIPELINE II=1
		// x_hat_1[i].real=*(x_hat_real + i);
		// x_hat_1[i].imag=*(x_hat_imag + i);
	// }

	
	/*��������*/
	like_float alpha;
	like_float dqam;
	like_float step_size;
	int i; int j; int k; int id;
	like_float log_pacc;
	like_float r_norm;/*the norm of residual vector*/
	like_float r_norm_survivor;
	like_float lr;/*learning rate*/
	like_float r_norm_prop;/*new norm*/
	like_float p_acc;
	like_float p_uni[9];

	like_float temp_3;

	int transA = 1;  // CblasConjTrans �ĵ�Чֵ����ʾ����ת��
	int transB = 0;  // CblasNoTrans �ĵ�Чֵ����ʾ��ת��

	///float* MHGD_start;
	///MHGD_start = (float*)malloc(time_size * sizeof(float));

    //transA = 1;
	#ifndef __SYNTHESIS__

	MHGD_start[(int)gradpre] = clock();

	#endif

	like_float numerator = 1.5;
    like_float denominator = hls::pow((ap_fixed<64,32>)2, (ap_fixed<64,32>)mu) - (ap_fixed<64,32>)1;
	like_float tttt = numerator/denominator;
	dqam = hls::sqrt((ap_fixed<64,32>)tttt);/*���ͷ���֮����С�����һ�룬�������㾭����һ��������Ľ��*/
	/*----------------------------------------------------------------*/
	c_eye_generate(sigma2eye, Nt, sigma2 / (float)(dqam * dqam));
	/*�����ݶ��½�������grad_preconditioner*/
	c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, HH_H);
	my_complex_add(Nt * Nt, HH_H, sigma2eye, grad_preconditioner);
	Inverse_LU(grad_preconditioner, Nt, Nt);

	//alpha = (like_float)1.0 / hls::pow((like_float)Nt / (like_float)8.0, (like_float)(1.0 / 3.0));
	like_float exponent = like_float(1) / like_float(3); // 避免浮点字面值隐式转换
	like_float exponent_1 = (like_float)Nt / like_float(8);
    alpha = like_float(1) / hls::pow<64,32>((ap_fixed<64,32>)exponent, (ap_fixed<64,32>)exponent_1);

	#ifndef __SYNTHESIS__

	MHGD_end[(int)gradpre] = clock();

	/*Э�������*/
	MHGD_start[covar_cal] = clock();
	#endif
	c_eye_generate(covar, Nt, 1.0f);
	#ifndef __SYNTHESIS__

	MHGD_end[covar_cal] = clock();

	MHGD_start[lr_cal] = clock();
	#endif
	/*For learning rate line search */
	if (!lr_approx)
	{
		c_matmultiple(H, transB, grad_preconditioner, transB, Nr, Nt, Nt, Nt, temp_NtNr);
		c_matmultiple(temp_NtNr, transB, H, transA, Nr, Nt, Nr, Nt, pmat);
	}
	else
	{
		for (i = 0; i < Nr * Nr; i++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nr_2_max min=Nr_2_min
			pmat[i].real = (like_float)0; pmat[i].imag = (like_float)0;
		}
	}
	#ifndef __SYNTHESIS__
	MHGD_end[lr_cal] = clock();

	/*initialize the estimate with size (np, nt, 1), np is for parallel samplers*/
	MHGD_start[mmse_and_r_init] = clock();
	#endif
	if (mmse_init)
	{
		/*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
		c_eye_generate(sigma2eye, Nt, sigma2);
		my_complex_add(Nt * Nt, sigma2eye, HH_H, temp_NtNt);
		Inverse_LU(temp_NtNt, Nt, Nt);
		c_matmultiple(H, transA, y, transB, Nr, Nt, Nr, transA, temp_Nt);
		c_matmultiple(temp_NtNt, transB, temp_Nt, transB, Nt, Nt, Nt, transA, x_mmse);
		/*ӳ�䵽��һ������ͼ�� xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map(mu, Nt, dqam, x_mmse, x_hat_1);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		generateUniformRandoms_int(Nt, x_init, mu);
		for (i = 0; i < Nt; i++)
		{   
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			x_hat_1[i] = constellation_norm[x_init[i]];
		}
	}
	/*����ʣ������r=y-Hx*/
	c_matmultiple(H, transB, x_hat_1, transB, Nr, Nt, Nt, transA, r);
	my_complex_sub(Nr, y, r, r);

	/*����ʣ�������ķ���������ģֵ��*/
	c_matmultiple(r, transA, r, transB, Nr, transA, Nr, transA, temp_1);
	r_norm = temp_1[0].real;
	my_complex_copy(Nt, x_hat_1, 1, x_survivor, 1);
	r_norm_survivor = r_norm;
	/*ȷ������ѧϰ��*/
	if (!lr_approx)
	{
		c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, transA, pr_prev);
		c_matmultiple(r, transA, pr_prev, transB, Nr, transA, Nr, transA, temp_1);
		c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, transA, Nr, transA, _temp_1);
		lr = temp_1[0].real / _temp_1[0].real;
	}
	else
	{
		lr = 1;
	}

	step_size = max(dqam, like_float(hls::sqrt(r_norm / (like_float)Nr))) * alpha;
	#ifndef __SYNTHESIS__
	MHGD_end[mmse_and_r_init] = clock();
	#endif
	/*================core==============*/
	int offset = 0;  // ���ļ��Ŀ�ʼλ�ÿ�ʼ��ȡ
	for (k = 0; k < iter; k++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=iter_max min=iter_min
	  	#pragma HLS PIPELINE II=1
		/*�����ݶ� z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
		#ifndef __SYNTHESIS__
		MHGD_start[zgrad] = clock();
		#endif
		c_matmultiple(H, transA, r, transB, Nr, Nt, Nr, transA, temp_Nt);
		c_matmultiple(grad_preconditioner, transB, temp_Nt, transB, Nt, Nt, Nt, transA, z_grad);
		my_complex_scal(Nt, lr, z_grad, 1); 
		my_complex_add(Nt, x_hat_1, z_grad, z_grad);
		#ifndef __SYNTHESIS__
		MHGD_end[zgrad] = clock();
		/*�����˹����Ŷ�*/
		MHGD_start[randomwalk] = clock();
		#endif
		//read_gaussian_data("/home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt", v, Nt, offset);
		for(int i = 0; i < Nt; i++){
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			#pragma HLS UNROLL
			v[i].real = v_tb[i+offset].real;
			v[i].imag = v_tb[i+offset].imag;
		}
		offset = offset + Nt;
		my_complex_scal(Nt, step_size, v, 1);
		my_complex_add(Nt, z_grad, v, z_prop);
		#ifndef __SYNTHESIS__
		MHGD_end[randomwalk] = clock();
		/*���ݶ�ӳ�䵽QAM�������� x_prop = constellation_norm[np.argmin(abs(z_prop * ones - constellation_norm), axis=2)].reshape(-1, nt, 1) */
		MHGD_start[map_op] = clock();
		#endif
		map(mu, Nt, dqam, z_prop, x_prop);
		#ifndef __SYNTHESIS__
		MHGD_end[map_op] = clock();
		

		/*�����µĲв�� calculate residual norm of the proposal*/
		MHGD_start[new_r] = clock();
		#endif
		c_matmultiple(H, transB, x_prop, transB, Nr, Nt, Nt, transA, temp_Nr);
		my_complex_sub(Nr, y, temp_Nr, r_prop);
		c_matmultiple(r_prop, transA, r_prop, transB, Nr, transA, Nr, transA, temp_1);
		r_norm_prop = temp_1[0].real;

		/*update the survivor*/
		if (r_norm_survivor > r_norm_prop)
		{
			my_complex_copy(Nt, x_prop, 1, x_survivor, 1);
			r_norm_survivor = r_norm_prop;
		}
		#ifndef __SYNTHESIS__
		MHGD_end[new_r] = clock();

		/*acceptance test*/
		MHGD_start[accept] = clock();
		#endif
		temp_3 = -(r_norm_prop - r_norm);
		log_pacc = min(like_float(0), temp_3);
		p_acc = hls::exp<64,32>((ap_fixed<64,32>)log_pacc);
		generateUniformRandoms_float(10, p_uni);
		if (p_acc > p_uni[5])/*������������ʱ��*/
		{
			my_complex_copy(Nt, x_prop, 1, x_hat_1, 1);
			my_complex_copy(Nt, r_prop, 1, r, 1);
			r_norm = r_norm_prop;
			/*update GD learning rate*/
			if (!lr_approx)
			{
				c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, transA, pr_prev);
				c_matmultiple(r, transA, pr_prev, transB, Nr, transA, Nr, transA, temp_1);
				c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, transA, Nr, transA, _temp_1);
				lr = temp_1[0].real / _temp_1[0].real;
			}
			/*update random walk size*/
			step_size = max(dqam, (like_float)hls::sqrt(r_norm / (like_float)Nr)) * alpha;
		}
		#ifndef __SYNTHESIS__
		MHGD_end[accept] = clock();

		for (i = zgrad; i <= accept; i++)
			MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
		#endif
	}
	/*����ʵ��ʱ��*/
	#ifndef __SYNTHESIS__
	MHGD_start[x_copy] = clock();
	#endif
	my_complex_copy(Nt, x_survivor, 1, x_hat_1, 1);
	#ifndef __SYNTHESIS__
	MHGD_end[x_copy] = clock();
	for (i = 0; i <= mmse_and_r_init; i++)
		MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
	MHGD_dur[x_copy] += (MHGD_end[x_copy] - MHGD_start[x_copy]) / (float)CLOCKS_PER_SEC;
	#endif
	for(int i= 0;i<8;i++){
		#pragma HLS UNROLL
		*(x_hat_real+i)=x_hat_1[i].real;
		*(x_hat_imag+i)=x_hat_1[i].imag;
	}

	//dqam = (like_float)hls::sqrt((like_float)1.5 / (like_float)((like_float)hls::pow((like_float)2.0, (like_float)mu) - (like_float)1));/*���ͷ���֮����С�����һ�룬�������㾭����һ��������Ľ��*/
	///*----------------------------------------------------------------*/
	//c_eye_generate(sigma2eye, Nt, sigma2 / (float)(dqam * dqam));
	///*�����ݶ��½�������grad_preconditioner*/
	//c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, HH_H);
	//my_complex_add(Nt * Nt, HH_H, sigma2eye, grad_preconditioner);
	//Inverse_LU(grad_preconditioner, Nt, Nt);
//
	////alpha = (like_float)1.0 / hls::pow((like_float)Nt / (like_float)8.0, (like_float)(1.0 / 3.0));
	//like_float exponent = like_float(1) / like_float(3); // 避免浮点字面值隐式转换
	//like_float exponent_1 = (like_float)Nt / like_float(8);
    //alpha = like_float(1) / hls::pow(exponent_1, exponent);
//
	///*协方差矩阵*/
	//c_eye_generate(covar, Nt, 1.0f);
//
	///*For learning rate line search */
	//if (!lr_approx)
	//{
	//	c_matmultiple(H, transB, grad_preconditioner, transB, Nr, Nt, Nt, Nt, temp_NtNr);
	//	c_matmultiple(temp_NtNr, transB, H, transA, Nr, Nt, Nr, Nt, pmat);
	//}
	//else
	//{
	//	for (i = 0; i < Nr * Nr; i++)
	//	{
	//		#pragma HLS LOOP_TRIPCOUNT max=Nr_2_max min=Nr_2_min
	//		pmat[i].real = (like_float)0; pmat[i].imag = (like_float)0;
	//	}
	//}
//
	///*initialize the estimate with size (np, nt, 1), np is for parallel samplers*/
	//if (mmse_init)
	//{
	//	/*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
	//	c_eye_generate(sigma2eye, Nt, sigma2);
	//	my_complex_add(Nt * Nt, sigma2eye, HH_H, temp_NtNt);
	//	Inverse_LU(temp_NtNt, Nt, Nt);
	//	c_matmultiple(H, transA, y, transB, Nr, Nt, Nr, 1, temp_Nt);
	//	c_matmultiple(temp_NtNt, transB, temp_Nt, transB, Nt, Nt, Nt, 1, x_mmse);
	//	/*ӳ�䵽��һ������ͼ�� xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
	//	map(mu, Nt, dqam, x_mmse, x_hat);
	//}
	//else
	//{
	//	/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
	//	generateUniformRandoms_int(Nt, x_init, mu);
	//	for (i = 0; i < Nt; i++)
	//	{   
	//		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
	//		x_hat[i] = constellation_norm[x_init[i]];
	//	}
	//}
	///*����ʣ������r=y-Hx*/
	//c_matmultiple(H, transB, x_hat, transB, Nr, Nt, Nt, 1, r);
	//my_complex_sub(Nr, y, r, r);
//
	///*����ʣ�������ķ���������ģֵ��*/
	//c_matmultiple(r, transA, r, transB, Nr, 1, Nr, 1, temp_1);
	//r_norm = temp_1[0].real;
	//my_complex_copy(Nt, x_hat, 1, x_survivor, 1);
	//r_norm_survivor = r_norm;
	///*ȷ������ѧϰ��*/
	//if (!lr_approx)
	//{
	//	c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, 1, pr_prev);
	//	c_matmultiple(r, transA, pr_prev, transB, Nr, 1, Nr, 1, temp_1);
	//	c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, 1, Nr, 1, _temp_1);
	//	lr = temp_1[0].real / _temp_1[0].real;
	//}
	//else
	//{
	//	lr = 1;
	//}
//
	//step_size = max(dqam, like_float(hls::sqrt(r_norm / (like_float)Nr))) * alpha;
	///*================core==============*/
	//int offset = 0;  // ���ļ��Ŀ�ʼλ�ÿ�ʼ��ȡ
//
	////uint32_t seed;
//
	//for (k = 0; k < iter; k++)
	//{
	//	#pragma HLS LOOP_TRIPCOUNT max=iter_max min=iter_min
	//	#pragma HLS PIPELINE II=1
	//	/*�����ݶ� z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
	//	c_matmultiple(H, transA, r, transB, Nr, Nt, Nr, 1, temp_Nt);
	//	c_matmultiple(grad_preconditioner, transB, temp_Nt, transB, Nt, Nt, Nt, 1, z_grad);
	//	my_complex_scal(Nt, lr, z_grad, 1); 
	//	my_complex_add(Nt, x_hat, z_grad, z_grad);
	//	/*�����˹����Ŷ�*/
	//	//read_gaussian_data("/home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt", v, Nt, offset);
	//	for(int i = 0; i < Nt; i++){
	//		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
	//		v[i].real = v_tb[i+offset].real;
	//		v[i].imag = v_tb[i+offset].imag;
	//	}
	//	offset = offset + Nt;
	//	my_complex_scal(Nt, step_size, v, 1);
	//	my_complex_add(Nt, z_grad, v, z_prop);
	//	/*���ݶ�ӳ�䵽QAM�������� x_prop = constellation_norm[np.argmin(abs(z_prop * ones - constellation_norm), axis=2)].reshape(-1, nt, 1) */
	//	map(mu, Nt, dqam, z_prop, x_prop);
//
	//	/*�����µĲв�� calculate residual norm of the proposal*/
//
	//	c_matmultiple(H, transB, x_prop, transB, Nr, Nt, Nt, 1, temp_Nr);
	//	my_complex_sub(Nr, y, temp_Nr, r_prop);
	//	c_matmultiple(r_prop, transA, r_prop, transB, Nr, 1, Nr, 1, temp_1);
	//	r_norm_prop = temp_1[0].real;
//
	//	/*update the survivor*/
	//	if (r_norm_survivor > r_norm_prop)
	//	{
	//		my_complex_copy(Nt, x_prop, 1, x_survivor, 1);
	//		r_norm_survivor = r_norm_prop;
	//	}
//
	//	/*acceptance test*/
	//	temp_3 = -(r_norm_prop - r_norm);
	//	log_pacc = min(like_float(0), temp_3);
	//	p_acc = (like_float)hls::exp(log_pacc);
	//	generateUniformRandoms_float(10, p_uni);//该函数给p_uni值
//
	//	if (p_acc > p_uni[5])/*������������ʱ��*/
	//	{
	//		my_complex_copy(Nt, x_prop, 1, x_hat, 1);
	//		my_complex_copy(Nt, r_prop, 1, r, 1);
	//		r_norm = r_norm_prop;
	//		/*update GD learning rate*/
	//		if (!lr_approx)
	//		{
	//			c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, 1, pr_prev);
	//			c_matmultiple(r, transA, pr_prev, transB, Nr, 1, Nr, 1, temp_1);
	//			c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, 1, Nr, 1, _temp_1);
	//			lr = temp_1[0].real / _temp_1[0].real;
	//		}
	//		/*update random walk size*/
	//		step_size = max(dqam, (like_float)hls::sqrt(r_norm / (like_float)Nr)) * alpha;
	//	}
	//}
	///*����ʵ��ʱ��*/
	//my_complex_copy(Nt, x_survivor, 1, x_hat, 1);
	////#endif

	//return 0;
}


void Hx_accel(int Nt, int Nr, int mu, float dqam, MyComplex* x, MyComplex* x_result)
{
	like_float dqam_1 = dqam;
	my_complex_scal(Nr, (like_float)0, x_result, 1);
	int i, j;
	int Hx_table_id;
#ifdef print_debug
	printf("x:\n");
	c_print_vector(x, Nt);
#endif
	/*������Ҫ���н������temp_real��temp_imag��*/
	for (i = 0; i < Nt; i++)
	{
		Hx_table_id = (x[i].real / dqam_1 + mu - 1) / 2;
#ifdef print_debug
		printf("Hx_table_id:%d  ", Hx_table_id);
#endif
		for (j = 0; j < Nr; j++)
		{
			temp_real[i][j] = table_real[Hx_table_id][Nt * j + i];
		}
		Hx_table_id = (x[i].imag / dqam_1 + mu - 1) / 2;
#ifdef print_debug
		printf("Hx_table_id:%d", Hx_table_id);
#endif
		for (j = 0; j < Nr; j++)
		{
			temp_imag[i][j] = table_imag[Hx_table_id][Nt * j + i];
		}
	}

#ifdef print_debug
	printf("\ntemp_real:\n");
	for (i = 0; i < Nt; i++)
	{
		c_print_vector(temp_real[i], Nr);
	}
	printf("temp_imag:\n");
	for (i = 0; i < Nt; i++)
	{
		c_print_vector(temp_imag[i], Nr);
	}
#endif
	for (i = 0; i < Nt; i++)
	{
		my_complex_add(Nr, x_result, temp_real[i], x_result);
		my_complex_add(Nr, x_result, temp_imag[i], x_result);
	}
	my_complex_scal(Nr, dqam_1, x_result, 1);
}

//void initial_part

//void QAM_Demodulation(MyComplex* x_hat, int Nt, int mu, int* bits_demod)
//{
//	switch (mu)
//	{
//	case 2:
//		QPSK_Demodulation(x_hat, Nt, bits_demod);
//		break;
//	case 4:
//		_16QAM_Demodulation(x_hat, Nt, bits_demod);
//		break;
//	case 6:
//		_64QAM_Demodulation(x_hat, Nt, bits_demod);
//		break;
//	default:
//		break;
//	}
//}

//void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
//{
//	int i = 0; int j = 0; int best_id = 0;
//	MyComplex* temp = (MyComplex*)malloc(4 * sizeof(MyComplex));
//	float distance[4] = { 0 };
//	/*for one in x_hat*/
//	for (i = 0; i < Nt; i++)
//	{
//		for (j = 0; j < 4; j++)
//		{
//			temp[j].real = x_hat[i].real * sqrt(2.0f) - QPSK_Constellation[j].real;
//			temp[j].imag = x_hat[i].imag * sqrt(2.0f) - QPSK_Constellation[j].imag;
//		}
//		for (j = 0; j < 4; j++)
//		{
//			distance[j] = sqrt(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
//		}
//		best_id = argmin(distance, 4);
//		switch (best_id)
//		{
//		case 0:
//			bits_demod[2 * i] = 0; bits_demod[2 * i + 1] = 0;
//			break;
//		case 1:
//			bits_demod[2 * i] = 0; bits_demod[2 * i + 1] = 1;
//			break;
//		case 2:
//			bits_demod[2 * i] = 1; bits_demod[2 * i + 1] = 0;
//			break;
//		case 3:
//			bits_demod[2 * i] = 1; bits_demod[2 * i + 1] = 1;
//			break;
//		default:
//			break;
//		}
//	}
//
//	free(temp);
//}
//
//void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
//{
//	int i = 0; int j = 0; int best_id = 0;
//	MyComplex* temp = (MyComplex*)malloc(16 * sizeof(MyComplex));
//	float distance[16] = { 0 };
//	/*for one in x_hat*/
//	for (i = 0; i < Nt; i++)
//	{
//		for (j = 0; j < 16; j++)
//		{
//			temp[j].real = x_hat[i].real * sqrt(10.0f) - _16QAM_Constellation[j].real;
//			temp[j].imag = x_hat[i].imag * sqrt(10.0f) - _16QAM_Constellation[j].imag;
//		}
//		/*find the closest one*/
//		for (j = 0; j < 16; j++)
//		{
//			distance[j] = sqrt(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
//		}
//		best_id = argmin(distance, 16);
//		bits_demod[4 * i] = (best_id >> 3) & 0x01;
//		bits_demod[4 * i + 1] = (best_id >> 2) & 0x01;
//		bits_demod[4 * i + 2] = (best_id >> 1) & 0x01;
//		bits_demod[4 * i + 3] = best_id & 0x01;
//	}
//
//	free(temp);
//}
//
//void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
//{
//	int i = 0; int j = 0; int best_id = 0;
//	MyComplex* temp = (MyComplex*)malloc(64 * sizeof(MyComplex));
//	float distance[64] = { 0 };
//	int int_real; int int_imag;
//	/*for one in x_hat*/
//	for (i = 0; i < Nt; i++)
//	{
//		for (j = 0; j < 64; j++)
//		{
//			temp[j].real = x_hat[i].real * sqrt(42.0f) - _64QAM_Constellation[j].real;
//			temp[j].imag = x_hat[i].imag * sqrt(42.0f) - _64QAM_Constellation[j].imag;
//			distance[j] = temp[j].real * temp[j].real + temp[j].imag * temp[j].imag;
//		}
//		best_id = argmin(distance, 64);
//		switch ((int)_64QAM_Constellation[best_id].real)
//		{
//		case -7:
//			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 0;
//			break;
//		case -5:
//			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 1;
//			break;
//		case -3:
//			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 1;
//			break;
//		case -1:
//			bits_demod[6 * i] = 0; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 0;
//			break;
//		case 1:
//			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 0;
//			break;
//		case 3:
//			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 1; bits_demod[6 * i + 2] = 1;
//			break;
//		case 5:
//			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 1;
//			break;
//		case 7:
//			bits_demod[6 * i] = 1; bits_demod[6 * i + 1] = 0; bits_demod[6 * i + 2] = 0;
//			break;
//		default:
//			break;
//		}
//		switch ((int)_64QAM_Constellation[best_id].imag)
//		{
//		case -7:
//			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 0;
//			break;
//		case -5:
//			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 1;
//			break;
//		case -3:
//			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 1;
//			break;
//		case -1:
//			bits_demod[6 * i + 3] = 0; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 0;
//			break;
//		case 1:
//			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 0;
//			break;
//		case 3:
//			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 1; bits_demod[6 * i + 5] = 1;
//			break;
//		case 5:
//			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 1;
//			break;
//		case 7:
//			bits_demod[6 * i + 3] = 1; bits_demod[6 * i + 4] = 0; bits_demod[6 * i + 5] = 0;
//			break;
//		default:
//			break;
//		}
//	}
//
//	free(temp);
//};
//void initialize_parameters(
//    int mu, int Nt, float sigma2,
//    float dqam, MyComplex sigma2eye[Nt*Nt]
//) {
//    #pragma HLS INLINE
//    // 计算dqam和sigma2eye
//    dqam = sqrt(1.5f / (pow(2.0f, (float)mu) - 1));
//    c_eye_generate(sigma2eye, Nt, sigma2 / (dqam * dqam));
//};
//
//void compute_grad_preconditioner(
//    MyComplex H[Nr][Nt],  // 需要改为hls::stream<MyComplex>
//    int Nr, int Nt,
//    MyComplex HH_H[Nt*Nt],
//    MyComplex grad_preconditioner[Nt*Nt]
//) {
//    #pragma HLS INLINE
//    // 计算H^H * H
//    c_matmultiple(H, 1, H, 0, Nr, Nt, Nr, Nt, HH_H);
//    // 矩阵求逆
//    my_complex_add(Nt*Nt, HH_H, sigma2eye, grad_preconditioner);
//    Inverse_LU(grad_preconditioner, Nt, Nt);
//};