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
// MyComplex* x_survivor;
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
// float* p_uni;
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
	// x_survivor = (MyComplex*)malloc(Nt * sizeof(MyComplex));
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
	// p_uni = (float*)malloc(10 * sizeof(float));
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
	// free(x_survivor);
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
	// free(p_uni);

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
void MHGD_detect_accel(
    Myreal* x_hat_real, 
    Myimage* x_hat_imag,
    int Nt, 
    int Nr, 
    int mu, 
    Myreal* H_real, 
    Myimage* H_imag, 
    Myreal* y_real, 
    Myimage* y_imag, 
    float sigma2, 
    int mmse_init, 
    int lr_approx, 
    int iter, 
    Myreal* v_tb_real, 
    Myimage* v_tb_imag
) {
    // AXI-Master 接口配置
    #pragma HLS INTERFACE mode=m_axi port=x_hat_real bundle=AXI_X depth=8 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=x_hat_imag bundle=AXI_X depth=8 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=H_real bundle=AXI_H depth=64 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=H_imag bundle=AXI_H depth=64 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=y_real bundle=AXI_Y depth=8 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=y_imag bundle=AXI_Y depth=8 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=v_tb_real bundle=AXI_V depth=256 offset=slave
    #pragma HLS INTERFACE mode=m_axi port=v_tb_imag bundle=AXI_V depth=256 offset=slave

    // AXI-Lite 控制寄存器绑定
    #pragma HLS INTERFACE mode=s_axilite port=Nt bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=Nr bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=mu bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=sigma2 bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=mmse_init bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=lr_approx bundle=CTRL
    #pragma HLS INTERFACE mode=s_axilite port=iter bundle=CTRL
// 
    // #pragma HLS DATA_PACK variable=x_hat_real
	// #pragma HLS DATA_PACK variable=x_hat_imag
	// #pragma HLS DATA_PACK variable=y_real
	// #pragma HLS DATA_PACK variable=y_imag
	// #pragma HLS DATA_PACK variable=H_real
	// #pragma HLS DATA_PACK variable=H_imag
	// #pragma HLS DATA_PACK variable=v_tb_real
	// #pragma HLS DATA_PACK variable=v_tb_imag
	
	int Nt_local = Nt;
	int Nr_local = Nr;
	int mu_local = mu;
	float sigma2_local = sigma2;
	int mmse_init_local = mmse_init;
	int lr_approx_local = lr_approx;
	int iter_local = iter;
	
	MyComplex x_hat_1[8];
	MyComplex y_local[8];
	MyComplex H_local[64];
	MyComplex v_tb_local[256];
	for(int i=0; i<8; i++){
		y_local[i].real = *(y_real+i);
		y_local[i].imag = *(y_imag+i);
	}
	for (int i = 0; i < 64; i++)
	{
		H_local[i].real = *(H_real+i);
		H_local[i].imag = *(H_imag+i);
	}
	for (int i = 0; i < 256; i++)
	{
		v_tb_local[i].real = *(v_tb_real+i);
		v_tb_local[i].imag = *(v_tb_imag+i);
	}

	//  /*变量定义*/
	//  int x_init[8];
	//  MyComplex HH_H[64];
	//  MyComplex grad_preconditioner[64];
	//  MyComplex constellation_norm[16];
	//  MyComplex x_mmse[8];
	//  MyComplex covar[64];
	//  MyComplex pmat[64];
	 MyComplex x_survivor[8];
	//  MyComplex r[8];
	//  MyComplex pr_prev[8];
	//  MyComplex z_grad[8];
	//  MyComplex v[8];
	//  MyComplex z_prop[8];
	//  MyComplex x_prop[8];
	//  MyComplex r_prop[8];
	//  MyComplex temp_NtNt[64];
	//  MyComplex temp_NtNr[64];
	//  MyComplex temp_1;
	//  MyComplex _temp_1;
	//  MyComplex temp_Nt[8];
	//  MyComplex temp_Nr[8];
	//  MyComplex sigma2eye[64];
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

	// /float* MHGD_start;
	// /MHGD_start = (float*)malloc(time_size * sizeof(float));

    // transA = 1;
	#ifndef __SYNTHESIS__

	MHGD_start[(int)gradpre] = clock();

	#endif

	like_float numerator = 1.5;
    like_float denominator = hls::pow((ap_fixed<64,32>)2, (ap_fixed<64,32>)mu_local) - (ap_fixed<64,32>)1;
	like_float tttt = numerator/denominator;
	dqam = hls::sqrt((ap_fixed<64,32>)tttt);/*���ͷ���֮����С�����һ�룬�������㾭����һ��������Ľ��*/
	/*----------------------------------------------------------------*/
	c_eye_generate(sigma2eye, Nt_local, sigma2_local / (float)(dqam * dqam));
	/*�����ݶ��½�������grad_preconditioner*/
	c_matmultiple(H_local, transA, H_local, transB, Nr_local, Nt_local, Nr_local, Nt_local, HH_H);
	my_complex_add(Nt_local * Nt_local, HH_H, sigma2eye, grad_preconditioner);
	Inverse_LU(grad_preconditioner, Nt_local, Nt_local);

	// alpha = (like_float)1.0 / hls::pow((like_float)Nt / (like_float)8.0, (like_float)(1.0 / 3.0));
	like_float exponent = like_float(1) / like_float(3); // 避免浮点字面值隐式转换
	like_float exponent_1 = (like_float)Nt_local / like_float(8);
    alpha = like_float(1) / hls::pow<64,32>((ap_fixed<64,32>)exponent, (ap_fixed<64,32>)exponent_1);

	#ifndef __SYNTHESIS__

	MHGD_end[(int)gradpre] = clock();

	/*Э�������*/
	MHGD_start[covar_cal] = clock();
	#endif
	c_eye_generate(covar, Nt_local, 1.0f);
	#ifndef __SYNTHESIS__

	MHGD_end[covar_cal] = clock();

	MHGD_start[lr_cal] = clock();
	#endif
	/*For learning rate line search */
	if (!lr_approx_local)
	{
		c_matmultiple(H_local, transB, grad_preconditioner, transB, Nr_local, Nt_local, Nt_local, Nt_local, temp_NtNr);
		c_matmultiple(temp_NtNr, transB, H_local, transA, Nr_local, Nt_local, Nr_local, Nt_local, pmat);
	}
	else
	{
		 for (i = 0; i < Nr_local * Nr_local; i++)
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
	if (mmse_init_local)
	{
		// /*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
		c_eye_generate(sigma2eye, Nt_local, sigma2_local);
		my_complex_add(Nt_local * Nt_local, sigma2eye, HH_H, temp_NtNt);
		Inverse_LU(temp_NtNt, Nt_local, Nt_local);
		c_matmultiple(H_local, transA, y_local, transB, Nr_local, Nt_local, Nr_local, transA, temp_Nt);
		c_matmultiple(temp_NtNt, transB, temp_Nt, transB, Nt_local, Nt_local, Nt_local, transA, x_mmse);
		/*ӳ�䵽��һ������ͼ�� xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map(mu_local, Nt_local, dqam, x_mmse, x_hat_1);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		generateUniformRandoms_int(Nt_local, x_init, mu_local);
		for (i = 0; i < Nt_local; i++)
		{   
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			x_hat_1[i] = constellation_norm[x_init[i]];
		}
	}
	// /*����ʣ������r=y-Hx*/
	c_matmultiple(H_local, transB, x_hat_1, transB, Nr_local, Nt_local, Nt_local, transA, r);
	my_complex_sub(Nr_local, y_local, r, r);

	/*计算剩余向量的范数（就是模值）*/
	c_matmultiple(r, transA, r, transB, Nr_local, transA, Nr_local, transA, temp_1);
	r_norm = temp_1[0].real;
	my_complex_copy(Nt_local, x_hat_1, 1, x_survivor, 1);
	r_norm_survivor = r_norm;
	/*确定最优学习率*/
	if (!lr_approx_local)
	{
		c_matmultiple(pmat, transB, r, transB, Nr_local, Nr_local, Nr_local, transA, pr_prev);
		c_matmultiple(r, transA, pr_prev, transB, Nr_local, transA, Nr_local, transA, temp_1);
		c_matmultiple(pr_prev, transA, pr_prev, transB, Nr_local, transA, Nr_local, transA, _temp_1);
		lr = temp_1[0].real / _temp_1[0].real;
	}
	else
	{
		lr = 1;
	}

	step_size = max(dqam, like_float(hls::sqrt(r_norm / (like_float)Nr_local))) * alpha;
	#ifndef __SYNTHESIS__
	MHGD_end[mmse_and_r_init] = clock();
	#endif
	/*================core==============*/
	int offset = 0;  // ���ļ��Ŀ�ʼλ�ÿ�ʼ��ȡ
	for (k = 0; k < iter_local; k++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=iter_max min=iter_min
//	  	#pragma HLS PIPELINE II=1
		/*�����ݶ� z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
		#ifndef __SYNTHESIS__
		MHGD_start[zgrad] = clock();
		#endif
		c_matmultiple(H_local, transA, r, transB, Nr_local, Nt_local, Nr_local, transA, temp_Nt);
		c_matmultiple(grad_preconditioner, transB, temp_Nt, transB, Nt_local, Nt_local, Nt_local, transA, z_grad);
		my_complex_scal(Nt_local, lr, z_grad, 1); 
		my_complex_add(Nt_local, x_hat_1, z_grad, z_grad);
		#ifndef __SYNTHESIS__
		MHGD_end[zgrad] = clock();
		/*�����˹����Ŷ�*/
		MHGD_start[randomwalk] = clock();
		#endif
//		// read_gaussian_data("/home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt", v, Nt, offset);
		for(int i = 0; i < Nt_local; i++){
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			v[i].real = v_tb_local[i+offset].real;
			v[i].imag = v_tb_local[i+offset].imag;
		}
		offset = offset + Nt_local;
		my_complex_scal(Nt_local, step_size, v, 1);
		my_complex_add(Nt_local, z_grad, v, z_prop);
		#ifndef __SYNTHESIS__
		MHGD_end[randomwalk] = clock();
		/*���ݶ�ӳ�䵽QAM�������� x_prop = constellation_norm[np.argmin(abs(z_prop * ones - constellation_norm), axis=2)].reshape(-1, nt, 1) */
		MHGD_start[map_op] = clock();
		#endif
		map(mu_local, Nt_local, dqam, z_prop, x_prop);
		#ifndef __SYNTHESIS__
		MHGD_end[map_op] = clock();
		

		/*�����µĲв�� calculate residual norm of the proposal*/
		MHGD_start[new_r] = clock();
		#endif
		c_matmultiple(H_local, transB, x_prop, transB, Nr_local, Nt_local, Nt_local, transA, temp_Nr);
		my_complex_sub(Nr_local, y_local, temp_Nr, r_prop);
		c_matmultiple(r_prop, transA, r_prop, transB, Nr_local, transA, Nr_local, transA, temp_1);
		r_norm_prop = temp_1[0].real;

		/*update the survivor*/
		if (r_norm_survivor > r_norm_prop)
		{
			my_complex_copy(Nt_local, x_prop, 1, x_survivor, 1);
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
			my_complex_copy(Nt_local, x_prop, 1, x_hat_1, 1);
			my_complex_copy(Nt_local, r_prop, 1, r, 1);
			r_norm = r_norm_prop;
			/*update GD learning rate*/
			if (!lr_approx_local)
			{
				c_matmultiple(pmat, transB, r, transB, Nr_local, Nr_local, Nr_local, transA, pr_prev);
				c_matmultiple(r, transA, pr_prev, transB, Nr_local, transA, Nr_local, transA, temp_1);
				c_matmultiple(pr_prev, transA, pr_prev, transB, Nr_local, transA, Nr_local, transA, _temp_1);
				lr = temp_1[0].real / _temp_1[0].real;
			}
			/*update random walk size*/
			step_size = max(dqam, (like_float)hls::sqrt(r_norm / (like_float)Nr_local)) * alpha;
		}
		#ifndef __SYNTHESIS__
		MHGD_end[accept] = clock();

		for (i = zgrad; i <= accept; i++)
			MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
		#endif
	}
	// /*����ʵ��ʱ��*/
	#ifndef __SYNTHESIS__
	MHGD_start[x_copy] = clock();
	#endif
	// my_complex_copy(Nt, x_survivor, 1, x_hat_1, 1);
	
	//for (int i = 0; i < Nt; i++ ){
	//	#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
	//	x_hat_1[i].real = x_survivor[i].real;  // ???????
	//	x_hat_1[i].imag = x_survivor[i].imag;  // ?????��
	//}

	#ifndef __SYNTHESIS__
	MHGD_end[x_copy] = clock();
	for (i = 0; i <= mmse_and_r_init; i++)
		MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
	MHGD_dur[x_copy] += (MHGD_end[x_copy] - MHGD_start[x_copy]) / (float)CLOCKS_PER_SEC;
	#endif
	for(int i= 0;i<8;i++){
		*(x_hat_real+i) = x_survivor[i].real;
		*(x_hat_imag+i) = x_survivor[i].imag;
	}
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