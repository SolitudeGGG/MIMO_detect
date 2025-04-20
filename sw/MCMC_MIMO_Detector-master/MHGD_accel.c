#include "MHGD_accel.h"
#include "MyComplex.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


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
float* p_uni;
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
	p_uni = (float*)malloc(10 * sizeof(float));
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
	free(p_uni);

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
float MHGD_detect_accel(MyComplex* x_hat, int Nt, int Nr, int mu, MyComplex* H, MyComplex* y, float sigma2, int mmse_init, int lr_approx, int iter)
{
	/*��������*/
	float alpha;
	float dqam;
	float step_size;
	int i; int j; int k; int id;
	float log_pacc;
	float r_norm;/*the norm of residual vector*/
	float r_norm_survivor;
	float lr;/*learning rate*/
	float r_norm_prop;/*new norm*/
	float p_acc;

	int transA = 1;  // CblasConjTrans �ĵ�Чֵ����ʾ����ת��
	int transB = 0;  // CblasNoTrans �ĵ�Чֵ����ʾ��ת��

	MHGD_start[(int)gradpre] = clock();

	dqam = sqrt(1.5f / (pow(2.0f, (float)mu) - 1));/*���ͷ���֮����С�����һ�룬�������㾭����һ��������Ľ��*/


	/*----------------------------------------------------------------*/
	c_eye_generate(sigma2eye, Nt, sigma2 / (dqam * dqam));

	/*�����ݶ��½�������grad_preconditioner*/
	c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, HH_H);
	my_complex_add(Nt * Nt, HH_H, sigma2eye, grad_preconditioner);
	Inverse_LU(grad_preconditioner, Nt, Nt);

	alpha = 1.0f / pow((float)Nt / 8.0f, 1.0f / 3.0f);

	MHGD_end[(int)gradpre] = clock();

	/*Э�������*/
	MHGD_start[covar_cal] = clock();
	c_eye_generate(covar, Nt, 1.0f);
	MHGD_end[covar_cal] = clock();

	MHGD_start[lr_cal] = clock();
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
			pmat[i].real = 0; pmat[i].imag = 0;
		}
	}
	MHGD_end[lr_cal] = clock();

	/*initialize the estimate with size (np, nt, 1), np is for parallel samplers*/
	MHGD_start[mmse_and_r_init] = clock();
	if (mmse_init)
	{
		/*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
		c_eye_generate(sigma2eye, Nt, sigma2);
		my_complex_add(Nt * Nt, sigma2eye, HH_H, temp_NtNt);
		Inverse_LU(temp_NtNt, Nt, Nt);
		c_matmultiple(H, transA, y, transB, Nr, Nt, Nr, 1, temp_Nt);
		c_matmultiple(temp_NtNt, transB, temp_Nt, transB, Nt, Nt, Nt, 1, x_mmse);
		/*ӳ�䵽��һ������ͼ�� xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map(mu, Nt, dqam, x_mmse, x_hat);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		generateUniformRandoms_int(Nt, x_init, mu);
		for (i = 0; i < Nt; i++)
			x_hat[i] = constellation_norm[x_init[i]];
	}
	/*����ʣ������r=y-Hx*/
	c_matmultiple(H, transB, x_hat, transB, Nr, Nt, Nt, 1, r);
	my_complex_sub(Nr, y, r, r);

	/*����ʣ�������ķ���������ģֵ��*/
	c_matmultiple(r, transA, r, transB, Nr, 1, Nr, 1, temp_1);
	r_norm = temp_1[0].real;
	my_complex_copy(Nt, x_hat, 1, x_survivor, 1);
	r_norm_survivor = r_norm;
	/*ȷ������ѧϰ��*/
	if (!lr_approx)
	{
		c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, 1, pr_prev);
		c_matmultiple(r, transA, pr_prev, transB, Nr, 1, Nr, 1, temp_1);
		c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, 1, Nr, 1, _temp_1);
		lr = temp_1[0].real / _temp_1[0].real;
	}
	else
	{
		lr = 1;
	}

	step_size = max(dqam, sqrt(r_norm / (float)Nr)) * alpha;
	MHGD_end[mmse_and_r_init] = clock();
	/*================core==============*/
	int offset = 0;  // ���ļ��Ŀ�ʼλ�ÿ�ʼ��ȡ
	for (k = 0; k < iter; k++)
	{
		/*�����ݶ� z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
		MHGD_start[zgrad] = clock();
		c_matmultiple(H, transA, r, transB, Nr, Nt, Nr, 1, temp_Nt);
		c_matmultiple(grad_preconditioner, transB, temp_Nt, transB, Nt, Nt, Nt, 1, z_grad);
		my_complex_scal(Nt, lr, z_grad, 1);
		my_complex_add(Nt, x_hat, z_grad, z_grad);
		MHGD_end[zgrad] = clock();
		/*�����˹����Ŷ�*/
		MHGD_start[randomwalk] = clock();
		read_gaussian_data("gaussian_random_values.txt", v, Nt, offset);
		offset = offset + Nt;
		my_complex_scal(Nt, step_size, v, 1);
		my_complex_add(Nt, z_grad, v, z_prop);
		MHGD_end[randomwalk] = clock();
		/*���ݶ�ӳ�䵽QAM�������� x_prop = constellation_norm[np.argmin(abs(z_prop * ones - constellation_norm), axis=2)].reshape(-1, nt, 1) */
		MHGD_start[map_op] = clock();
		map(mu, Nt, dqam, z_prop, x_prop);
		MHGD_end[map_op] = clock();

		/*�����µĲв�� calculate residual norm of the proposal*/

		MHGD_start[new_r] = clock();
		c_matmultiple(H, transB, x_prop, transB, Nr, Nt, Nt, 1, temp_Nr);
		my_complex_sub(Nr, y, temp_Nr, r_prop);
		c_matmultiple(r_prop, transA, r_prop, transB, Nr, 1, Nr, 1, temp_1);
		r_norm_prop = temp_1[0].real;

		/*update the survivor*/
		if (r_norm_survivor > r_norm_prop)
		{
			my_complex_copy(Nt, x_prop, 1, x_survivor, 1);
			r_norm_survivor = r_norm_prop;
		}
		MHGD_end[new_r] = clock();

		/*acceptance test*/
		MHGD_start[accept] = clock();
		log_pacc = min(0, -(r_norm_prop - r_norm));
		p_acc = exp(log_pacc);
		generateUniformRandoms_float(10, p_uni);
		if (p_acc > p_uni[5])/*������������ʱ��*/
		{
			my_complex_copy(Nt, x_prop, 1, x_hat, 1);
			my_complex_copy(Nt, r_prop, 1, r, 1);
			r_norm = r_norm_prop;
			/*update GD learning rate*/
			if (!lr_approx)
			{
				c_matmultiple(pmat, transB, r, transB, Nr, Nr, Nr, 1, pr_prev);
				c_matmultiple(r, transA, pr_prev, transB, Nr, 1, Nr, 1, temp_1);
				c_matmultiple(pr_prev, transA, pr_prev, transB, Nr, 1, Nr, 1, _temp_1);
				lr = temp_1[0].real / _temp_1[0].real;
			}
			/*update random walk size*/
			step_size = max(dqam, sqrt(r_norm / (float)Nr)) * alpha;
		}
		MHGD_end[accept] = clock();

		for (i = zgrad; i <= accept; i++)
			MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
	}
	/*����ʵ��ʱ��*/
	MHGD_start[x_copy] = clock();
	my_complex_copy(Nt, x_survivor, 1, x_hat, 1);
	MHGD_end[x_copy] = clock();
	for (i = 0; i <= mmse_and_r_init; i++)
		MHGD_dur[i] += (MHGD_end[i] - MHGD_start[i]) / (float)CLOCKS_PER_SEC;
	MHGD_dur[x_copy] += (MHGD_end[x_copy] - MHGD_start[x_copy]) / (float)CLOCKS_PER_SEC;

	return 0;
}


void Hx_accel(int Nt, int Nr, int mu, float dqam, MyComplex* x, MyComplex* x_result)
{
	my_complex_scal(Nr, 0, x_result, 1);
	int i, j;
	int Hx_table_id;
#ifdef print_debug
	printf("x:\n");
	c_print_vector(x, Nt);
#endif
	/*������Ҫ���н������temp_real��temp_imag��*/
	for (i = 0; i < Nt; i++)
	{
		Hx_table_id = (x[i].real / dqam + mu - 1) / 2;
#ifdef print_debug
		printf("Hx_table_id:%d  ", Hx_table_id);
#endif
		for (j = 0; j < Nr; j++)
		{
			temp_real[i][j] = table_real[Hx_table_id][Nt * j + i];
		}
		Hx_table_id = (x[i].imag / dqam + mu - 1) / 2;
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
	my_complex_scal(Nr, dqam, x_result, 1);
}