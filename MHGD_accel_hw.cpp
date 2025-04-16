#include "MyComplex_1.h"
#include "hls_math.h"
#include "MHGD_accel_hw.h"
#include <stdio.h>

MyComplex QPSK_Constellation_hw[4] = {{-1,-1},{-1,1},{1,-1},{1,1}};
MyComplex _16QAM_Constellation_hw[16] = {{-3,-3},{-3,-1},{-3,3},{-3,1},{-1,-3},{-1,-1},{-1,3},{-1,1},{3,-3},{3,-1},{3,3},{3,1},{1,-3},{1,-1},{1,3},{1,1}};
MyComplex _64QAM_Constellation_hw[64] = {{-7.0, 7.0},  {-5.0, 7.0},  {-3.0, 7.0}, {-1.0, 7.0},
                                         {1.0, 7.0},   {3.0, 7.0},   {5.0, 7.0},  {7.0, 7.0},
                                         {-7.0, 5.0},  {-5.0, 5.0},  {-3.0, 5.0}, {-1.0, 5.0},
                                         {1.0, 5.0},   {3.0, 5.0},   {5.0, 5.0},  {7.0, 5.0},
                                         {-7.0, 3.0},  {-5.0, 3.0},  {-3.0, 3.0}, {-1.0, 3.0},
                                         {1.0, 3.0},   {3.0, 3.0},   {5.0, 3.0},  {7.0, 3.0},
                                         {-7.0, 1.0},  {-5.0, 1.0},  {-3.0, 1.0}, {-1.0, 1.0},
                                         {1.0, 1.0},   {3.0, 1.0},   {5.0, 1.0},  {7.0, 1.0},
                                         {-7.0, -1.0}, {-5.0, -1.0}, {-3.0, -1.0},{-1.0, -1.0},
                                         {1.0, -1.0},  {3.0, -1.0},  {5.0, -1.0}, {7.0, -1.0},
                                         {-7.0, -3.0}, {-5.0, -3.0}, {-3.0, -3.0},{-1.0, -3.0},
                                         {1.0, -3.0},  {3.0, -3.0},  {5.0, -3.0}, {7.0, -3.0},
                                         {-7.0, -5.0}, {-5.0, -5.0}, {-3.0, -5.0},{-1.0, -5.0},
                                         {1.0, -5.0},  {3.0, -5.0},  {5.0, -5.0}, {7.0, -5.0},
                                         {-7.0, -7.0}, {-5.0, -7.0}, {-3.0, -7.0},{-1.0, -7.0},
                                         {1.0, -7.0},  {3.0, -7.0},  {5.0, -7.0}, {7.0, -7.0}};

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/***************************************部分供tb调用的C代码***********************************/
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
// 读取已经生成的高斯随机变量数据
void read_gaussian_data_hw(const char* filename, MyComplex* array, int n, int offset) {
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening file!\n");
		return;  // 错误打开文件
	}

	int i;
	float temp_real, temp_imag;

	// 跳过文件中的前offset个数据，继续读取
	for (i = 0; i < offset; i++) {
		// 读取每一行的格式 "Real: <value>, Imaginary: <value>"
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data at offset %d!\n", i);
			fclose(f);
			return;  // 错误读取数据
		}
	}

	// 读取n个数据点
	for (i = 0; i < n; i++) {
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data!\n");
			fclose(f);
			return;  // 返回已成功读取的数据数量
		}

		// 将读取的实部和虚部数据存入数组
		array[i].real = temp_real / (float)sqrt(2.0f);
		array[i].imag = temp_imag / (float)sqrt(2.0f);
	}

	fclose(f);
}

void QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int mu, int* bits_demod)
{
	switch (mu)
	{
	case 2:
		QPSK_Demodulation_hw(x_hat, Nt, bits_demod);
		break;
	case 4:
		_16QAM_Demodulation_hw(x_hat, Nt, bits_demod);
		break;
	case 6:
		_64QAM_Demodulation_hw(x_hat, Nt, bits_demod);
		break;
	default:
		break;
	}
}

/*QPSK解调*/
void QPSK_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex* temp = (MyComplex*)malloc(4 * sizeof(MyComplex));
	float distance[4] = { 0 };
	/*for one in x_hat*/
	for (i = 0; i < Nt; i++)
	{
		for (j = 0; j < 4; j++)
		{
			temp[j].real = x_hat[i].real * (like_float)hls::sqrt(2.0) - QPSK_Constellation_hw[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)hls::sqrt(2.0) - QPSK_Constellation_hw[j].imag;
		}
		for (j = 0; j < 4; j++)
		{
			distance[j] = std::sqrt((float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag));
		}
		best_id = argmin_hw(distance, 4);
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
void _16QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex* temp = (MyComplex*)malloc(16 * sizeof(MyComplex));
	float distance[16] = { 0 };
	/*for one in x_hat*/
	for (i = 0; i < Nt; i++)
	{
		for (j = 0; j < 16; j++)
		{
			temp[j].real = x_hat[i].real * (like_float)hls::sqrt(10.0) - _16QAM_Constellation_hw[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)hls::sqrt(10.0) - _16QAM_Constellation_hw[j].imag;
		}
		/*find the closest one*/
		for (j = 0; j < 16; j++)
		{
			distance[j] = std::sqrt((float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag));
		}
		best_id = argmin_hw(distance, 16);
		bits_demod[4 * i] = (best_id >> 3) & 0x01;
		bits_demod[4 * i + 1] = (best_id >> 2) & 0x01;
		bits_demod[4 * i + 2] = (best_id >> 1) & 0x01;
		bits_demod[4 * i + 3] = best_id & 0x01;
	}

	free(temp);
}

/*64QAM解调*/
void _64QAM_Demodulation_hw(MyComplex* x_hat, int Nt, int* bits_demod)
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
			temp[j].real = x_hat[i].real * (like_float)sqrt(42.0) - _64QAM_Constellation_hw[j].real;
			temp[j].imag = x_hat[i].imag * (like_float)sqrt(42.0) - _64QAM_Constellation_hw[j].imag;
			distance[j] = (float)(temp[j].real * temp[j].real + temp[j].imag * temp[j].imag);
		}
		best_id = argmin_hw(distance, 64);
		switch ((int)_64QAM_Constellation_hw[best_id].real)
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
		switch ((int)_64QAM_Constellation_hw[best_id].imag)
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
int argmin_hw(float* array, int n)
{
	int i = 0; int best_id = 0;
	like_float minval = (like_float)array[0];
	for (i = 1; i < n; i++)
	{
		if (minval > (like_float)array[i])
		{
			minval = (like_float)array[i];
			best_id = i;
		}
	}
	return best_id;
}
/*计算误bit数*/
int unequal_times_hw(int* array1, int* array2, int n)
{
	int unequal = 0; int i = 0;
	for (i = 0; i < n; i++)
	{
		if (array1[i] != array2[i])
			unequal++;
	}
	return unequal;
}
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/***************************************功能函数***********************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
// 线性同余生成器
unsigned int lcg_rand_hw() {
    unsigned int lcg_seed = 1;
	lcg_seed = (LCG_A * lcg_seed + LCG_C) % LIMIT_MAX;;
	
    return lcg_seed;
}
like_float lcg_rand_1_hw() {
    unsigned int lcg_seed = 1;
	lcg_seed = (LCG_A * lcg_seed + LCG_C) % LIMIT_MAX;
    return (like_float)(lcg_seed / LIMIT_MAX_F);
}
// 复数取共轭
MyComplex complex_conjugate_hw(MyComplex a) {
	MyComplex result;
	result.real = a.real;   // 实部不变
	result.imag = -a.imag;  // 虚部取负
	return result;
}
void get_dqam_hw(like_float dqam){
    like_float numerator = 1.5;
    like_float denominator = hls::pow((ap_fixed<64,32>)2, (ap_fixed<64,32>)mu_1) - (ap_fixed<64,32>)1;
	like_float tttt = numerator/denominator;
	dqam = hls::sqrt((ap_fixed<64,32>)tttt);
}
void c_eye_generate_hw(MyComplex* Mat, float val)
{
	like_float val_1 = val;
    int i;
	for (i = 0; i < Ntr_2; i++)
	{
		Mat[i].real = (like_float)0.0;
		Mat[i].imag = (like_float)0.0;
	}
	for (i = 0; i < Ntr_1; i++)
	{
		Mat[i * Ntr_1 + i].real = val;
	}
}
void data_local(MyComplex* H, MyComplex* y, MyComplex* v_tb, Myreal* H_real, Myimage* H_imag, Myreal* y_real, Myimage* y_imag, Myreal* v_tb_real, Myimage* v_tb_imag)
{
    for(int i=0; i<Ntr_1; i++){
		y[i].real = *(y_real+i);
		y[i].imag = *(y_imag+i);
	}
	for (int i = 0; i < Ntr_1*Ntr_1; i++)
	{
		H[i].real = *(H_real+i);
		H[i].imag = *(H_imag+i);
	}
	for (int i = 0; i < Ntr_1*iter_1; i++)
	{
		v_tb[i].real = *(v_tb_real+i);
		v_tb[i].imag = *(v_tb_imag+i);
	}
}
void c_matmultiple_hw(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res)
{
	int m, n, k;
	// 根据转置的情况计算 m, n, k
	if (transA == 0) {  // CblasNoTrans
		m = ma;
		k = na;
	}
	else {
		k = ma;
		m = na;
	}
	if (transB == 0) {  // CblasNoTrans
		n = nb;
	}
	else {
		n = mb;
	}
	// 逐个计算矩阵元素
	for (int i = 0; i < m; i++) {
        #pragma HLS LOOP_TRIPCOUNT max=Ntr_1
		for (int j = 0; j < n; j++) {
            #pragma HLS LOOP_TRIPCOUNT max=Ntr_1
			MyComplex sum =  { (like_float)0.0, (like_float)0.0 }; // 初始化为复数零
			MyComplex temp = { (like_float)0.0, (like_float)0.0 };
			for (int l = 0; l < k; l++) {
				MyComplex a_element, b_element;
				// 获取矩阵A的元素
				if (transA == 1) {  // 共轭转置
					a_element = matA[l * na + i];  // 获取矩阵A[l, i]
					a_element = complex_conjugate_hw(a_element);  // 对元素取共轭
				}
				else if (transA == 0) {  // 普通转置
					a_element = matA[i * na + l];  // 获取矩阵A[i, l]
				}
				// 获取矩阵B的元素
				if (transB == 0) {  // No Transpose
					// 问题点1-1
					b_element = matB[l * nb + j];  // 获取矩阵B[l, j]
				}
				else {  // 转置
					b_element = matB[j * nb + l];  // 获取矩阵B[j, l]
					b_element = complex_conjugate_hw(b_element);  // 对元素取共轭
				}
				// 复数乘法并累加
				temp = complex_multiply_hw(a_element, b_element);
				sum = complex_add_hw(sum, temp);
			}
			// 将结果赋值到矩阵C（res）
			res[i * nb + j].real = sum.real;
            res[i * nb + j].imag = sum.imag;
		}
	}
}
void my_complex_add_hw(const MyComplex a[], const MyComplex b[], MyComplex r[])
{
	for (int i = 0; i < Ntr_2; i++) {
		r[i].real = a[i].real + b[i].real; 
		r[i].imag = a[i].imag + b[i].imag;  
	}
}
void my_complex_add_hw_1(const MyComplex a[], const MyComplex b[], MyComplex r[])
{
	for (int i = 0; i < Ntr_1; i++) {
		r[i].real = a[i].real + b[i].real; 
		r[i].imag = a[i].imag + b[i].imag;  
	}
}
// 初始化复数矩阵
void initMatrix_hw(MyComplex* A)
{
    int i, j;
	for (i = 0; i < Ntr_1; i++)
	{
		for (j = 0; j < Ntr_1; j++)
		{
			A[i * Ntr_1 + j].real = (like_float)0.0;
			A[i * Ntr_1 + j].imag = (like_float)0.0;
		}
	}
}
// 复数除法
MyComplex complex_divide_hw(MyComplex a, MyComplex b)
{
    MyComplex result;
    like_float denominator = b.real * b.real + b.imag * b.imag;
    result.real = (a.real * b.real + a.imag * b.imag) / denominator;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
    return result;
}
// 复数乘法
MyComplex complex_multiply_hw(MyComplex a, MyComplex b)
{
    MyComplex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}
// 复数加法
MyComplex complex_add_hw(MyComplex a, MyComplex b)
{
    MyComplex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}
// 复数减法
MyComplex complex_subtract_hw(MyComplex a, MyComplex b)
{
    MyComplex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}
// 复数数组减法
void my_complex_sub_hw(const MyComplex a[], const MyComplex b[], MyComplex r[])
{
	for (int i = 0; i < Ntr_1; i++) {
		r[i].real = a[i].real - b[i].real;  
		r[i].imag = a[i].imag - b[i].imag; 
	}
}
// 矩阵乘法（用于求逆函数中调用）
void MulMatrix_hw(const MyComplex* A, const MyComplex* B, MyComplex* C)
{
	int i, j, k;
	MyComplex sum, tmp;
	for (i = 0; i < Ntr_1; i++)
	{
        #pragma HLS pipeline
		for (j = 0; j < Ntr_1; j++)
		{
            #pragma HLS pipeline
			sum.real = 0.0;
			sum.imag = 0.0;
			for (k = 0; k < Ntr_1; k++)
			{
                #pragma HLS pipeline
				sum = complex_add_hw(sum, complex_multiply_hw(A[i * Ntr_1 + k], B[k * Ntr_1 + j]));
			}
			C[i * Ntr_1 + j].real = sum.real;
			C[i * Ntr_1 + j].imag = sum.imag;
		}
	}
}
// 复数矩阵求逆
void Inverse_LU_hw(MyComplex* A)
{
    MyComplex _L[Ntr_2];//row * col
	MyComplex _U[Ntr_2];//row * col
	MyComplex _L_Inverse[Ntr_2];//row * col
	MyComplex _U_Inverse[Ntr_2];//row * col

	MyComplex* L = &_L[0];
	MyComplex* U = &_U[0];
	MyComplex* L_Inverse = &_L_Inverse[0];
	MyComplex* U_Inverse = &_U_Inverse[0];
    initMatrix_hw(L);
	initMatrix_hw(U);
	initMatrix_hw(L_Inverse);
	initMatrix_hw(U_Inverse);
    int i, j, k, t;
	MyComplex s;
	for (i = 0; i < Ntr_1; i++)//??????
	{
		L[i * Ntr_1 + i].real = (like_float)1.0;
		L[i * Ntr_1 + i].imag = (like_float)0.0;
	}
	for (j = 0; j < Ntr_1; j++)
	{
		U[j] = A[j];
	}
	for (i = 1; i < Ntr_1; i++)
	{
		L[i * Ntr_1] = complex_divide_hw(A[i * Ntr_1], U[0]);

	}
	for (k = 1; k < Ntr_1; k++)
	{
		for (j = k; j < Ntr_1; j++)
		{
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (t = 0; t < k; t++) {
				s = complex_add_hw(s, complex_multiply_hw(L[k * Ntr_1 + t], U[t * Ntr_1 + j]));
			}
			U[k * Ntr_1 + j] = complex_subtract_hw(A[k * Ntr_1 + j], s);
		}
		for (i = k; i < Ntr_1; i++)
		{
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (t = 0; t < k; t++)
			{
				s = complex_add_hw(s, complex_multiply_hw(L[i * Ntr_1 + t], U[t * Ntr_1 + k]));
			}
			L[i * Ntr_1 + k] = complex_divide_hw(complex_subtract_hw(A[i * Ntr_1 + k], s), U[k * Ntr_1 + k]);
		}
	}

	for (i = 0; i < Ntr_1; i++)
	{
		L_Inverse[i * Ntr_1 + i].imag = (like_float)0.0;
		L_Inverse[i * Ntr_1 + i].real = (like_float)1.0;
	}
	for (j = 0; j < Ntr_1; j++)
	{
		for (i = j + 1; i < Ntr_1; i++)
		{
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (k = j; k < i; k++) {
				s = complex_add_hw(s, complex_multiply_hw(L[i * Ntr_1 + k], L_Inverse[k * Ntr_1 + j]));
			}
			L_Inverse[i * Ntr_1 + j].real = (-1) * complex_multiply_hw(L_Inverse[j * Ntr_1 + j], s).real;
			L_Inverse[i * Ntr_1 + j].imag = (-1) * complex_multiply_hw(L_Inverse[j * Ntr_1 + j], s).imag;
		}
	}
	MyComplex di;
	di.real = (like_float)1.0;
	di.imag = (like_float)0.0;
	for (i = 0; i < Ntr_1; i++) //按列序，列内按照从下到上，计算u的逆矩阵
	{
		U_Inverse[i * Ntr_1 + i] = complex_divide_hw(di, U[i * Ntr_1 + i]);

	}
	for (j = 0; j < Ntr_1; j++)
	{
		for (i = j - 1; i >= 0; i--)
		{
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (k = i + 1; k <= j; k++)
			{
				s = complex_add_hw(s, complex_multiply_hw(U[i * Ntr_1 + k], U_Inverse[k * Ntr_1 + j]));
			}
			s.imag = -s.imag;
			s.real = -s.real;
			U_Inverse[i * Ntr_1 + j] = complex_divide_hw(s, U[i * Ntr_1 + i]);
		}
	}
	MulMatrix_hw(U_Inverse, L_Inverse, A);
}
// 缩放复数数组
void my_complex_scal_hw(const like_float alpha, MyComplex* X, const int incX)
{
	for (int i = 0; i < Ntr_1; i++) {
		X[i * incX].real *= alpha;  // 实部乘以 alpha
		X[i * incX].imag *= alpha;  // 虚部乘以 alpha
	}
}
void my_complex_scal_hw_1(const like_float alpha, MyComplex* X, const int incX)
{
	for (int i = 0; i < 16; i++) {
		X[i * incX].real *= alpha;  // 实部乘以 alpha
		X[i * incX].imag *= alpha;  // 虚部乘以 alpha
	}
}
// 星座点映射
void map_hw(like_float dqam, MyComplex* x, MyComplex* x_hat)
{
	int i;
	for (i = 0; i < Ntr_1; i++)
	{
		x_hat[i].real = ((like_float)2.0 * (like_float)hls::floor(x[i].real / dqam / (like_float)2.0) + (like_float)1.0);
		x_hat[i].imag = ((like_float)2.0 * (like_float)hls::floor(x[i].imag / dqam / (like_float)2.0) + (like_float)1.0);
		switch (mu_1)
		{
		case 2:
			x_hat[i].real = (x_hat[i].real > (like_float) 1.0  ? (like_float)1.0 : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < (like_float)-1.0 ? (like_float)-1.0 : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > (like_float) 1.0  ? (like_float)1.0 : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < (like_float)-1.0 ? (like_float)-1.0 : x_hat[i].imag);
			break;
		case 4:
			x_hat[i].real = (x_hat[i].real > (like_float) 3.0  ? (like_float)3.0 : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < (like_float)-3.0 ? (like_float)-3.0 : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > (like_float) 3.0  ? (like_float)3.0 : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < (like_float)-3.0 ? (like_float)-3.0 : x_hat[i].imag);
			break;
		case 6:
			x_hat[i].real = (x_hat[i].real > (like_float) 7.0  ? (like_float)7.0 : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < (like_float)-7.0 ? (like_float)-7.0 : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > (like_float) 7.0  ? (like_float)7.0 : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < (like_float)-7.0 ? (like_float)-7.0 : x_hat[i].imag);
			break;
		default:
			break;
		}
	}
	my_complex_scal_hw(dqam, x_hat, 1);
}
// 生成均匀分布随机整数
void generateUniformRandoms_int_hw(int* x_init)
{
	// 计算 2^mu
	double upper_bound = pow(2, mu_1);  // 计算 2^mu
	// 生成 Nt 个随机数
	for (int i = 0; i < Ntr_1; ++i) {
		// 使用 rand() 生成一个 [0, RAND_MAX] 范围内的整数，然后缩放到 [0, upper_bound]
		double rand_val = ((double)lcg_rand_hw() / LIMIT_MAX) * upper_bound;
		// 将生成的随机数存储到 x_init 数组中
		x_init[i] = (int)rand_val;
	}
}
// 复制复数数组
void my_complex_copy_hw(const MyComplex* X, const int incX, MyComplex* Y, const int incY)
{
	for (int i = 0, j = 0; i < Ntr_1; i++, j++) {
		Y[j * incY].real = X[i * incX].real;  // 复制实部
		Y[j * incY].imag = X[i * incX].imag;  // 复制虚部
	}
}
// 复制复数数组
void my_complex_copy_hw_1(const MyComplex* X, const int incX, MyComplex* Y, const int incY)
{
	for (int i = 0, j = 0; i < 16; i++, j++) {
		Y[j * incY].real = X[i * incX].real;  // 复制实部
		Y[j * incY].imag = X[i * incX].imag;  // 复制虚部
	}
}
// 生成均匀分布随机变量
void generateUniformRandoms_float_hw(like_float* p_uni)
{//一个随机数生成的函数 其中lcg_rand用于生成随机数，/LIMIT_MAX用于限幅
	for (int i = 0; i < 10; i++) {
		p_uni[i] = lcg_rand_1_hw();
	}
}
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/*********************************顶层函数分布计算函数******** **********************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/*计算grad_preconditioner(梯度更新的预条件矩阵)*/
void grad_preconditioner_updater_hw(MyComplex* H, int transA, int transB, MyComplex* HH_H, MyComplex* sigma2eye, MyComplex* grad_preconditioner, int Nt, int Nr)
{
    c_matmultiple_hw(H, transA, H, transB, Nr, Nt, Nr, Nt, HH_H);
    my_complex_add_hw(HH_H, sigma2eye, grad_preconditioner);
    Inverse_LU_hw(grad_preconditioner);
}
/*For learning rate line search*/
void learning_rate_line_search_hw(int lr_approx, MyComplex* H, int transA, int transB, MyComplex* grad_preconditioner, int Nr, int Nt, MyComplex* temp_NtNr, MyComplex* pmat)
{
    if (!lr_approx)
	{
		c_matmultiple_hw(H, transB, grad_preconditioner, transB, Nr, Nt, Nt, Nt, temp_NtNr);
		c_matmultiple_hw(temp_NtNr, transB, H, transA, Nr, Nt, Nr, Nt, pmat);
	}
	else
	{
		 for (int i = 0; i < Ntr_2; i++)
		{
			pmat[i].real = (like_float)0; pmat[i].imag = (like_float)0;
		}
	}
}
void x_initialize_hw(int mmse_init, MyComplex* sigma2eye, int Nt, int Nr, float sigma2, MyComplex* HH_H, MyComplex* temp_NtNt, MyComplex* H, int transA, int transB, MyComplex* y, MyComplex* temp_Nt, MyComplex* x_mmse, like_float dqam, MyComplex* x_hat, int* x_init, MyComplex* constellation_norm)
{
    if (mmse_init)
	{
		/*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
		c_eye_generate_hw(sigma2eye, sigma2);
		my_complex_add_hw(sigma2eye, HH_H, temp_NtNt);
		Inverse_LU_hw(temp_NtNt);
		c_matmultiple_hw(H, transA, y, transB, Nr, Nt, Nr, transA, temp_Nt);
		c_matmultiple_hw(temp_NtNt, transB, temp_Nt, transB, Nt, Nt, Nt, transA, x_mmse);
		/*映射到归一化星座图中 xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map_hw(dqam, x_mmse, x_hat);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		generateUniformRandoms_int_hw(x_init);
		for (int i = 0; i < Ntr_1; i++)
			x_hat[i] = constellation_norm[x_init[i]];
	}
}
/*计算剩余向量r=y-Hx*/
void r_hw(MyComplex* H, int transB, int transA, MyComplex* x_hat, int Nr, int Nt, MyComplex* r, MyComplex* y)
{
    c_matmultiple_hw(H, transB, x_hat, transB, Nr, Nt, Nt, transA, r);
	my_complex_sub_hw(y, r, r);
}
/*计算剩余向量的范数（就是模值）*/
void r_cal_hw(MyComplex* r, int transA, int transB, int Nr, int Nt, MyComplex* temp_1, MyComplex* x_hat, MyComplex* x_survivor, like_float r_norm, like_float r_norm_survivor)
{
    c_matmultiple_hw(r, transA, r, transB, Nr, transA, Nr, transA, temp_1);
	r_norm = temp_1->real;
	my_complex_copy_hw(x_hat, 1, x_survivor, 1);
	r_norm_survivor = r_norm;
}
/*确定最优学习率*/
void lr_hw(int lr_approx, MyComplex* pmat, int transB, int transA, MyComplex* r, int Nr, int Nt, MyComplex* pr_prev, MyComplex* temp_1, MyComplex* _temp_1, like_float lr)
{
    if (!lr_approx)
	{
		c_matmultiple_hw(pmat, transB, r, transB, Nr, Nr, Nr, transA, pr_prev);
		c_matmultiple_hw(r, transA, pr_prev, transB, Nr, transA, Nr, transA, temp_1);
		c_matmultiple_hw(pr_prev, transA, pr_prev, transB, Nr, transA, Nr, transA, _temp_1);
		lr = temp_1->real / _temp_1->real;
	}
	else
	{
		lr = 1;
	}
}
/*更新梯度 z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
void z_grad_hw(MyComplex* H, int transA, int transB, MyComplex* temp_Nt, MyComplex* grad_preconditioner, MyComplex* z_grad, like_float lr, MyComplex* x_hat, MyComplex* r)
{
    c_matmultiple_hw(H, transA, r, transB, Ntr_1, Ntr_1, Ntr_1, transA, temp_Nt);
	c_matmultiple_hw(grad_preconditioner, transB, temp_Nt, transB, Ntr_1, Ntr_1, Ntr_1, transA, z_grad);
	my_complex_scal_hw(lr, z_grad, 1); 
	my_complex_add_hw_1(x_hat, z_grad, z_grad);
}
/*加入高斯随机扰动*/
void gauss_add_hw(MyComplex* v, MyComplex* v_tb, int offset, like_float step_size, MyComplex* z_grad, MyComplex* z_prop)
{
    for(int i = 0; i < Ntr_1; i++){
        v[i].real = v_tb[i+offset].real;
        v[i].imag = v_tb[i+offset].imag;
    }
    offset = offset + Ntr_1;
    my_complex_scal_hw(step_size, v, 1);
    my_complex_add_hw_1(z_grad, v, z_prop);
}
/*计算新的残差范数 calculate residual norm of the proposal*/
void r_newnorm_hw(MyComplex* H, int transB, MyComplex* x_prop, int transA, MyComplex* temp_Nr, MyComplex* y, MyComplex* r_prop, MyComplex* temp_1, like_float r_norm_prop)
{
    c_matmultiple_hw(H, transB, x_prop, transB, Ntr_1, Ntr_1, Ntr_1, transA, temp_Nr);
	my_complex_sub_hw(y, temp_Nr, r_prop);
	c_matmultiple_hw(r_prop, transA, r_prop, transB, Ntr_1, transA, Ntr_1, transA, temp_1);
	r_norm_prop = temp_1->real;
}
/*update the survivor*/
void survivor_hw(like_float r_norm_survivor, like_float r_norm_prop, MyComplex* x_prop, MyComplex* x_survivor)
{
    if (r_norm_survivor > r_norm_prop)
	{
		my_complex_copy_hw(x_prop, 1, x_survivor, 1);
		r_norm_survivor = r_norm_prop;
	}
}
/*acceptance test＆update GD learning rate＆update random walk size*/
void acceptance_hw(int transB, int transA, like_float r_norm_prop, like_float r_norm, like_float log_pacc, like_float p_acc, like_float* p_uni, MyComplex* x_prop , MyComplex* x_hat, MyComplex* r_prop, MyComplex* r, MyComplex* pmat, MyComplex* pr_prev, MyComplex* temp_1, MyComplex* _temp_1, like_float lr, like_float step_size, like_float dqam, like_float alpha)
{
    like_float temp_3 = -(r_norm_prop - r_norm);
	log_pacc = min(like_float(0), temp_3);
	p_acc = hls::exp<64,32>((ap_fixed<64,32>)log_pacc);
	generateUniformRandoms_float_hw(p_uni);
	if (p_acc > p_uni[5])/*概率满足条件时候*/
	{
		my_complex_copy_hw(x_prop, 1, x_hat, 1);
		my_complex_copy_hw(r_prop, 1, r, 1);
		r_norm = r_norm_prop;
		/*update GD learning rate*/
		if (!lr_approx_1)
		{
			c_matmultiple_hw(pmat, transB, r, transB, Ntr_1, Ntr_1, Ntr_1, transA, pr_prev);
			c_matmultiple_hw(r, transA, pr_prev, transB, Ntr_1, transA, Ntr_1, transA, temp_1);
			c_matmultiple_hw(pr_prev, transA, pr_prev, transB, Ntr_1, transA, Ntr_1, transA, _temp_1);
			lr = temp_1->real / _temp_1->real;
		}
		/*update random walk size*/
		step_size = max(dqam, (like_float)hls::sqrt(r_norm / (like_float)Ntr_1)) * alpha;
	}
}
/*输出*/
void out_hw(const MyComplex* X, const int incX, Myreal* Y_real, Myimage* Y_imag, int incY)
{
	for (int i = 0, j = 0; i < Ntr_1; i++, j++) {
		Y_real[j * incY] = X[i * incX].real;  // 复制实部
		Y_imag[j * incY] = X[i * incX].imag;  // 复制虚部
	}
}

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/***************************************顶层函数***********************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

void MHGD_detect_accel_hw(
    Myreal* x_hat_real, Myimage* x_hat_imag, 
    Myreal* H_real, Myimage* H_imag, 
    Myreal* y_real, Myimage* y_imag, 
    float sigma2,
    Myreal* v_tb_real, Myimage* v_tb_imag
){
    /*hls interface＆接口数据本地化*/
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
    #pragma HLS INTERFACE mode=s_axilite port=sigma2 bundle=CTRL
	
	float sigma2_local = sigma2;
	
	MyComplex x_hat_1[8];
	MyComplex y_local[8];
	MyComplex H_local[64];
	MyComplex v_tb_local[256];
	data_local(H_local, y_local, v_tb_local, H_real, H_imag, y_real, y_imag, v_tb_real, v_tb_imag);

    /********************变量定义*******************/
    like_float alpha;
	like_float dqam;
	like_float step_size;
	int i; int j; int k; int id;int x_init[8];int offset = 0;
    int transA = 1;  // CblasConjTrans 的等效值，表示共轭转置
	int transB = 0;  // CblasNoTrans 的等效值，表示不转置
	like_float log_pacc;
	like_float r_norm;/*the norm of residual vector*/
	like_float r_norm_survivor;
	like_float lr;/*learning rate*/
	like_float r_norm_prop;/*new norm*/
	like_float p_acc;
	like_float p_uni[10];
	MyComplex HH_H[64];
	MyComplex grad_preconditioner[64];
	MyComplex constellation_norm[16];/*depend on 2^mu*/
	MyComplex x_mmse[8];
	MyComplex covar[64];
	MyComplex pmat[64];
	MyComplex x_survivor[8];
    MyComplex r[8];
    MyComplex pr_prev[8];
    MyComplex z_grad[8];
    MyComplex v[8];
    MyComplex z_prop[8];
    MyComplex x_prop[8];
    MyComplex r_prop[8];
    MyComplex temp_NtNt[64];
    MyComplex temp_NtNr[64];
    MyComplex temp_1;
    MyComplex _temp_1;
    MyComplex temp_Nt[8];
    MyComplex temp_Nr[8];
    MyComplex sigma2eye[64]; 
    /****************************核心迭代计算前部分数据处理*************************/
    /*定义发送符号之间最小距离的一半，是星座点经过归一化处理后的结果*/
    get_dqam_hw(dqam);
    /*初始化constellation_norm*/
    my_complex_copy_hw_1(_16QAM_Constellation_hw, 1, constellation_norm, 1);
	my_complex_scal_hw_1(dqam, constellation_norm, 1);
    /*生成部分中间量供后续使用*/
    c_eye_generate_hw(sigma2eye, sigma2_local / (float)(dqam * dqam));
	/*二阶梯度下降，计算grad_preconditioner(梯度更新的预条件矩阵)*/
    grad_preconditioner_updater_hw(H_local, transA, transB, HH_H, sigma2eye, grad_preconditioner, Ntr_1, Ntr_1);
    like_float exponent = like_float(1) / like_float(3); // 避免浮点字面值隐式转换
	like_float exponent_1 = (like_float)Ntr_1 / like_float(8);
    alpha = like_float(1) / hls::pow<64,32>((ap_fixed<64,32>)exponent, (ap_fixed<64,32>)exponent_1);
    c_eye_generate_hw(covar, 1.0f); //源代码有的covar，似乎并未使用
    /*For learning rate line search */
    learning_rate_line_search_hw(lr_approx_1, H_local, transA, transB, grad_preconditioner, Ntr_1, Ntr_1, temp_NtNr, pmat);
    /*x的初始化*/
    x_initialize_hw(mmse_init_1, sigma2eye, Ntr_1, Ntr_1, sigma2_local, HH_H, temp_NtNt, H_local, transA, transB, y_local, temp_Nt, x_mmse, dqam, x_hat_1, x_init, constellation_norm);
    /*计算剩余向量r=y-Hx*/
    r_hw(H_local, transB, transA, x_hat_1, Ntr_1, Ntr_1, r, y_local);
    /*计算剩余向量的范数（就是模值）*/
    r_cal_hw(r, transA, transB, Ntr_1, Ntr_1, &temp_1, x_hat_1, x_survivor, r_norm, r_norm_survivor);
    /*确定最优学习率*/
    lr_hw(lr_approx_1, pmat, transB, transA, r, Ntr_1, Ntr_1, pr_prev, &temp_1, &_temp_1, lr);
    /*补偿初始化*/
    step_size = max(dqam, like_float(hls::sqrt(r_norm / (like_float)Ntr_1))) * alpha;
    /****************************核心迭代计算*********************************/
    for (k = 0; k < iter_1; k++)
	{
        /*更新梯度 z_grad = xhat + lr * (grad_preconditioner @ (AH @ r))*/
        z_grad_hw(H_local, transA, transB, temp_Nt, grad_preconditioner, z_grad, lr, x_hat_1, r);
        /*加入高斯随机扰动*/
        gauss_add_hw(v, v_tb_local, offset, step_size, z_grad, z_prop);
        /*将梯度映射到QAM星座点中 x_prop = constellation_norm[np.argmin(abs(z_prop * ones - constellation_norm), axis=2)].reshape(-1, nt, 1) */
        map_hw(dqam, z_prop, x_prop);
        /*计算新的残差范数 calculate residual norm of the proposal*/
        r_newnorm_hw(H_local, transB, x_prop, transA, temp_Nr, y_local, r_prop, &temp_1, r_norm_prop);
        /*update the survivor*/
        survivor_hw(r_norm_survivor, r_norm_prop, x_prop, x_survivor);
        /*acceptance test＆update GD learning rate＆update random walk size*/
        acceptance_hw(transB, transA, r_norm_prop, r_norm, log_pacc, p_acc, p_uni, x_prop , x_hat_1, r_prop, r, pmat, pr_prev, &temp_1, &_temp_1, lr, step_size, dqam, alpha);
    }
    /****************************迭代结束x_survivor写入输出口*********************************/
    out_hw(x_survivor, 1, x_hat_real, x_hat_imag, 1);
}