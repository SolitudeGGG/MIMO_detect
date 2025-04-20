#include "MyComplex.h"
#include "hls_math.h"
#include "MHGD_accel_sw.h"
#include <stdio.h>
#include "string.h"
#include "time.h"
#include "math.h"
#include <stdlib.h>

MyComplex QPSK_Constellation[4] = {{-1,-1},{-1,1},{1,-1},{1,1}};
MyComplex _16QAM_Constellation[16] = {{-3.0,-3.0},{-3.0,-1.0},{-3.0,3.0},{-3.0,1.0},{-1.0,-3.0},{-1.0,-1.0},{-1.0,3.0},{-1.0,1.0},{3.0,-3.0},{3.0,-1.0},{3.0,3.0},{3.0,1.0},{1.0,-3.0},{1.0,-1.0},{1.0,3.0},{1.0,1.0}};
MyComplex _64QAM_Constellation[64] = {{-7.0, 7.0},  {-5.0, 7.0},  {-3.0, 7.0}, {-1.0, 7.0}, 
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
void read_gaussian_data(const char* filename, MyComplex* array, int n, int offset) {
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening file!\n");
		return;  // 锟斤拷锟斤拷锟斤拷募锟�
	}

	int i;
	float temp_real, temp_imag;

	// 锟斤拷锟斤拷锟侥硷拷锟叫碉拷前offset锟斤拷锟斤拷锟捷ｏ拷锟斤拷锟斤拷锟斤拷取
	for (i = 0; i < offset; i++) {
		// 锟斤拷取每一锟叫的革拷式 "Real: <value>, Imaginary: <value>"
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data at offset %d!\n", i);
			fclose(f);
			return;  // 锟斤拷锟斤拷锟饺★拷锟斤拷锟�
		}
	}

	// 锟斤拷取n锟斤拷锟斤拷锟捷碉拷
	for (i = 0; i < n; i++) {
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data!\n");
			fclose(f);
			return;  // 锟斤拷锟斤拷锟窖成癸拷锟斤拷取锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷
		}

		// 锟斤拷锟斤拷取锟斤拷实锟斤拷锟斤拷锟介部锟斤拷锟捷达拷锟斤拷锟斤拷锟斤拷
		array[i].real = temp_real / (float)sqrt(2.0f);
		array[i].imag = temp_imag / (float)sqrt(2.0f);
	}

	fclose(f);
}
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

/*QPSK���*/
void QPSK_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex temp[4];
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
}

/*16QAM���*/
void _16QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex temp[16];
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
}

/*64QAM���*/
void _64QAM_Demodulation(MyComplex* x_hat, int Nt, int* bits_demod)
{
	int i = 0; int j = 0; int best_id = 0;
	MyComplex temp[64];
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
}
int argmin(float* array, int n)
{
	int i = 0; int best_id = 0;
	float minval = array[0];
	for (i = 1; i < n; i++)
	{
		if (minval > array[i])
		{
			minval = array[i];
			best_id = i;
		}
	}
	return best_id;
}
int unequal_times(int* array1, int* array2, int n)
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
MyComplex complex_conjugate(MyComplex a) {
	MyComplex result;
	result.real = a.real;   // 实锟斤拷锟斤拷锟斤拷
	result.imag = -a.imag;  // 锟介部取锟斤拷
	return result;
}
void c_eye_generate(MyComplex* Mat, int n, float val)
{
	int i;
	for (i = 0; i < n * n; i++)
	{
		Mat[i].real = 0.0f;
		Mat[i].imag = 0.0f;
	}
	for (i = 0; i < n; i++)
	{
		Mat[i * n + i].real = val;
	}
}
void c_matmultiple(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res) {
	int m, n, k;

	// 锟斤拷锟斤拷转锟矫碉拷锟斤拷锟斤拷锟斤拷锟� m, n, k
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

	// 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷元锟斤拷
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			MyComplex sum = { 0.0f, 0.0f };  // 锟斤拷始锟斤拷为锟斤拷锟斤拷锟斤拷
			MyComplex temp = { 0.0f, 0.0f };
			for (int l = 0; l < k; l++) {
				MyComplex a_element, b_element;

				// 锟斤拷取锟斤拷锟斤拷A锟斤拷元锟斤拷
				if (transA == 1) {  // 锟斤拷锟斤拷转锟斤拷
					a_element = matA[l * na + i];  // 锟斤拷取锟斤拷锟斤拷A[l, i]
					a_element = complex_conjugate(a_element);  // 锟斤拷元锟斤拷取锟斤拷锟斤拷
				}
				else if (transA == 0) {  // 锟斤拷通转锟斤拷
					a_element = matA[i * na + l];  // 锟斤拷取锟斤拷锟斤拷A[i, l]
				}

				// 锟斤拷取锟斤拷锟斤拷B锟斤拷元锟斤拷
				if (transB == 0) {  // No Transpose
					// 锟斤拷锟斤拷锟�1-1
					b_element = matB[l * nb + j];  // 锟斤拷取锟斤拷锟斤拷B[l, j]
				}
				else {  // 转锟斤拷
					b_element = matB[j * nb + l];  // 锟斤拷取锟斤拷锟斤拷B[j, l]
					b_element = complex_conjugate(b_element);  // 锟斤拷元锟斤拷取锟斤拷锟斤拷
				}

				// 锟斤拷锟斤拷锟剿凤拷锟斤拷锟桔硷拷
				temp = complex_multiply(a_element, b_element);
				sum = complex_add(sum, temp);
			}

			// 锟斤拷锟斤拷锟斤拷锟街碉拷锟斤拷锟斤拷锟紺锟斤拷res锟斤拷
			res[i * nb + j] = sum;
		}
	}
}
void my_complex_add(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		r[i].real = a[i].real + b[i].real; 
		r[i].imag = a[i].imag + b[i].imag;  
	}
}
void initMatrix(MyComplex* A, int row, int col)
{
	int i, j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			A[i * col + j].real = 0.0;
			A[i * col + j].imag = 0.0;
		}
	}
	//锟斤拷锟斤拷写锟斤拷锟角斤拷锟斤拷锟介看锟缴撅拷锟斤拷锟斤拷锟斤拷锟叫讹拷取锟斤拷i锟叫碉拷j锟斤拷锟斤拷i*col+j锟斤拷
}
MyComplex complex_divide(MyComplex a, MyComplex b) {
    MyComplex result;
    float denominator = b.real * b.real + b.imag * b.imag;
    result.real = (a.real * b.real + a.imag * b.imag) / denominator;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
    return result;
}
MyComplex complex_multiply(MyComplex a, MyComplex b) {
    MyComplex result;
    float local_temp_1;
	float local_temp_2;
	float local_temp_3;
	float local_temp_4;
    local_temp_1 = a.real * b.real;
    local_temp_2 = a.imag * b.imag;
    local_temp_3 = a.real * b.imag;
    local_temp_4 = a.imag * b.real;
    result.real = local_temp_1 - local_temp_2;
    result.imag = local_temp_3 + local_temp_4;
    return result;
}
MyComplex complex_add(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}
MyComplex complex_subtract(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}
void my_complex_sub(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		r[i].real = a[i].real - b[i].real;  
		r[i].imag = a[i].imag - b[i].imag; 
	}
}
void MulMatrix(const MyComplex* A, const MyComplex* B, MyComplex* C, int row, int AB_rc, int col)
{
	int i, j, k;
	MyComplex sum, tmp;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{

			sum.real = 0.0;
			sum.imag = 0.0;
			for (k = 0; k < AB_rc; k++)
			{
				sum = complex_add(sum, complex_multiply(A[i * AB_rc + k], B[k * col + j]));
			}
			C[i * col + j].real = sum.real;
			C[i * col + j].imag = sum.imag;
		}
	}
}
void Inverse_LU(MyComplex* A, int row, int col)
{
	//MyComplex* L;
	//MyComplex* U;
	//MyComplex* L_Inverse;
	//MyComplex* U_Inverse;
	//L = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	//U = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	//L_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	//U_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));

    MyComplex L[row * col];//row * col
	MyComplex U[row * col];//row * col
	MyComplex L_Inverse[row * col];//row * col
	MyComplex U_Inverse[row * col];//row * col

	// MyComplex* L = &_L[0];
	// MyComplex* U = &_U[0];
	// MyComplex* L_Inverse = &_L_Inverse[0];
	// MyComplex* U_Inverse = &_U_Inverse[0];

	initMatrix(L, row, col);
	initMatrix(U, row, col);
	initMatrix(L_Inverse, row, col);
	initMatrix(U_Inverse, row, col);

	int i, j, k, t;
	MyComplex s;
	for (i = 0; i < row; i++)//锟皆斤拷元锟斤拷
	{
		L[i * col + i].real = 1.0;
		L[i * col + i].imag = 0.0;
	}
	for (j = 0; j < col; j++)
	{
		U[j] = A[j];
	}
	for (i = 1; i < col; i++)
	{
		L[i * col] = complex_divide(A[i * col], U[0]);
	}
	for (k = 1; k < row; k++)
	{
		for (j = k; j < col; j++)
		{
			s.imag = 0.0;
			s.real = 0.0;
			for (t = 0; t < k; t++) {
				s = complex_add(s, complex_multiply(L[k * col + t], U[t * col + j]));
			}
			U[k * col + j] = complex_subtract(A[k * col + j], s);
		}
		for (i = k; i < col; i++)
		{
			s.imag = 0.0;
			s.real = 0.0;
			for (t = 0; t < k; t++)
			{
				s = complex_add(s, complex_multiply(L[i * col + t], U[t * col + k]));
			}
			L[i * col + k] = complex_divide(complex_subtract(A[i * col + k], s), U[k * col + k]);
		}
	}

	for (i = 0; i < row; i++)
	{
		L_Inverse[i * col + i].imag = 0.0;
		L_Inverse[i * col + i].real = 1.0;
	}
	for (j = 0; j < col; j++)
	{
		for (i = j + 1; i < row; i++)
		{
			s.imag = 0.0;
			s.real = 0.0;
			for (k = j; k < i; k++) {
				s = complex_add(s, complex_multiply(L[i * col + k], L_Inverse[k * col + j]));
			}
			L_Inverse[i * col + j].real = (-1) * complex_multiply(L_Inverse[j * col + j], s).real;
			L_Inverse[i * col + j].imag = (-1) * complex_multiply(L_Inverse[j * col + j], s).imag;
		}
	}
	MyComplex di;
	di.real = 1.0;
	di.imag = 0.0;
	for (i = 0; i < col; i++) //锟斤拷锟斤拷锟斤拷锟斤拷锟节帮拷锟秸达拷锟铰碉拷锟较ｏ拷锟斤拷锟斤拷u锟斤拷锟斤拷锟斤拷锟�
	{
		U_Inverse[i * col + i] = complex_divide(di, U[i * col + i]);

	}
	for (j = 0; j < col; j++)
	{
		for (i = j - 1; i >= 0; i--)
		{
			s.imag = 0.0;
			s.real = 0.0;
			for (k = i + 1; k <= j; k++)
			{
				s = complex_add(s, complex_multiply(U[i * col + k], U_Inverse[k * col + j]));

			}
			s.imag = -s.imag;
			s.real = -s.real;
			U_Inverse[i * col + j] = complex_divide(s, U[i * col + i]);
		}
	}
	MulMatrix(U_Inverse, L_Inverse, A, row, row, col);
}
void my_complex_scal(const int N, const float alpha, MyComplex* X, const int incX) {
	for (int i = 0; i < N; i++) {
		X[i * incX].real *= alpha;  // 实锟斤拷锟斤拷锟斤拷 alpha
		X[i * incX].imag *= alpha;  // 锟介部锟斤拷锟斤拷 alpha
	}
}
void map(int mu, int Nt, float dqam, MyComplex* x, MyComplex* x_hat)
{
	int i;
	for (i = 0; i < Nt; i++)
	{
		x_hat[i].real = (2.0f * floor(x[i].real / dqam / 2.0f) + 1.0f);
		x_hat[i].imag = (2.0f * floor(x[i].imag / dqam / 2.0f) + 1.0f);
		switch (mu)
		{
		case 2:
			x_hat[i].real = (x_hat[i].real > 1.0f ? 1.0f : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < -1.0f ? -1.0f : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > 1.0f ? 1.0f : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < -1.0f ? -1.0f : x_hat[i].imag);
			break;
		case 4:
			x_hat[i].real = (x_hat[i].real > 3.0f ? 3.0f : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < -3.0f ? -3.0f : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > 3.0f ? 3.0f : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < -3.0f ? -3.0f : x_hat[i].imag);
			break;
		case 6:
			x_hat[i].real = (x_hat[i].real > 7.0f ? 7.0f : x_hat[i].real);
			x_hat[i].real = (x_hat[i].real < -7.0f ? -7.0f : x_hat[i].real);
			x_hat[i].imag = (x_hat[i].imag > 7.0f ? 7.0f : x_hat[i].imag);
			x_hat[i].imag = (x_hat[i].imag < -7.0f ? -7.0f : x_hat[i].imag);
			break;
		default:
			break;
		}
	}
	my_complex_scal(Nt, dqam, x_hat, 1);
}
void generateUniformRandoms_int(int Nt, int* x_init, int mu) {
	// 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷樱锟绞癸拷玫锟角笆憋拷锟�
	srand((unsigned int)time(NULL));

	// 锟斤拷锟斤拷 2^mu
	double upper_bound = pow(2, mu);  // 锟斤拷锟斤拷 2^mu

	// 锟斤拷锟斤拷 Nt 锟斤拷锟斤拷锟斤拷锟�
	for (int i = 0; i < Nt; ++i) {
		// 使锟斤拷 rand() 锟斤拷锟斤拷一锟斤拷 [0, RAND_MAX] 锟斤拷围锟节碉拷锟斤拷锟斤拷锟斤拷然锟斤拷锟斤拷锟脚碉拷 [0, upper_bound]
		double rand_val = ((double)rand() / RAND_MAX) * upper_bound;

		// 锟斤拷锟斤拷锟缴碉拷锟斤拷锟斤拷锟斤拷娲�锟斤拷 x_init 锟斤拷锟斤拷锟斤拷
		x_init[i] = (int)rand_val;
	}
}
void my_complex_copy(const int N, const MyComplex* X, const int incX, MyComplex* Y, const int incY) {
	for (int i = 0, j = 0; i < N; i++, j++) {
		Y[j * incY].real = X[i * incX].real;  // 锟斤拷锟斤拷实锟斤拷
		Y[j * incY].imag = X[i * incX].imag;  // 锟斤拷锟斤拷锟介部
	}
}
void generateUniformRandoms_float(int Nt, float* p_uni) {
	// 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷樱锟绞癸拷玫锟角笆憋拷锟�
	srand((unsigned int)time(NULL));

	// 锟斤拷锟斤拷 Nt 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟轿э拷锟� [0.0, 1.0] 之锟斤拷
	for (int i = 0; i < Nt; ++i) {
		// 使锟斤拷 rand() 锟斤拷锟斤拷一锟斤拷 [0, RAND_MAX] 锟斤拷围锟斤拷锟斤拷锟斤拷锟斤拷然锟斤拷锟斤拷锟脚碉拷 [0.0, 1.0]
		p_uni[i] = (float)rand() / (float)RAND_MAX;
	}
}
void out_hw(const MyComplex* X, const int incX, float* Y_real, float* Y_imag, int incY)
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
void MHGD_detect_accel_sw(
    float* x_hat_real, float* x_hat_imag, 
    float* H_real, float* H_imag, 
    float* y_real, float* y_imag, 
    float sigma2,
    float* v_tb_real, float* v_tb_imag
){
    MyComplex x_hat[8];
	MyComplex y[8];
	MyComplex H[64];
	MyComplex v_tb[256];
    int i;
	for(i = 0; i < Ntr_1; i++){
		y[i].real = y_real[i];
		y[i].imag = y_imag[i];
	}
	for (i = 0; i < Ntr_1*Ntr_1; i++)
	{
		H[i].real = H_real[i];
		H[i].imag = H_imag[i];
	}
	for (i = 0; i < Ntr_1*iter_1; i++)
	{
		v_tb[i].real = v_tb_real[i];
		v_tb[i].imag = v_tb_imag[i];
	}
    /********************变量定义*******************/
    float alpha;
	float dqam;
	float step_size;
	int x_init[8] = {1, 11, 3, 6, 13, 11, 15, 9};
	int j; int k; int offset = 0;
    int transA = 1;  // CblasConjTrans 的等效值，表示共轭转置
	int transB = 0;  // CblasNoTrans 的等效值，表示不转置
	float log_pacc;
	float r_norm;/*the norm of residual vector*/
	float r_norm_survivor;
	float lr;/*learning rate*/
	float r_norm_prop;/*new norm*/
	float p_acc;
	float p_uni[10] = {0.0546594001, 0.372195959, 0.999145865, 0.0510859713, 0.411008626, 0.0656750798, 0.0993093923, 0.258695126, 0.443532109, 0.960875571};
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
    MyComplex temp_1[1];
    MyComplex _temp_1[1];
    MyComplex temp_Nt[8];
    MyComplex temp_Nr[8];
    MyComplex sigma2eye[64]; 
	float local_temp_1;
	float local_temp_2;
    dqam = sqrt(1.5f / (pow(2.0f, (float)mu_1) - 1));
    my_complex_copy(pow(2, mu_1), _16QAM_Constellation, 1, constellation_norm, 1);
    my_complex_scal(pow(2, mu_1), dqam, constellation_norm, 1);

    c_eye_generate(sigma2eye, Ntr_1, sigma2 / (dqam * dqam));

    c_matmultiple(H, transA, H, transB, Ntr_1, Ntr_1, Ntr_1, Ntr_1, HH_H);
	my_complex_add(Ntr_1 * Ntr_1, HH_H, sigma2eye, grad_preconditioner);
	Inverse_LU(grad_preconditioner, Ntr_1, Ntr_1);

    alpha = 1.0f / pow((float)Ntr_1 / 8.0f, 1.0f / 3.0f);

    c_eye_generate(covar, Ntr_1, 1.0f);

    if (!lr_approx_1)
	{
		c_matmultiple(H, transB, grad_preconditioner, transB, Ntr_1, Ntr_1, Ntr_1, Ntr_1, temp_NtNr);
		c_matmultiple(temp_NtNr, transB, H, transA, Ntr_1, Ntr_1, Ntr_1, Ntr_1, pmat);
	}
	else
	{
		for (i = 0; i < Ntr_1 * Ntr_1; i++)
		{
			pmat[i].real = 0; pmat[i].imag = 0;
		}
	}

    if (mmse_init_1)
	{
		/*x_mmse = la.inv(AHA + noise_var * np.eye(nt)) @ AH @ y*/
		c_eye_generate(sigma2eye, Ntr_1, sigma2);
		my_complex_add(Ntr_1 * Ntr_1, sigma2eye, HH_H, temp_NtNt);
		Inverse_LU(temp_NtNt, Ntr_1, Ntr_1);
		c_matmultiple(H, transA, y, transB, Ntr_1, Ntr_1, Ntr_1, 1, temp_Nt);
		c_matmultiple(temp_NtNt, transB, temp_Nt, transB, Ntr_1, Ntr_1, Ntr_1, 1, x_mmse);
		/*ӳ�䵽��һ������ͼ�� xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map(mu_1, Ntr_1, dqam, x_mmse, x_hat);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		//generateUniformRandoms_int(Ntr_1, x_init, mu_1);
		for (i = 0; i < Ntr_1; i++)
			x_hat[i] = constellation_norm[x_init[i]];
	}

    c_matmultiple(H, transB, x_hat, transB, Ntr_1, Ntr_1, Ntr_1, 1, r);
	my_complex_sub(Ntr_1, y, r, r);

    c_matmultiple(r, transA, r, transB, Ntr_1, 1, Ntr_1, 1, temp_1);
	r_norm = temp_1[0].real;
	my_complex_copy(Ntr_1, x_hat, 1, x_survivor, 1);
	r_norm_survivor = r_norm;

    if (!lr_approx_1)
	{
		c_matmultiple(pmat, transB, r, transB, Ntr_1, Ntr_1, Ntr_1, 1, pr_prev);
		c_matmultiple(r, transA, pr_prev, transB, Ntr_1, 1, Ntr_1, 1, temp_1);
		c_matmultiple(pr_prev, transA, pr_prev, transB, Ntr_1, 1, Ntr_1, 1, _temp_1);
		lr = temp_1[0].real / _temp_1[0].real;
	}
	else
	{
		lr = 1;
	}

    step_size = fmax(dqam, sqrt(r_norm / (float)Ntr_1)) * alpha;

    for (k = 0; k < iter_1; k++)
    {
        c_matmultiple(H, transA, r, transB, Ntr_1, Ntr_1, Ntr_1, 1, temp_Nt);
		c_matmultiple(grad_preconditioner, transB, temp_Nt, transB, Ntr_1, Ntr_1, Ntr_1, 1, z_grad);
		my_complex_scal(Ntr_1, lr, z_grad, 1);
		my_complex_add(Ntr_1, x_hat, z_grad, z_grad);

        read_gaussian_data("gaussian_random_values.txt", v, Ntr_1, offset);
		offset = offset + Ntr_1;
		my_complex_scal(Ntr_1, step_size, v, 1);
		my_complex_add(Ntr_1, z_grad, v, z_prop);

        map(mu_1, Ntr_1, dqam, z_prop, x_prop);

        c_matmultiple(H, transB, x_prop, transB, Ntr_1, Ntr_1, Ntr_1, 1, temp_Nr);
		my_complex_sub(Ntr_1, y, temp_Nr, r_prop);
		c_matmultiple(r_prop, transA, r_prop, transB, Ntr_1, 1, Ntr_1, 1, temp_1);
		r_norm_prop = temp_1[0].real;

        if (r_norm_survivor > r_norm_prop)
		{
			my_complex_copy(Ntr_1, x_prop, 1, x_survivor, 1);
			r_norm_survivor = r_norm_prop;
		}

        log_pacc = fmin(0, -(r_norm_prop - r_norm));
		p_acc = exp(log_pacc);
		//generateUniformRandoms_float(10, p_uni);

        if (p_acc > p_uni[5])/*������������ʱ��*/
		{
			my_complex_copy(Ntr_1, x_prop, 1, x_hat, 1);
			my_complex_copy(Ntr_1, r_prop, 1, r, 1);
			r_norm = r_norm_prop;
			/*update GD learning rate*/
			if (!lr_approx_1)
			{
				c_matmultiple(pmat, transB, r, transB, Ntr_1, Ntr_1, Ntr_1, 1, pr_prev);
				c_matmultiple(r, transA, pr_prev, transB, Ntr_1, 1, Ntr_1, 1, temp_1);
				c_matmultiple(pr_prev, transA, pr_prev, transB, Ntr_1, 1, Ntr_1, 1, _temp_1);
				lr = temp_1[0].real / _temp_1[0].real;
			}
			/*update random walk size*/
            local_temp_1 = sqrt(r_norm / (float)Ntr_1);    
            local_temp_2 = fmax(dqam, local_temp_1);
			step_size = local_temp_2 * alpha;
		}
    }
    /****************************迭代结束x_survivor写入输出口*********************************/
    out_hw(x_survivor, 1, x_hat_real, x_hat_imag, 1);
}