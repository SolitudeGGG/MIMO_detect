#include "util_1.h"
#include "MyComplex_1.h"
#include <stdio.h>
#include "MHGD_accel_1.h"
#include "hls_math.h"

#define LCG_A 1664525
#define LCG_C 1013904223
//#define LCG_M 4294967296  // 2^32
#define LIMIT_MAX 0x7fffffff // ???????????????????
#define LIMIT_MAX_F 2147483647.0f
#define M_PI 3.14159265358979323846 // pi
//#define  __SYNTHESIS__
//#define NO_SYNTH

FILE* f;

// 线性同余生成器
unsigned int lcg_rand() {
    unsigned int lcg_seed = 1;
	lcg_seed = (LCG_A * lcg_seed + LCG_C) % LIMIT_MAX;;
	
    return lcg_seed;
}
like_float lcg_rand_1() {
    unsigned int lcg_seed = 1;
	lcg_seed = (LCG_A * lcg_seed + LCG_C) % LIMIT_MAX;
	//lcg_seed = (LCG_A * lcg_seed + LCG_C) % LIMIT_MAX;
    return (like_float)(lcg_seed / LIMIT_MAX_F);
}

// ?????????
void c_eye_generate(MyComplex* Mat, int n, float val)
{
	like_float val_1 = val;
	int i;
	for (i = 0; i < n * n; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_2_max min=Nt_2_min
		//printf("i=%d",i);
		Mat[i].real = (like_float)0.0;
		Mat[i].imag = (like_float)0.0;
	}
	for (i = 0; i < n; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		Mat[i * n + i].real = val_1;
	}
}

// ?????1????
void c_ones_generate(MyComplex* Mat, int row, int col)
{
	int i;
	for (i = 0; i < row * col; i++)
	{
		Mat[i].real = (like_float)1.0; Mat[i].imag = (like_float)0.0;
	}
}

// ???????
MyComplex complex_multiply(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

// ???????
MyComplex complex_add(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// ????????
MyComplex complex_subtract(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// ????????
MyComplex complex_divide(MyComplex a, MyComplex b) {
    MyComplex result;
    like_float denominator = b.real * b.real + b.imag * b.imag;
    result.real = (a.real * b.real + a.imag * b.imag) / denominator;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
    return result;
}

// ??????????
like_float complex_abs(MyComplex a) {
	return (like_float)hls::sqrt(a.real * a.real + a.imag * a.imag);
}

// ?????????????
void initMatrix(MyComplex* A, int row, int col)
{
	int i, j;
	for (i = 0; i < row; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		for (j = 0; j < col; j++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			A[i * col + j].real = (like_float)0.0;
			A[i * col + j].imag = (like_float)0.0;
		}
	}
	//????��???????????????????��????i?��?j????i*col+j??
}

// ???????????????��???��????
void MulMatrix(const MyComplex* A, const MyComplex* B, MyComplex* C, int row, int AB_rc, int col)
{
	int i, j, k;
	MyComplex sum, tmp;
	for (i = 0; i < row; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		for (j = 0; j < col; j++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			sum.real = (like_float)0.0;
			sum.imag = (like_float)0.0;
			for (k = 0; k < AB_rc; k++)
			{
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				sum = complex_add(sum, complex_multiply(A[i * AB_rc + k], B[k * col + j]));
			}
			C[i * col + j].real = sum.real;
			C[i * col + j].imag = sum.imag;
		}
	}
}

// ????????????
void Inverse_LU(MyComplex* A, int row, int col)
{
	#ifdef NO_SYNTH
	MyComplex* L;
	MyComplex* U;
	MyComplex* L_Inverse;
	MyComplex* U_Inverse;
	L = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	U = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	L_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	U_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	#else
	MyComplex _L[64];//row * col
	MyComplex _U[64];//row * col
	MyComplex _L_Inverse[64];//row * col
	MyComplex _U_Inverse[64];//row * col

	MyComplex* L = &_L[0];
	MyComplex* U = &_U[0];
	MyComplex* L_Inverse = &_L_Inverse[0];
	MyComplex* U_Inverse = &_U_Inverse[0];
	#endif

	initMatrix(L, row, col);
	initMatrix(U, row, col);
	initMatrix(L_Inverse, row, col);
	initMatrix(U_Inverse, row, col);

	int i, j, k, t;
	MyComplex s;
	#pragma HLS PIPELINE II=1
	for (i = 0; i < row; i++)//??????
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS UNROLL
		L[i * col + i].real = (like_float)1.0;
		L[i * col + i].imag = (like_float)0.0;
	}
	for (j = 0; j < col; j++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS UNROLL
		U[j] = A[j];
	}
	for (i = 1; i < col; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS UNROLL
		L[i * col] = complex_divide(A[i * col], U[0]);

	}
	for (k = 1; k < row; k++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS PIPELINE II=1
		for (j = k; j < col; j++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			#pragma HLS PIPELINE II=1
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (t = 0; t < k; t++) {
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				#pragma HLS UNROLL
				s = complex_add(s, complex_multiply(L[k * col + t], U[t * col + j]));
			}
			U[k * col + j] = complex_subtract(A[k * col + j], s);
		}
		for (i = k; i < col; i++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			#pragma HLS PIPELINE II=1
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (t = 0; t < k; t++)
			{
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				#pragma HLS UNROLL
				s = complex_add(s, complex_multiply(L[i * col + t], U[t * col + k]));
			}
			L[i * col + k] = complex_divide(complex_subtract(A[i * col + k], s), U[k * col + k]);
		}
	}

	for (i = 0; i < row; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS UNROLL
		L_Inverse[i * col + i].imag = (like_float)0.0;
		L_Inverse[i * col + i].real = (like_float)1.0;
	}
	for (j = 0; j < col; j++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS PIPELINE II=1
		for (i = j + 1; i < row; i++)
		{
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			#pragma HLS PIPELINE II=1
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (k = j; k < i; k++) {
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				#pragma HLS UNROLL
				s = complex_add(s, complex_multiply(L[i * col + k], L_Inverse[k * col + j]));
			}
			L_Inverse[i * col + j].real = (-1) * complex_multiply(L_Inverse[j * col + j], s).real;
			L_Inverse[i * col + j].imag = (-1) * complex_multiply(L_Inverse[j * col + j], s).imag;
		}
	}
	MyComplex di;
	di.real = (like_float)1.0;
	di.imag = (like_float)0.0;
	for (i = 0; i < col; i++) //????????????????????????u???????
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS UNROLL
		U_Inverse[i * col + i] = complex_divide(di, U[i * col + i]);

	}
	for (j = 0; j < col; j++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS PIPELINE II=1
		for (i = j - 1; i >= 0; i--)
		{
			#pragma HLS PIPELINE II=1
			s.imag = (like_float)0.0;
			s.real = (like_float)0.0;
			for (k = i + 1; k <= j; k++)
			{
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				#pragma HLS UNROLL
				s = complex_add(s, complex_multiply(U[i * col + k], U_Inverse[k * col + j]));

			}
			s.imag = -s.imag;
			s.real = -s.real;
			U_Inverse[i * col + j] = complex_divide(s, U[i * col + i]);
		}
	}
	MulMatrix(U_Inverse, L_Inverse, A, row, row, col);
}

// ??????????
void my_complex_add(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		r[i].real = a[i].real + b[i].real; 
		r[i].imag = a[i].imag + b[i].imag;  
	}
}

// ???????????
void my_complex_sub(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		#pragma HLS LOOP_TRIPCOUNT max=Nr_max min=Nr_min
		r[i].real = a[i].real - b[i].real;  
		r[i].imag = a[i].imag - b[i].imag; 
	}
}

// ?????????
MyComplex complex_conjugate(MyComplex a) {
	MyComplex result;
	result.real = a.real;   // ???????
	result.imag = -a.imag;  // ?��???
	return result;
}

// ???????????
void my_complex_copy(const int N, const MyComplex* X, const int incX, MyComplex* Y, const int incY) {
	for (int i = 0, j = 0; i < N; i++, j++) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		Y[j * incY].real = X[i * incX].real;  // ???????
		Y[j * incY].imag = X[i * incX].imag;  // ?????��
	}
}

// ???????????
void my_complex_scal(const int N, const like_float alpha, MyComplex* X, const int incX) {
	//const like_float alpha_1 = alpha;
	for (int i = 0; i < N; i++) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		X[i * incX].real *= alpha;  // ??????? alpha
		X[i * incX].imag *= alpha;  // ?��???? alpha
	}
}

// ?????????
void c_matmultiple(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res) {
	int m, n, k;

	// ??????????????? m, n, k
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

	// ?????????????
	for (int i = 0; i < m; i++) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS PIPELINE II=1
		for (int j = 0; j < n; j++) {
			#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
			MyComplex sum =  { (like_float)0.0, (like_float)0.0 };  // ????????????
			MyComplex temp = { (like_float)0.0, (like_float)0.0 };
			for (int l = 0; l < k; l++) {
				#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
				#pragma HLS UNROLL
				MyComplex a_element;
                MyComplex b_element;

				// ???????A?????
				if (transA == 1) {  // ???????
					a_element = matA[l * na + i];  // ???????A[l, i]
					a_element = complex_conjugate(a_element);  // ??????????
				}
				else if (transA == 0) {  // ??????
					a_element = matA[i * na + l];  // ???????A[i, l]
				}else {
                // 处理非法 transA 值（可选：报错或断言）
                // 例如：fprintf(stderr, "Invalid transA value: %d\n", transA);
                // 此处暂时赋默认值，避免未初始化警告
                a_element = (MyComplex){(like_float)0.0, (like_float)0.0};
                }

				// ???????B?????
				if (transB == 0) {  // No Transpose
					// ?????1-1
					b_element = matB[l * nb + j];  // ???????B[l, j]
				}
				else {  // ???
					b_element = matB[j * nb + l];  // ???????B[j, l]
					b_element = complex_conjugate(b_element);  // ??????????
				}
                

				// ????????????
				temp = complex_multiply(a_element, b_element);
				sum = complex_add(sum, temp);
			}

			// ??????????????C??res??
			res[i * nb + j] = sum;
		}
	}
}

// ????????????????????????
void read_gaussian_data(const char* filename, MyComplex* array, int n, int offset) {
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening file!\n");
		return;  // ????????
	}
	int i;
	float temp_real, temp_imag;
	// 跳过前offset个数据
	for (i = 0; i < offset; i++) {
		// ??????��??? "Real: <value>, Imaginary: <value>"
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			
			printf("Error reading data at offset %d!\n", i);
			fclose(f);
			
			return;  // ??????????
		}
	}
	
	// ???n???????
	for (i = 0; i < n; i++) {
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data!\n");
			fclose(f);
			return;  // ?????????????????????
		}
		// ?????????????��???????????
		array[i].real = temp_real / (float)sqrt(2.0f);
		array[i].imag = temp_imag / (float)sqrt(2.0f);
	}
	fclose(f);
}

// ????????????????
void generateUniformRandoms_int(int Nt, int* x_init, int mu) {
	// 设置随机数种子，使用当前时间
	//srand((unsigned int)time(NULL));

	// 计算 2^mu
	double upper_bound = pow(2, mu);  // 计算 2^mu

	// 生成 Nt 个随机数
	for (int i = 0; i < Nt; ++i) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		// 使用 rand() 生成一个 [0, LIMIT_MAX] 范围内的整数，然后缩放到 [0, upper_bound]
		double rand_val = ((double)lcg_rand() / LIMIT_MAX) * upper_bound;


		// 将生成的随机数存储到 x_init 数组中
		x_init[i] = (int)rand_val;
	}
}

// ????????????????
void generateUniformRandoms_float(int Nt, like_float p_uni[9]) {//一个随机数生成的函数 其中lcg_rand用于生成随机数，/LIMIT_MAX用于限幅
	for (int i = 0; i < Nt; ++i) {
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		#pragma HLS PIPELINE II=1
		p_uni[i] = lcg_rand_1();
	}
}

// ??????????��?
int argmin(float* array, int n)
{
	int i = 0; int best_id = 0;
	like_float minval = (like_float)array[0];
	for (i = 1; i < n; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		if (minval > (like_float)array[i])
		{
			minval = (like_float)array[i];
			best_id = i;
		}
	}
	return best_id;
}

// ?????????????????��???????????
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

// ?????????
void map(int mu, int Nt, like_float& dqam, MyComplex* x, MyComplex* x_hat)
{
	int i;
	for (i = 0; i < Nt; i++)
	{
		#pragma HLS LOOP_TRIPCOUNT max=Nt_max min=Nt_min
		x_hat[i].real = ((like_float)2.0 * (like_float)hls::floor(x[i].real / dqam / (like_float)2.0) + (like_float)1.0);
		x_hat[i].imag = ((like_float)2.0 * (like_float)hls::floor(x[i].imag / dqam / (like_float)2.0) + (like_float)1.0);
		switch (mu)
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
	my_complex_scal(Nt, dqam, x_hat, 1);
}