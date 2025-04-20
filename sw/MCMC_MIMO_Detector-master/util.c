#include "util.h"
#include "MyComplex.h"
#include "time.h"

#define RAND_MAX 0x7fffffff // ��������������?�õ���
#define M_PI 3.14159265358979323846 // pi

FILE* f;


// ���ɶԽǾ���
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

// ����ȫ1����
void c_ones_generate(MyComplex* Mat, int row, int col)
{
	int i;
	for (i = 0; i < row * col; i++)
	{
		Mat[i].real = 1.0f; Mat[i].imag = 0.0f;
	}
}

// �����˷�
MyComplex complex_multiply(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

// �����ӷ�
MyComplex complex_add(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// ��������
MyComplex complex_subtract(MyComplex a, MyComplex b) {
    MyComplex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// ��������
MyComplex complex_divide(MyComplex a, MyComplex b) {
    MyComplex result;
    float denominator = b.real * b.real + b.imag * b.imag;
    result.real = (a.real * b.real + a.imag * b.imag) / denominator;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
    return result;
}

// �����ľ���ֵ
float complex_abs(MyComplex a) {
	return sqrt(a.real * a.real + a.imag * a.imag);
}

// ��ʼ����������
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
	//����д���ǽ����鿴�ɾ��������ж�ȡ��i�е�j����i*col+j��
}

// ����˷�����������?���е��ã�
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

// ������������
void Inverse_LU(MyComplex* A, int row, int col)
{
	MyComplex* L;
	MyComplex* U;
	MyComplex* L_Inverse;
	MyComplex* U_Inverse;
	L = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	U = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	L_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));
	U_Inverse = (MyComplex*)malloc(row * col * sizeof(MyComplex));

	initMatrix(L, row, col);
	initMatrix(U, row, col);
	initMatrix(L_Inverse, row, col);
	initMatrix(U_Inverse, row, col);

	int i, j, k, t;
	MyComplex s;
	for (i = 0; i < row; i++)//�Խ�Ԫ��
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
	for (i = 0; i < col; i++) //���������ڰ��մ��µ��ϣ�����u�������?
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

// ��������ӷ�?
void my_complex_add(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		r[i].real = a[i].real + b[i].real; 
		r[i].imag = a[i].imag + b[i].imag;  
	}
}

// �����������?
void my_complex_sub(const int n, const MyComplex a[], const MyComplex b[], MyComplex r[]) {
	for (int i = 0; i < n; i++) {
		r[i].real = a[i].real - b[i].real;  
		r[i].imag = a[i].imag - b[i].imag; 
	}
}

// ����ȡ����
MyComplex complex_conjugate(MyComplex a) {
	MyComplex result;
	result.real = a.real;   // ʵ������
	result.imag = -a.imag;  // �鲿ȡ��
	return result;
}

// ���Ƹ�������
void my_complex_copy(const int N, const MyComplex* X, const int incX, MyComplex* Y, const int incY) {
	for (int i = 0, j = 0; i < N; i++, j++) {
		Y[j * incY].real = X[i * incX].real;  // ����ʵ��
		Y[j * incY].imag = X[i * incX].imag;  // �����鲿
	}
}

// ���Ÿ�������
void my_complex_scal(const int N, const float alpha, MyComplex* X, const int incX) {
	for (int i = 0; i < N; i++) {
		X[i * incX].real *= alpha;  // ʵ������ alpha
		X[i * incX].imag *= alpha;  // �鲿���� alpha
	}
}

// �����µĳ˷�
void c_matmultiple(MyComplex* matA, int transA, MyComplex* matB, int transB, int ma, int na, int mb, int nb, MyComplex* res) {
	int m, n, k;

	// ����ת�õ��������? m, n, k
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

	// ����������Ԫ��
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			MyComplex sum = { 0.0f, 0.0f };  // ��ʼ��Ϊ������
			MyComplex temp = { 0.0f, 0.0f };
			for (int l = 0; l < k; l++) {
				MyComplex a_element, b_element;

				// ��ȡ����A��Ԫ��
				if (transA == 1) {  // ����ת��
					a_element = matA[l * na + i];  // ��ȡ����A[l, i]
					a_element = complex_conjugate(a_element);  // ��Ԫ��ȡ����
				}
				else if (transA == 0) {  // ��ͨת��
					a_element = matA[i * na + l];  // ��ȡ����A[i, l]
				}

				// ��ȡ����B��Ԫ��
				if (transB == 0) {  // No Transpose
					// �����?1-1
					b_element = matB[l * nb + j];  // ��ȡ����B[l, j]
				}
				else {  // ת��
					b_element = matB[j * nb + l];  // ��ȡ����B[j, l]
					b_element = complex_conjugate(b_element);  // ��Ԫ��ȡ����
				}

				// �����˷����ۼ�
				temp = complex_multiply(a_element, b_element);
				sum = complex_add(sum, temp);
			}

			// �������ֵ������C��res��
			res[i * nb + j] = sum;
		}
	}
}

// ��ȡ�Ѿ����ɵĸ�˹�����������?
void read_gaussian_data(const char* filename, MyComplex* array, int n, int offset) {
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening file!\n");
		return;  // ������ļ�?
	}

	int i;
	float temp_real, temp_imag;

	// �����ļ��е�ǰoffset�����ݣ�������ȡ
	for (i = 0; i < offset; i++) {
		// ��ȡÿһ�еĸ�ʽ "Real: <value>, Imaginary: <value>"
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data at offset %d!\n", i);
			fclose(f);
			return;  // �����ȡ����?
		}
	}

	// ��ȡn�����ݵ�
	for (i = 0; i < n; i++) {
		if (fscanf(f, "Real: %f, Imaginary: %f\n", &temp_real, &temp_imag) != 2) {
			printf("Error reading data!\n");
			fclose(f);
			return;  // �����ѳɹ���ȡ����������
		}

		// ����ȡ��ʵ�����鲿���ݴ�������
		array[i].real = temp_real / (float)sqrt(2.0f);
		array[i].imag = temp_imag / (float)sqrt(2.0f);
	}

	fclose(f);
}

// ���ɾ��ȷֲ��������?
void generateUniformRandoms_int(int Nt, int* x_init, int mu) {
	// ������������ӣ�ʹ�õ�ǰʱ��?
	srand((unsigned int)time(NULL));

	// ���� 2^mu
	double upper_bound = pow(2, mu);  // ���� 2^mu

	// ���� Nt �������?
	for (int i = 0; i < Nt; ++i) {
		// ʹ�� rand() ����һ�� [0, RAND_MAX] ��Χ�ڵ�������Ȼ�����ŵ� [0, upper_bound]
		double rand_val = ((double)rand() / RAND_MAX) * upper_bound;

		// �����ɵ��������?�� x_init ������
		x_init[i] = (int)rand_val;
	}
}

// ���ɾ��ȷֲ��������?
void generateUniformRandoms_float(int Nt, float* p_uni) {
	// ������������ӣ�ʹ�õ�ǰʱ��?
	srand((unsigned int)time(NULL));

	// ���� Nt �����������Χ��? [0.0, 1.0] ֮��
	for (int i = 0; i < Nt; ++i) {
		// ʹ�� rand() ����һ�� [0, RAND_MAX] ��Χ��������Ȼ�����ŵ� [0.0, 1.0]
		p_uni[i] = (float)rand() / (float)RAND_MAX;
	}
}

// ����������Сֵ
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

// �����������������в�����?�صĸ���
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

// ������ӳ��
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