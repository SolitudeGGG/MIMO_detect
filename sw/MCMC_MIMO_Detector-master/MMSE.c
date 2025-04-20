#include "MMSE.h"
#include "MyComplex.h"

// MyComplex* mmse_Nt;
// MyComplex* mmse_NtNt;
// MyComplex* sigma2eye;

void MMSE_Init(int Nt, int Nr)
{
	// mmse_Nt = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	// mmse_NtNt = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	// sigma2eye = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
}

float MMSE_detect(int Nt, int Nr, MyComplex* y, MyComplex* H, MyComplex* x_hat, float sigma2)
{
	int transA = 1;  // CblasConjTrans �ĵ�Чֵ����ʾ����ת��
	int transB = 0;  // CblasNoTrans �ĵ�Чֵ����ʾ��ת��
	// c_matmultiple(H, transA, y, transB, Nr, Nt, Nr, 1, mmse_Nt);
	// c_eye_generate(sigma2eye, Nt, sigma2 / 2.0f);
	// c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, mmse_NtNt);
	// my_complex_add(Nt * Nt, mmse_NtNt, sigma2eye, mmse_NtNt);
	// Inverse_LU(mmse_NtNt, Nt, Nt);
	// c_matmultiple(mmse_NtNt, transB, mmse_Nt, transB, Nt, Nt, Nt, 1, x_hat);

	return 0;
}

void MMSE_free()
{
	// free(mmse_Nt);
	// free(mmse_NtNt);
	// free(sigma2eye);
}