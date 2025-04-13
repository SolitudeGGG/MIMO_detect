#include "MMSE_1.h"
#include "MyComplex_1.h"
#include "util_1.h"

MyComplex* mmse_Nt;
MyComplex* mmse_NtNt;
//MyComplex* sigma2eye;

void MMSE_Init(int Nt, int Nr)
{
	//MyComplex* sigma2eye;
	mmse_Nt = (MyComplex*)malloc(Nt * sizeof(MyComplex));
	mmse_NtNt = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	//sigma2eye = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
}

float MMSE_detect(int Nt, int Nr, MyComplex* y, MyComplex* H, MyComplex x_hat[8], float sigma2)
{
	MyComplex* sigma2eye;
	sigma2eye = (MyComplex*)malloc(Nt * Nt * sizeof(MyComplex));
	int transA = 1;  // CblasConjTrans 的等效值，表示共轭转置
	int transB = 0;  // CblasNoTrans 的等效值，表示不转置
	c_matmultiple(H, transA, y, transB, Nr, Nt, Nr, 1, mmse_Nt);
	c_eye_generate(sigma2eye, Nt, sigma2 / 2.0f);
	c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, mmse_NtNt);
	my_complex_add(Nt * Nt, mmse_NtNt, sigma2eye, mmse_NtNt);
	Inverse_LU(mmse_NtNt, Nt, Nt);
	c_matmultiple(mmse_NtNt, transB, mmse_Nt, transB, Nt, Nt, Nt, 1, x_hat);

	return 0;
}

void MMSE_free()
{
	free(mmse_Nt);
	free(mmse_NtNt);
	//free(sigma2eye);
}