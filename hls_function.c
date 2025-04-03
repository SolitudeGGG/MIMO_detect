#include "MyComplex.h"
#include "util.h"
#include "MHGD_accel.h"

void initialize_parameters() {
   /*----------------------------------------------------------------*/
	c_eye_generate(sigma2eye, Nt, sigma2 / (dqam * dqam));

	/*二阶梯度下降，计算grad_preconditioner(梯度更新的预条件矩阵)*/
	c_matmultiple(H, transA, H, transB, Nr, Nt, Nr, Nt, HH_H);
	my_complex_add(Nt * Nt, HH_H, sigma2eye, grad_preconditioner);
	Inverse_LU(grad_preconditioner, Nt, Nt);

	alpha = 1.0f / pow((float)Nt / 8.0f, 1.0f / 3.0f);

	MHGD_end[(int)gradpre] = clock();

	/*协方差矩阵*/
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
		/*映射到归一化星座图中 xhat = constellation_norm[np.argmin(abs(x_mmse * np.ones(nt, 2 * *mu) - constellation_norm), axis = 1)].reshape(-1, 1)*/
		map(mu, Nt, dqam, x_mmse, x_hat);
	}
	else
	{
		/*xhat = constellation_norm[np.random.randint(low=0, high=2 ** mu, size=(samplers, nt, 1))].copy()*/
		generateUniformRandoms_int(Nt, x_init, mu);
		for (i = 0; i < Nt; i++)
			x_hat[i] = constellation_norm[x_init[i]];
	}
	/*计算剩余向量r=y-Hx*/
	c_matmultiple(H, transB, x_hat, transB, Nr, Nt, Nt, 1, r);
	my_complex_sub(Nr, y, r, r);

	/*计算剩余向量的范数（就是模值）*/
	c_matmultiple(r, transA, r, transB, Nr, 1, Nr, 1, temp_1);
	r_norm = temp_1[0].real;
	my_complex_copy(Nt, x_hat, 1, x_survivor, 1);
	r_norm_survivor = r_norm;
	/*确定最优学习率*/
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
};