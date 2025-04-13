#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "util_1.h"
#include "math.h"
#include "time.h"
#include "MyComplex_1.h"

void MMSE_Init(int Nt, int Nr);
float MMSE_detect(int Nt, int Nr, MyComplex* y, MyComplex* H, MyComplex x_hat[8], float sigma2);
void MMSE_free();