#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "math.h"
#include "time.h"
#include "MyComplex.h"

void MMSE_Init(int Nt, int Nr);
float MMSE_detect(int Nt, int Nr, MyComplex* y, MyComplex* H, MyComplex* x_hat, float sigma2);
void MMSE_free();