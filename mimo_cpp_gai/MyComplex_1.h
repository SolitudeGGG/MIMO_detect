#pragma once
#include "ap_fixed.h"	//ap_fixed<18,6,AP_TRN_ZERO,AP_SAT>		<W,I,Q,O,N>
#include "ap_int.h"	//ap_int<N> or ap_uint<N>, 1<=N<=1024
#include "hls_math.h"	//data_t s = hls::sinf(angle);
#include "hls_stream.h"

typedef ap_fixed<64, 32>Myreal;
typedef ap_fixed<64, 32>Myimage;
typedef ap_fixed<64, 32>like_float;

// 自定义复数结构体
typedef struct {
    Myreal real;  // 实部
    Myimage imag;  // 虚部
} MyComplex;


