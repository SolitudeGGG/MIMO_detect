#pragma once

typedef ap_fixed<64, 32>Myreal;
typedef ap_fixed<64, 32>Myimage;
typedef ap_fixed<64, 32>like_float;

// 自定义复数结构体
typedef struct {
    Myreal real;  // 实部
    Myimage imag;  // 虚部
} MyComplex;


