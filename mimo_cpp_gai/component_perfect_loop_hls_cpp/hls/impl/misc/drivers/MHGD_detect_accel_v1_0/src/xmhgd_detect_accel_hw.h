// ==============================================================
// Vitis HLS - High-Level Synthesis from C, C++ and OpenCL v2024.2 (64-bit)
// Tool Version Limit: 2024.11
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// 
// ==============================================================
// control
// 0x00 : reserved
// 0x04 : reserved
// 0x08 : reserved
// 0x0c : reserved
// 0x10 : Data signal of x_hat
//        bit 31~0 - x_hat[31:0] (Read/Write)
// 0x14 : Data signal of x_hat
//        bit 31~0 - x_hat[63:32] (Read/Write)
// 0x18 : reserved
// 0x1c : Data signal of H
//        bit 31~0 - H[31:0] (Read/Write)
// 0x20 : Data signal of H
//        bit 31~0 - H[63:32] (Read/Write)
// 0x24 : reserved
// 0x28 : Data signal of y
//        bit 31~0 - y[31:0] (Read/Write)
// 0x2c : Data signal of y
//        bit 31~0 - y[63:32] (Read/Write)
// 0x30 : reserved
// 0x34 : Data signal of v_tb
//        bit 31~0 - v_tb[31:0] (Read/Write)
// 0x38 : Data signal of v_tb
//        bit 31~0 - v_tb[63:32] (Read/Write)
// 0x3c : reserved
// (SC = Self Clear, COR = Clear on Read, TOW = Toggle on Write, COH = Clear on Handshake)

#define XMHGD_DETECT_ACCEL_CONTROL_ADDR_X_HAT_DATA 0x10
#define XMHGD_DETECT_ACCEL_CONTROL_BITS_X_HAT_DATA 64
#define XMHGD_DETECT_ACCEL_CONTROL_ADDR_H_DATA     0x1c
#define XMHGD_DETECT_ACCEL_CONTROL_BITS_H_DATA     64
#define XMHGD_DETECT_ACCEL_CONTROL_ADDR_Y_DATA     0x28
#define XMHGD_DETECT_ACCEL_CONTROL_BITS_Y_DATA     64
#define XMHGD_DETECT_ACCEL_CONTROL_ADDR_V_TB_DATA  0x34
#define XMHGD_DETECT_ACCEL_CONTROL_BITS_V_TB_DATA  64

