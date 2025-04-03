// ==============================================================
// Vitis HLS - High-Level Synthesis from C, C++ and OpenCL v2024.2 (64-bit)
// Tool Version Limit: 2024.11
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// 
// ==============================================================
/***************************** Include Files *********************************/
#include "xmhgd_detect_accel.h"

/************************** Function Implementation *************************/
#ifndef __linux__
int XMhgd_detect_accel_CfgInitialize(XMhgd_detect_accel *InstancePtr, XMhgd_detect_accel_Config *ConfigPtr) {
    Xil_AssertNonvoid(InstancePtr != NULL);
    Xil_AssertNonvoid(ConfigPtr != NULL);

    InstancePtr->Control_BaseAddress = ConfigPtr->Control_BaseAddress;
    InstancePtr->IsReady = XIL_COMPONENT_IS_READY;

    return XST_SUCCESS;
}
#endif

void XMhgd_detect_accel_Set_x_hat(XMhgd_detect_accel *InstancePtr, u64 Data) {
    Xil_AssertVoid(InstancePtr != NULL);
    Xil_AssertVoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_X_HAT_DATA, (u32)(Data));
    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_X_HAT_DATA + 4, (u32)(Data >> 32));
}

u64 XMhgd_detect_accel_Get_x_hat(XMhgd_detect_accel *InstancePtr) {
    u64 Data;

    Xil_AssertNonvoid(InstancePtr != NULL);
    Xil_AssertNonvoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    Data = XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_X_HAT_DATA);
    Data += (u64)XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_X_HAT_DATA + 4) << 32;
    return Data;
}

void XMhgd_detect_accel_Set_H(XMhgd_detect_accel *InstancePtr, u64 Data) {
    Xil_AssertVoid(InstancePtr != NULL);
    Xil_AssertVoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_H_DATA, (u32)(Data));
    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_H_DATA + 4, (u32)(Data >> 32));
}

u64 XMhgd_detect_accel_Get_H(XMhgd_detect_accel *InstancePtr) {
    u64 Data;

    Xil_AssertNonvoid(InstancePtr != NULL);
    Xil_AssertNonvoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    Data = XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_H_DATA);
    Data += (u64)XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_H_DATA + 4) << 32;
    return Data;
}

void XMhgd_detect_accel_Set_y(XMhgd_detect_accel *InstancePtr, u64 Data) {
    Xil_AssertVoid(InstancePtr != NULL);
    Xil_AssertVoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_Y_DATA, (u32)(Data));
    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_Y_DATA + 4, (u32)(Data >> 32));
}

u64 XMhgd_detect_accel_Get_y(XMhgd_detect_accel *InstancePtr) {
    u64 Data;

    Xil_AssertNonvoid(InstancePtr != NULL);
    Xil_AssertNonvoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    Data = XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_Y_DATA);
    Data += (u64)XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_Y_DATA + 4) << 32;
    return Data;
}

void XMhgd_detect_accel_Set_v_tb(XMhgd_detect_accel *InstancePtr, u64 Data) {
    Xil_AssertVoid(InstancePtr != NULL);
    Xil_AssertVoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_V_TB_DATA, (u32)(Data));
    XMhgd_detect_accel_WriteReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_V_TB_DATA + 4, (u32)(Data >> 32));
}

u64 XMhgd_detect_accel_Get_v_tb(XMhgd_detect_accel *InstancePtr) {
    u64 Data;

    Xil_AssertNonvoid(InstancePtr != NULL);
    Xil_AssertNonvoid(InstancePtr->IsReady == XIL_COMPONENT_IS_READY);

    Data = XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_V_TB_DATA);
    Data += (u64)XMhgd_detect_accel_ReadReg(InstancePtr->Control_BaseAddress, XMHGD_DETECT_ACCEL_CONTROL_ADDR_V_TB_DATA + 4) << 32;
    return Data;
}

