// ==============================================================
// Vitis HLS - High-Level Synthesis from C, C++ and OpenCL v2024.2 (64-bit)
// Tool Version Limit: 2024.11
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// 
// ==============================================================
#ifndef XMHGD_DETECT_ACCEL_H
#define XMHGD_DETECT_ACCEL_H

#ifdef __cplusplus
extern "C" {
#endif

/***************************** Include Files *********************************/
#ifndef __linux__
#include "xil_types.h"
#include "xil_assert.h"
#include "xstatus.h"
#include "xil_io.h"
#else
#include <stdint.h>
#include <assert.h>
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stddef.h>
#endif
#include "xmhgd_detect_accel_hw.h"

/**************************** Type Definitions ******************************/
#ifdef __linux__
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
#else
typedef struct {
#ifdef SDT
    char *Name;
#else
    u16 DeviceId;
#endif
    u64 Control_BaseAddress;
} XMhgd_detect_accel_Config;
#endif

typedef struct {
    u64 Control_BaseAddress;
    u32 IsReady;
} XMhgd_detect_accel;

typedef u32 word_type;

/***************** Macros (Inline Functions) Definitions *********************/
#ifndef __linux__
#define XMhgd_detect_accel_WriteReg(BaseAddress, RegOffset, Data) \
    Xil_Out32((BaseAddress) + (RegOffset), (u32)(Data))
#define XMhgd_detect_accel_ReadReg(BaseAddress, RegOffset) \
    Xil_In32((BaseAddress) + (RegOffset))
#else
#define XMhgd_detect_accel_WriteReg(BaseAddress, RegOffset, Data) \
    *(volatile u32*)((BaseAddress) + (RegOffset)) = (u32)(Data)
#define XMhgd_detect_accel_ReadReg(BaseAddress, RegOffset) \
    *(volatile u32*)((BaseAddress) + (RegOffset))

#define Xil_AssertVoid(expr)    assert(expr)
#define Xil_AssertNonvoid(expr) assert(expr)

#define XST_SUCCESS             0
#define XST_DEVICE_NOT_FOUND    2
#define XST_OPEN_DEVICE_FAILED  3
#define XIL_COMPONENT_IS_READY  1
#endif

/************************** Function Prototypes *****************************/
#ifndef __linux__
#ifdef SDT
int XMhgd_detect_accel_Initialize(XMhgd_detect_accel *InstancePtr, UINTPTR BaseAddress);
XMhgd_detect_accel_Config* XMhgd_detect_accel_LookupConfig(UINTPTR BaseAddress);
#else
int XMhgd_detect_accel_Initialize(XMhgd_detect_accel *InstancePtr, u16 DeviceId);
XMhgd_detect_accel_Config* XMhgd_detect_accel_LookupConfig(u16 DeviceId);
#endif
int XMhgd_detect_accel_CfgInitialize(XMhgd_detect_accel *InstancePtr, XMhgd_detect_accel_Config *ConfigPtr);
#else
int XMhgd_detect_accel_Initialize(XMhgd_detect_accel *InstancePtr, const char* InstanceName);
int XMhgd_detect_accel_Release(XMhgd_detect_accel *InstancePtr);
#endif


void XMhgd_detect_accel_Set_x_hat(XMhgd_detect_accel *InstancePtr, u64 Data);
u64 XMhgd_detect_accel_Get_x_hat(XMhgd_detect_accel *InstancePtr);
void XMhgd_detect_accel_Set_H(XMhgd_detect_accel *InstancePtr, u64 Data);
u64 XMhgd_detect_accel_Get_H(XMhgd_detect_accel *InstancePtr);
void XMhgd_detect_accel_Set_y(XMhgd_detect_accel *InstancePtr, u64 Data);
u64 XMhgd_detect_accel_Get_y(XMhgd_detect_accel *InstancePtr);
void XMhgd_detect_accel_Set_v_tb(XMhgd_detect_accel *InstancePtr, u64 Data);
u64 XMhgd_detect_accel_Get_v_tb(XMhgd_detect_accel *InstancePtr);

#ifdef __cplusplus
}
#endif

#endif
