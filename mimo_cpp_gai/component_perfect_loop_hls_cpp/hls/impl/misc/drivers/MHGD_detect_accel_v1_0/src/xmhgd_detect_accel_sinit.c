// ==============================================================
// Vitis HLS - High-Level Synthesis from C, C++ and OpenCL v2024.2 (64-bit)
// Tool Version Limit: 2024.11
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// 
// ==============================================================
#ifndef __linux__

#include "xstatus.h"
#ifdef SDT
#include "xparameters.h"
#endif
#include "xmhgd_detect_accel.h"

extern XMhgd_detect_accel_Config XMhgd_detect_accel_ConfigTable[];

#ifdef SDT
XMhgd_detect_accel_Config *XMhgd_detect_accel_LookupConfig(UINTPTR BaseAddress) {
	XMhgd_detect_accel_Config *ConfigPtr = NULL;

	int Index;

	for (Index = (u32)0x0; XMhgd_detect_accel_ConfigTable[Index].Name != NULL; Index++) {
		if (!BaseAddress || XMhgd_detect_accel_ConfigTable[Index].Control_BaseAddress == BaseAddress) {
			ConfigPtr = &XMhgd_detect_accel_ConfigTable[Index];
			break;
		}
	}

	return ConfigPtr;
}

int XMhgd_detect_accel_Initialize(XMhgd_detect_accel *InstancePtr, UINTPTR BaseAddress) {
	XMhgd_detect_accel_Config *ConfigPtr;

	Xil_AssertNonvoid(InstancePtr != NULL);

	ConfigPtr = XMhgd_detect_accel_LookupConfig(BaseAddress);
	if (ConfigPtr == NULL) {
		InstancePtr->IsReady = 0;
		return (XST_DEVICE_NOT_FOUND);
	}

	return XMhgd_detect_accel_CfgInitialize(InstancePtr, ConfigPtr);
}
#else
XMhgd_detect_accel_Config *XMhgd_detect_accel_LookupConfig(u16 DeviceId) {
	XMhgd_detect_accel_Config *ConfigPtr = NULL;

	int Index;

	for (Index = 0; Index < XPAR_XMHGD_DETECT_ACCEL_NUM_INSTANCES; Index++) {
		if (XMhgd_detect_accel_ConfigTable[Index].DeviceId == DeviceId) {
			ConfigPtr = &XMhgd_detect_accel_ConfigTable[Index];
			break;
		}
	}

	return ConfigPtr;
}

int XMhgd_detect_accel_Initialize(XMhgd_detect_accel *InstancePtr, u16 DeviceId) {
	XMhgd_detect_accel_Config *ConfigPtr;

	Xil_AssertNonvoid(InstancePtr != NULL);

	ConfigPtr = XMhgd_detect_accel_LookupConfig(DeviceId);
	if (ConfigPtr == NULL) {
		InstancePtr->IsReady = 0;
		return (XST_DEVICE_NOT_FOUND);
	}

	return XMhgd_detect_accel_CfgInitialize(InstancePtr, ConfigPtr);
}
#endif

#endif

