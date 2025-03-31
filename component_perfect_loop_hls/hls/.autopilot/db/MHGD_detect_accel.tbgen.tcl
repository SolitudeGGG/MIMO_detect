set moduleName MHGD_detect_accel
set isTopModule 1
set isCombinational 0
set isDatapathOnly 0
set isPipelined 0
set pipeline_type none
set FunctionProtocol ap_ctrl_hs
set isOneStateSeq 0
set ProfileFlag 0
set StallSigGenFlag 0
set isEnableWaveformDebug 1
set hasInterrupt 0
set DLRegFirstOffset 0
set DLRegItemOffset 0
set svuvm_can_support 1
set cdfgNum 3
set C_modelName {MHGD_detect_accel}
set C_modelType { float 32 }
set ap_memory_interface_dict [dict create]
set C_modelArgList {
	{ x_hat int 64 unused {pointer 0}  }
	{ Nt int 32 regular  }
	{ Nr int 32 unused  }
	{ mu int 32 unused  }
	{ H int 64 unused {pointer 0}  }
	{ y int 64 unused {pointer 0}  }
	{ sigma2 float 32 unused  }
	{ mmse_init int 32 regular  }
	{ lr_approx int 32 unused  }
	{ iter int 32 regular  }
	{ v_tb int 64 unused {pointer 0}  }
}
set hasAXIMCache 0
set l_AXIML2Cache [list]
set AXIMCacheInstDict [dict create]
set C_modelArgMapList {[ 
	{ "Name" : "x_hat", "interface" : "wire", "bitwidth" : 64, "direction" : "READONLY"} , 
 	{ "Name" : "Nt", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "Nr", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "mu", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "H", "interface" : "wire", "bitwidth" : 64, "direction" : "READONLY"} , 
 	{ "Name" : "y", "interface" : "wire", "bitwidth" : 64, "direction" : "READONLY"} , 
 	{ "Name" : "sigma2", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "mmse_init", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "lr_approx", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "iter", "interface" : "wire", "bitwidth" : 32, "direction" : "READONLY"} , 
 	{ "Name" : "v_tb", "interface" : "wire", "bitwidth" : 64, "direction" : "READONLY"} , 
 	{ "Name" : "ap_return", "interface" : "wire", "bitwidth" : 32} ]}
# RTL Port declarations: 
set portNum 18
set portList { 
	{ ap_clk sc_in sc_logic 1 clock -1 } 
	{ ap_rst sc_in sc_logic 1 reset -1 active_high_sync } 
	{ ap_start sc_in sc_logic 1 start -1 } 
	{ ap_done sc_out sc_logic 1 predone -1 } 
	{ ap_idle sc_out sc_logic 1 done -1 } 
	{ ap_ready sc_out sc_logic 1 ready -1 } 
	{ x_hat sc_in sc_lv 64 signal 0 } 
	{ Nt sc_in sc_lv 32 signal 1 } 
	{ Nr sc_in sc_lv 32 signal 2 } 
	{ mu sc_in sc_lv 32 signal 3 } 
	{ H sc_in sc_lv 64 signal 4 } 
	{ y sc_in sc_lv 64 signal 5 } 
	{ sigma2 sc_in sc_lv 32 signal 6 } 
	{ mmse_init sc_in sc_lv 32 signal 7 } 
	{ lr_approx sc_in sc_lv 32 signal 8 } 
	{ iter sc_in sc_lv 32 signal 9 } 
	{ v_tb sc_in sc_lv 64 signal 10 } 
	{ ap_return sc_out sc_lv 32 signal -1 } 
}
set NewPortList {[ 
	{ "name": "ap_clk", "direction": "in", "datatype": "sc_logic", "bitwidth":1, "type": "clock", "bundle":{"name": "ap_clk", "role": "default" }} , 
 	{ "name": "ap_rst", "direction": "in", "datatype": "sc_logic", "bitwidth":1, "type": "reset", "bundle":{"name": "ap_rst", "role": "default" }} , 
 	{ "name": "ap_start", "direction": "in", "datatype": "sc_logic", "bitwidth":1, "type": "start", "bundle":{"name": "ap_start", "role": "default" }} , 
 	{ "name": "ap_done", "direction": "out", "datatype": "sc_logic", "bitwidth":1, "type": "predone", "bundle":{"name": "ap_done", "role": "default" }} , 
 	{ "name": "ap_idle", "direction": "out", "datatype": "sc_logic", "bitwidth":1, "type": "done", "bundle":{"name": "ap_idle", "role": "default" }} , 
 	{ "name": "ap_ready", "direction": "out", "datatype": "sc_logic", "bitwidth":1, "type": "ready", "bundle":{"name": "ap_ready", "role": "default" }} , 
 	{ "name": "x_hat", "direction": "in", "datatype": "sc_lv", "bitwidth":64, "type": "signal", "bundle":{"name": "x_hat", "role": "default" }} , 
 	{ "name": "Nt", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "Nt", "role": "default" }} , 
 	{ "name": "Nr", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "Nr", "role": "default" }} , 
 	{ "name": "mu", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "mu", "role": "default" }} , 
 	{ "name": "H", "direction": "in", "datatype": "sc_lv", "bitwidth":64, "type": "signal", "bundle":{"name": "H", "role": "default" }} , 
 	{ "name": "y", "direction": "in", "datatype": "sc_lv", "bitwidth":64, "type": "signal", "bundle":{"name": "y", "role": "default" }} , 
 	{ "name": "sigma2", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "sigma2", "role": "default" }} , 
 	{ "name": "mmse_init", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "mmse_init", "role": "default" }} , 
 	{ "name": "lr_approx", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "lr_approx", "role": "default" }} , 
 	{ "name": "iter", "direction": "in", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "iter", "role": "default" }} , 
 	{ "name": "v_tb", "direction": "in", "datatype": "sc_lv", "bitwidth":64, "type": "signal", "bundle":{"name": "v_tb", "role": "default" }} , 
 	{ "name": "ap_return", "direction": "out", "datatype": "sc_lv", "bitwidth":32, "type": "signal", "bundle":{"name": "ap_return", "role": "default" }}  ]}

set RtlHierarchyInfo {[
	{"ID" : "0", "Level" : "0", "Path" : "`AUTOTB_DUT_INST", "Parent" : "", "Child" : ["1", "4", "5"],
		"CDFG" : "MHGD_detect_accel",
		"Protocol" : "ap_ctrl_hs",
		"ControlExist" : "1", "ap_start" : "1", "ap_ready" : "1", "ap_done" : "1", "ap_continue" : "0", "ap_idle" : "1", "real_start" : "0",
		"Pipeline" : "None", "UnalignedPipeline" : "0", "RewindPipeline" : "0", "ProcessNetwork" : "0",
		"II" : "0",
		"VariableLatency" : "1", "ExactLatency" : "-1", "EstimateLatencyMin" : "1", "EstimateLatencyMax" : "30018",
		"Combinational" : "0",
		"Datapath" : "0",
		"ClockEnable" : "0",
		"HasSubDataflow" : "0",
		"InDataflowNetwork" : "0",
		"HasNonBlockingOperation" : "0",
		"IsBlackBox" : "0",
		"Port" : [
			{"Name" : "x_hat", "Type" : "None", "Direction" : "I"},
			{"Name" : "Nt", "Type" : "None", "Direction" : "I"},
			{"Name" : "Nr", "Type" : "None", "Direction" : "I"},
			{"Name" : "mu", "Type" : "None", "Direction" : "I"},
			{"Name" : "H", "Type" : "None", "Direction" : "I"},
			{"Name" : "y", "Type" : "None", "Direction" : "I"},
			{"Name" : "sigma2", "Type" : "None", "Direction" : "I"},
			{"Name" : "mmse_init", "Type" : "None", "Direction" : "I"},
			{"Name" : "lr_approx", "Type" : "None", "Direction" : "I"},
			{"Name" : "iter", "Type" : "None", "Direction" : "I"},
			{"Name" : "v_tb", "Type" : "None", "Direction" : "I"}]},
	{"ID" : "1", "Level" : "1", "Path" : "`AUTOTB_DUT_INST.grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3_fu_82", "Parent" : "0", "Child" : ["2", "3"],
		"CDFG" : "MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3",
		"Protocol" : "ap_ctrl_hs",
		"ControlExist" : "1", "ap_start" : "1", "ap_ready" : "1", "ap_done" : "1", "ap_continue" : "0", "ap_idle" : "1", "real_start" : "0",
		"Pipeline" : "None", "UnalignedPipeline" : "0", "RewindPipeline" : "0", "ProcessNetwork" : "0",
		"II" : "0",
		"VariableLatency" : "1", "ExactLatency" : "-1", "EstimateLatencyMin" : "10003", "EstimateLatencyMax" : "30003",
		"Combinational" : "0",
		"Datapath" : "0",
		"ClockEnable" : "0",
		"HasSubDataflow" : "0",
		"InDataflowNetwork" : "0",
		"HasNonBlockingOperation" : "0",
		"IsBlackBox" : "0",
		"Port" : [
			{"Name" : "iter", "Type" : "None", "Direction" : "I"},
			{"Name" : "bitcast_ln524", "Type" : "None", "Direction" : "I"},
			{"Name" : "empty", "Type" : "None", "Direction" : "I"},
			{"Name" : "p_acc", "Type" : "None", "Direction" : "I"},
			{"Name" : "icmp_ln449", "Type" : "None", "Direction" : "I"}],
		"Loop" : [
			{"Name" : "VITIS_LOOP_483_3", "PipelineType" : "UPC",
				"LoopDec" : {"FSMBitwidth" : "2", "FirstState" : "ap_ST_fsm_state1", "FirstStateIter" : "", "FirstStateBlock" : "ap_ST_fsm_state1_blk", "LastState" : "ap_ST_fsm_state2", "LastStateIter" : "", "LastStateBlock" : "ap_ST_fsm_state2_blk", "QuitState" : "ap_ST_fsm_state2", "QuitStateIter" : "", "QuitStateBlock" : "ap_ST_fsm_state2_blk", "OneDepthLoop" : "1", "has_ap_ctrl" : "1", "has_continue" : "0"}}]},
	{"ID" : "2", "Level" : "2", "Path" : "`AUTOTB_DUT_INST.grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3_fu_82.fcmp_32ns_32ns_1_2_no_dsp_1_U1", "Parent" : "1"},
	{"ID" : "3", "Level" : "2", "Path" : "`AUTOTB_DUT_INST.grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3_fu_82.flow_control_loop_pipe_sequential_init_U", "Parent" : "1"},
	{"ID" : "4", "Level" : "1", "Path" : "`AUTOTB_DUT_INST.fptrunc_64ns_32_2_no_dsp_1_U8", "Parent" : "0"},
	{"ID" : "5", "Level" : "1", "Path" : "`AUTOTB_DUT_INST.dexp_64ns_64ns_64_13_full_dsp_1_U9", "Parent" : "0"}]}


set ArgLastReadFirstWriteLatency {
	MHGD_detect_accel {
		x_hat {Type I LastRead -1 FirstWrite -1}
		Nt {Type I LastRead 0 FirstWrite -1}
		Nr {Type I LastRead -1 FirstWrite -1}
		mu {Type I LastRead -1 FirstWrite -1}
		H {Type I LastRead -1 FirstWrite -1}
		y {Type I LastRead -1 FirstWrite -1}
		sigma2 {Type I LastRead -1 FirstWrite -1}
		mmse_init {Type I LastRead 0 FirstWrite -1}
		lr_approx {Type I LastRead -1 FirstWrite -1}
		iter {Type I LastRead 0 FirstWrite -1}
		v_tb {Type I LastRead -1 FirstWrite -1}}
	MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3 {
		iter {Type I LastRead 0 FirstWrite -1}
		bitcast_ln524 {Type I LastRead 0 FirstWrite -1}
		empty {Type I LastRead 0 FirstWrite -1}
		p_acc {Type I LastRead 0 FirstWrite -1}
		icmp_ln449 {Type I LastRead 0 FirstWrite -1}}}

set hasDtUnsupportedChannel 0

set PerformanceInfo {[
	{"Name" : "Latency", "Min" : "1", "Max" : "30018"}
	, {"Name" : "Interval", "Min" : "2", "Max" : "30019"}
]}

set PipelineEnableSignalInfo {[
]}

set Spec2ImplPortList { 
	x_hat { ap_none {  { x_hat in_data 0 64 } } }
	Nt { ap_none {  { Nt in_data 0 32 } } }
	Nr { ap_none {  { Nr in_data 0 32 } } }
	mu { ap_none {  { mu in_data 0 32 } } }
	H { ap_none {  { H in_data 0 64 } } }
	y { ap_none {  { y in_data 0 64 } } }
	sigma2 { ap_none {  { sigma2 in_data 0 32 } } }
	mmse_init { ap_none {  { mmse_init in_data 0 32 } } }
	lr_approx { ap_none {  { lr_approx in_data 0 32 } } }
	iter { ap_none {  { iter in_data 0 32 } } }
	v_tb { ap_none {  { v_tb in_data 0 64 } } }
}

set maxi_interface_dict [dict create]

# RTL port scheduling information:
set fifoSchedulingInfoList { 
}

# RTL bus port read request latency information:
set busReadReqLatencyList { 
}

# RTL bus port write response latency information:
set busWriteResLatencyList { 
}

# RTL array port load latency information:
set memoryLoadLatencyList { 
}
