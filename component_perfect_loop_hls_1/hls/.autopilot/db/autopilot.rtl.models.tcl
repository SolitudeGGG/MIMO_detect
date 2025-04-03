set SynModuleInfo {
  {SRCNAME MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3 MODELNAME MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3 RTLNAME MHGD_detect_accel_MHGD_detect_accel_Pipeline_VITIS_LOOP_483_3
    SUBMODULES {
      {MODELNAME MHGD_detect_accel_fcmp_32ns_32ns_1_2_no_dsp_1 RTLNAME MHGD_detect_accel_fcmp_32ns_32ns_1_2_no_dsp_1 BINDTYPE op TYPE fcmp IMPL auto LATENCY 1 ALLOW_PRAGMA 1}
      {MODELNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init RTLNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init BINDTYPE interface TYPE internal_upc_flow_control INSTNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init_U}
    }
  }
  {SRCNAME MHGD_detect_accel MODELNAME MHGD_detect_accel RTLNAME MHGD_detect_accel IS_TOP 1
    SUBMODULES {
      {MODELNAME MHGD_detect_accel_fptrunc_64ns_32_2_no_dsp_1 RTLNAME MHGD_detect_accel_fptrunc_64ns_32_2_no_dsp_1 BINDTYPE op TYPE fptrunc IMPL auto LATENCY 1 ALLOW_PRAGMA 1}
      {MODELNAME MHGD_detect_accel_dexp_64ns_64ns_64_13_full_dsp_1 RTLNAME MHGD_detect_accel_dexp_64ns_64ns_64_13_full_dsp_1 BINDTYPE op TYPE dexp IMPL fulldsp LATENCY 12 ALLOW_PRAGMA 1}
    }
  }
}
