set SynModuleInfo {
  {SRCNAME map.4 MODELNAME map_4 RTLNAME MHGD_detect_accel_map_4
    SUBMODULES {
      {MODELNAME MHGD_detect_accel_sdiv_96ns_1ns_96_100_1 RTLNAME MHGD_detect_accel_sdiv_96ns_1ns_96_100_1 BINDTYPE op TYPE sdiv IMPL auto LATENCY 99 ALLOW_PRAGMA 1}
    }
  }
  {SRCNAME MHGD_detect_accel_Pipeline_VITIS_LOOP_489_3 MODELNAME MHGD_detect_accel_Pipeline_VITIS_LOOP_489_3 RTLNAME MHGD_detect_accel_MHGD_detect_accel_Pipeline_VITIS_LOOP_489_3
    SUBMODULES {
      {MODELNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init RTLNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init BINDTYPE interface TYPE internal_upc_flow_control INSTNAME MHGD_detect_accel_flow_control_loop_pipe_sequential_init_U}
    }
  }
  {SRCNAME MHGD_detect_accel MODELNAME MHGD_detect_accel RTLNAME MHGD_detect_accel IS_TOP 1
    SUBMODULES {
      {MODELNAME MHGD_detect_accel_gmem_m_axi RTLNAME MHGD_detect_accel_gmem_m_axi BINDTYPE interface TYPE adapter IMPL m_axi}
      {MODELNAME MHGD_detect_accel_control_s_axi RTLNAME MHGD_detect_accel_control_s_axi BINDTYPE interface TYPE interface_s_axilite}
    }
  }
}
