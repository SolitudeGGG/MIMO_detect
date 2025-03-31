#include "hls_design_meta.h"
const Port_Property HLS_Design_Meta::port_props[]={
	Port_Property("ap_start", 1, hls_in, -1, "", "", 1),
	Port_Property("ap_done", 1, hls_out, -1, "", "", 1),
	Port_Property("ap_idle", 1, hls_out, -1, "", "", 1),
	Port_Property("ap_ready", 1, hls_out, -1, "", "", 1),
	Port_Property("x_hat", 64, hls_in, 0, "ap_none", "in_data", 1),
	Port_Property("Nt", 32, hls_in, 1, "ap_none", "in_data", 1),
	Port_Property("Nr", 32, hls_in, 2, "ap_none", "in_data", 1),
	Port_Property("mu", 32, hls_in, 3, "ap_none", "in_data", 1),
	Port_Property("H", 64, hls_in, 4, "ap_none", "in_data", 1),
	Port_Property("y", 64, hls_in, 5, "ap_none", "in_data", 1),
	Port_Property("sigma2", 32, hls_in, 6, "ap_none", "in_data", 1),
	Port_Property("mmse_init", 32, hls_in, 7, "ap_none", "in_data", 1),
	Port_Property("lr_approx", 32, hls_in, 8, "ap_none", "in_data", 1),
	Port_Property("iter", 32, hls_in, 9, "ap_none", "in_data", 1),
	Port_Property("v_tb", 64, hls_in, 10, "ap_none", "in_data", 1),
	Port_Property("ap_return", 32, hls_out, -1, "", "", 1),
	Port_Property("ap_rst", 1, hls_in, -1, "", "", 1),
};
const char* HLS_Design_Meta::dut_name = "MHGD_detect_accel";
