#
# Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
# Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Create a project

open_component -reset component_perfect_loop_hls_cpp_cosim -flow_target vivado

# Add design files
add_files MHGD_accel_1.cpp
add_files MHGD_accel_1.h
add_files MIMO_simulation_1.cpp
add_files MIMO_simulation_1.h
# add_files MyComplex.cpp
add_files MyComplex_1.h
add_files sys_config_1.h
add_files util_1.cpp
add_files util_1.h
add_files MMSE_1.h
add_files MMSE_1.cpp

# Add test bench & files
add_files -tb tb_1.cpp
# add_files -tb result.golden.dat
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt -cflags -I.

add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt -cflags -I.

add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt -cflags -I.

add_files -tb /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt -cflags -I.
add_files -tb /home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt -cflags -I.

# Set the top-level function
set_top MHGD_detect_accel

# ########################################################
# Create a solution
# Define technology and clock rate
set_part  {xck26-sfvc784-2LV-c}
create_clock -period 10

# Set variable to select which steps to execute
set hls_exec 2


config_interface -m_axi_max_widen_bitwidth 512
config_interface -m_axi_alignment_byte_size 64

 csim_design
# Set any optimization directives
set_directive_pipeline loop_perfect/LOOP_J
# End of directives
if {$hls_exec == 1} {
	# Run Synthesis and Exit
	csynth_design
	
} elseif {$hls_exec == 2} {
	# Run Synthesis, RTL Simulation and Exit
	csynth_design
	
	cosim_design
} elseif {$hls_exec == 3} { 
	# Run Synthesis, RTL Simulation, RTL implementation and Exit
	csynth_design
	
	cosim_design
	export_design
} else {
	# Default is to exit after setup
	csynth_design
}

exit