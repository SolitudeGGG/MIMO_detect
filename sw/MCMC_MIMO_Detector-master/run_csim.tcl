 open_component -reset MHGD_sw -flow_target vivado

set_top MIMO_detect

# Add design files
add_files MHGD_accel_sw.cpp
add_files MHGD_accel_sw.h
add_files MyComplex.h


# Add test bench & files
add_files -tb main_sw.cpp
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
add_files -tb /home/ggg_wufuqi/hls/MHGD/gaussian_random_values.txt -cflags -I.

set_top MHGD_detect_accel_sw

set_part  {xck26-sfvc784-2LV-c}
create_clock -period 10

csim_design
exit