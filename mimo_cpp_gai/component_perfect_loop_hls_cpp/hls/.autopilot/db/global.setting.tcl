
set TopModule "MHGD_detect_accel"
set ClockPeriod 10
set ClockList ap_clk
set AxiliteClockList {}
set HasVivadoClockPeriod 0
set CombLogicFlag 0
set PipelineFlag 0
set DataflowTaskPipelineFlag 1
set TrivialPipelineFlag 0
set noPortSwitchingFlag 0
set FloatingPointFlag 0
set FftOrFirFlag 0
set NbRWValue 0
set intNbAccess 0
set NewDSPMapping 1
set HasDSPModule 0
set ResetLevelFlag 0
set ResetStyle control
set ResetSyncFlag 1
set ResetRegisterFlag 0
set ResetVariableFlag 0
set ResetRegisterNum 0
set FsmEncStyle onehot
set MaxFanout 0
set RtlPrefix {}
set RtlSubPrefix MHGD_detect_accel_
set ExtraCCFlags {}
set ExtraCLdFlags {}
set SynCheckOptions {}
set PresynOptions {}
set PreprocOptions {}
set SchedOptions {}
set BindOptions {}
set RtlGenOptions {}
set RtlWriterOptions {}
set CbcGenFlag {}
set CasGenFlag {}
set CasMonitorFlag {}
set AutoSimOptions {}
set ExportMCPathFlag 0
set SCTraceFileName mytrace
set SCTraceFileFormat vcd
set SCTraceOption all
set TargetInfo xck26:-sfvc784:-2LV-c
set SourceFiles {sc {} c {../../MMSE_1.cpp ../../util_1.cpp ../../MIMO_simulation_1.cpp ../../MHGD_accel_1.cpp}}
set SourceFlags {sc {} c {{} {} {} {}}}
set DirectiveFile {}
set TBFiles {verilog {/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt ../../tb_1.cpp} bc {/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt ../../tb_1.cpp} sc {/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt ../../tb_1.cpp} vhdl {/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt ../../tb_1.cpp} c {} cas {/home/ggg_wufuqi/hls/MHGD/_16QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/_64QAM_Constellation.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/y_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/input_file/H_SNR=5.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=25.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=20.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=15.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=10.txt /home/ggg_wufuqi/hls/MHGD/MHGD/output_file/bits_output_SNR=5.txt ../../tb_1.cpp}}
set SpecLanguage C
set TVInFiles {bc {} c {} sc {} cas {} vhdl {} verilog {}}
set TVOutFiles {bc {} c {} sc {} cas {} vhdl {} verilog {}}
set TBTops {verilog {} bc {} sc {} vhdl {} c {} cas {}}
set TBInstNames {verilog {} bc {} sc {} vhdl {} c {} cas {}}
set XDCFiles {}
set ExtraGlobalOptions {"area_timing" 1 "clock_gate" 1 "impl_flow" map "power_gate" 0}
set TBTVFileNotFound {}
set AppFile {}
set ApsFile hls.aps
set AvePath ../../.
set DefaultPlatform DefaultPlatform
set multiClockList {}
set SCPortClockMap {}
set intNbAccess 0
set PlatformFiles {{DefaultPlatform {xilinx/zynquplus/zynquplus}}}
set HPFPO 0
