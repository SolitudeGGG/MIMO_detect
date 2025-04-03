; ModuleID = '/home/ggg_wufuqi/hls/MIMO_detect-main/mimo_cpp_gai/component_perfect_loop_hls_cpp/hls/.autopilot/db/a.g.ld.5.gdce.bc'
source_filename = "llvm-link"
target datalayout = "e-m:e-i64:64-i128:128-i256:256-i512:512-i1024:1024-i2048:2048-i4096:4096-n8:16:32:64-S128-v16:16-v24:32-v32:32-v48:64-v96:128-v192:256-v256:256-v512:512-v1024:1024"
target triple = "fpga64-xilinx-none"

%struct.MyComplex = type { %"struct.ap_fixed<64, 32, AP_TRN, AP_WRAP, 0>", %"struct.ap_fixed<64, 32, AP_TRN, AP_WRAP, 0>" }
%"struct.ap_fixed<64, 32, AP_TRN, AP_WRAP, 0>" = type { %"struct.ap_fixed_base<64, 32, true, AP_TRN, AP_WRAP, 0>" }
%"struct.ap_fixed_base<64, 32, true, AP_TRN, AP_WRAP, 0>" = type { %"struct.ssdm_int<64, true>" }
%"struct.ssdm_int<64, true>" = type { i64 }

; Function Attrs: noinline willreturn
define float @apatb_MHGD_detect_accel_ir(%struct.MyComplex* noalias nonnull "maxi" %x_hat, i32 %Nt, i32 %Nr, i32 %mu, %struct.MyComplex* noalias nocapture nonnull readonly "maxi" %H, %struct.MyComplex* noalias nocapture nonnull readonly "maxi" %y, float %sigma2, i32 %mmse_init, i32 %lr_approx, i32 %iter, %struct.MyComplex* noalias nocapture nonnull readonly "maxi" %v_tb) local_unnamed_addr #0 {
entry:
  %x_hat_copy = alloca %struct.MyComplex, align 512
  %H_copy = alloca %struct.MyComplex, align 512
  %y_copy = alloca %struct.MyComplex, align 512
  %v_tb_copy = alloca %struct.MyComplex, align 512
  call fastcc void @copy_in(%struct.MyComplex* nonnull %x_hat, %struct.MyComplex* nonnull align 512 %x_hat_copy, %struct.MyComplex* nonnull %H, %struct.MyComplex* nonnull align 512 %H_copy, %struct.MyComplex* nonnull %y, %struct.MyComplex* nonnull align 512 %y_copy, %struct.MyComplex* nonnull %v_tb, %struct.MyComplex* nonnull align 512 %v_tb_copy)
  %0 = call float @apatb_MHGD_detect_accel_hw(%struct.MyComplex* %x_hat_copy, i32 %Nt, i32 %Nr, i32 %mu, %struct.MyComplex* %H_copy, %struct.MyComplex* %y_copy, float %sigma2, i32 %mmse_init, i32 %lr_approx, i32 %iter, %struct.MyComplex* %v_tb_copy)
  call void @copy_back(%struct.MyComplex* %x_hat, %struct.MyComplex* %x_hat_copy, %struct.MyComplex* %H, %struct.MyComplex* %H_copy, %struct.MyComplex* %y, %struct.MyComplex* %y_copy, %struct.MyComplex* %v_tb, %struct.MyComplex* %v_tb_copy)
  ret float %0
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_in(%struct.MyComplex* noalias readonly, %struct.MyComplex* noalias align 512, %struct.MyComplex* noalias readonly, %struct.MyComplex* noalias align 512, %struct.MyComplex* noalias readonly, %struct.MyComplex* noalias align 512, %struct.MyComplex* noalias readonly, %struct.MyComplex* noalias align 512) unnamed_addr #1 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* align 512 %1, %struct.MyComplex* %0)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* align 512 %3, %struct.MyComplex* %2)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* align 512 %5, %struct.MyComplex* %4)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* align 512 %7, %struct.MyComplex* %6)
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* noalias align 512 %dst, %struct.MyComplex* noalias readonly %src) unnamed_addr #2 {
entry:
  %0 = icmp eq %struct.MyComplex* %dst, null
  %1 = icmp eq %struct.MyComplex* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %src.0.0.0.05 = getelementptr %struct.MyComplex, %struct.MyComplex* %src, i64 0, i32 0, i32 0, i32 0, i32 0
  %dst.0.0.0.06 = getelementptr %struct.MyComplex, %struct.MyComplex* %dst, i64 0, i32 0, i32 0, i32 0, i32 0
  %3 = load i64, i64* %src.0.0.0.05, align 8
  store i64 %3, i64* %dst.0.0.0.06, align 512
  %src.1.0.0.011 = getelementptr %struct.MyComplex, %struct.MyComplex* %src, i64 0, i32 1, i32 0, i32 0, i32 0
  %dst.1.0.0.012 = getelementptr %struct.MyComplex, %struct.MyComplex* %dst, i64 0, i32 1, i32 0, i32 0, i32 0
  %4 = load i64, i64* %src.1.0.0.011, align 8
  store i64 %4, i64* %dst.1.0.0.012, align 8
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_out(%struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512) unnamed_addr #3 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %0, %struct.MyComplex* align 512 %1)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %2, %struct.MyComplex* align 512 %3)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %4, %struct.MyComplex* align 512 %5)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %6, %struct.MyComplex* align 512 %7)
  ret void
}

declare float @apatb_MHGD_detect_accel_hw(%struct.MyComplex*, i32, i32, i32, %struct.MyComplex*, %struct.MyComplex*, float, i32, i32, i32, %struct.MyComplex*)

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_back(%struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512, %struct.MyComplex* noalias, %struct.MyComplex* noalias readonly align 512) unnamed_addr #3 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %0, %struct.MyComplex* align 512 %1)
  ret void
}

declare float @MHGD_detect_accel_hw_stub(%struct.MyComplex* noalias nonnull, i32, i32, i32, %struct.MyComplex* noalias nocapture nonnull readonly, %struct.MyComplex* noalias nocapture nonnull readonly, float, i32, i32, i32, %struct.MyComplex* noalias nocapture nonnull readonly)

define float @MHGD_detect_accel_hw_stub_wrapper(%struct.MyComplex*, i32, i32, i32, %struct.MyComplex*, %struct.MyComplex*, float, i32, i32, i32, %struct.MyComplex*) #4 {
entry:
  call void @copy_out(%struct.MyComplex* null, %struct.MyComplex* %0, %struct.MyComplex* null, %struct.MyComplex* %4, %struct.MyComplex* null, %struct.MyComplex* %5, %struct.MyComplex* null, %struct.MyComplex* %10)
  %11 = call float @MHGD_detect_accel_hw_stub(%struct.MyComplex* %0, i32 %1, i32 %2, i32 %3, %struct.MyComplex* %4, %struct.MyComplex* %5, float %6, i32 %7, i32 %8, i32 %9, %struct.MyComplex* %10)
  call void @copy_in(%struct.MyComplex* null, %struct.MyComplex* %0, %struct.MyComplex* null, %struct.MyComplex* %4, %struct.MyComplex* null, %struct.MyComplex* %5, %struct.MyComplex* null, %struct.MyComplex* %10)
  ret float %11
}

attributes #0 = { noinline willreturn "fpga.wrapper.func"="wrapper" }
attributes #1 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyin" }
attributes #2 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="onebyonecpy_hls" }
attributes #3 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyout" }
attributes #4 = { "fpga.wrapper.func"="stub" }

!llvm.dbg.cu = !{}
!llvm.ident = !{!0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0}
!llvm.module.flags = !{!1, !2, !3}
!blackbox_cfg = !{!4}

!0 = !{!"clang version 7.0.0 "}
!1 = !{i32 2, !"Dwarf Version", i32 4}
!2 = !{i32 2, !"Debug Info Version", i32 3}
!3 = !{i32 1, !"wchar_size", i32 4}
!4 = !{}
