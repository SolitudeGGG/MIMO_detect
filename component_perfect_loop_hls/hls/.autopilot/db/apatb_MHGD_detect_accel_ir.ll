; ModuleID = '/home/ggg_wufuqi/hls/MHGD/component_perfect_loop_hls/hls/.autopilot/db/a.g.ld.5.gdce.bc'
source_filename = "llvm-link"
target datalayout = "e-m:e-i64:64-i128:128-i256:256-i512:512-i1024:1024-i2048:2048-i4096:4096-n8:16:32:64-S128-v16:16-v24:32-v32:32-v48:64-v96:128-v192:256-v256:256-v512:512-v1024:1024"
target triple = "fpga64-xilinx-none"

%struct.MyComplex = type { float, float }

; Function Attrs: noinline willreturn
define float @apatb_MHGD_detect_accel_ir(%struct.MyComplex* noalias nonnull %x_hat, i32 %Nt, i32 %Nr, i32 %mu, %struct.MyComplex* noalias nocapture nonnull readonly %H, %struct.MyComplex* noalias nocapture nonnull readonly %y, float %sigma2, i32 %mmse_init, i32 %lr_approx, i32 %iter, %struct.MyComplex* noalias nocapture nonnull readonly %v_tb) local_unnamed_addr #0 {
entry:
  %x_hat_copy = alloca i64, align 512
  %H_copy = alloca i64, align 512
  %y_copy = alloca i64, align 512
  %v_tb_copy = alloca i64, align 512
  call fastcc void @copy_in(%struct.MyComplex* nonnull %x_hat, i64* nonnull align 512 %x_hat_copy, %struct.MyComplex* nonnull %H, i64* nonnull align 512 %H_copy, %struct.MyComplex* nonnull %y, i64* nonnull align 512 %y_copy, %struct.MyComplex* nonnull %v_tb, i64* nonnull align 512 %v_tb_copy)
  %0 = call float @apatb_MHGD_detect_accel_hw(i64* %x_hat_copy, i32 %Nt, i32 %Nr, i32 %mu, i64* %H_copy, i64* %y_copy, float %sigma2, i32 %mmse_init, i32 %lr_approx, i32 %iter, i64* %v_tb_copy)
  call void @copy_back(%struct.MyComplex* %x_hat, i64* %x_hat_copy, %struct.MyComplex* %H, i64* %H_copy, %struct.MyComplex* %y, i64* %y_copy, %struct.MyComplex* %v_tb, i64* %v_tb_copy)
  ret float %0
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_in(%struct.MyComplex* noalias readonly, i64* noalias align 512, %struct.MyComplex* noalias readonly, i64* noalias align 512, %struct.MyComplex* noalias readonly, i64* noalias align 512, %struct.MyComplex* noalias readonly, i64* noalias align 512) unnamed_addr #1 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex.333(i64* align 512 %1, %struct.MyComplex* %0)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex.333(i64* align 512 %3, %struct.MyComplex* %2)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex.333(i64* align 512 %5, %struct.MyComplex* %4)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex.333(i64* align 512 %7, %struct.MyComplex* %6)
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_out(%struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512) unnamed_addr #2 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %0, i64* align 512 %1)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %2, i64* align 512 %3)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %4, i64* align 512 %5)
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %6, i64* align 512 %7)
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* noalias %dst, i64* noalias readonly align 512 %src) unnamed_addr #3 {
entry:
  %0 = icmp eq %struct.MyComplex* %dst, null
  %1 = icmp eq i64* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %dst.0 = getelementptr %struct.MyComplex, %struct.MyComplex* %dst, i64 0, i32 0
  %3 = load i64, i64* %src, align 512
  %.partselect1 = trunc i64 %3 to i32
  %4 = call float @_llvm.fpga.unpack.bits.f32.i32(i32 %.partselect1)
  store float %4, float* %dst.0, align 512
  %dst.1 = getelementptr %struct.MyComplex, %struct.MyComplex* %dst, i64 0, i32 1
  %5 = lshr i64 %3, 32
  %.partselect = trunc i64 %5 to i32
  %6 = call float @_llvm.fpga.unpack.bits.f32.i32(i32 %.partselect)
  store float %6, float* %dst.1, align 4
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: alwaysinline nounwind readnone willreturn
define internal float @_llvm.fpga.unpack.bits.f32.i32(i32 %A) #4 {
  %A.cast = bitcast i32 %A to float
  ret float %A.cast
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0struct.MyComplex.333(i64* noalias align 512 %dst, %struct.MyComplex* noalias readonly %src) unnamed_addr #3 {
entry:
  %0 = icmp eq i64* %dst, null
  %1 = icmp eq %struct.MyComplex* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %src.0 = getelementptr %struct.MyComplex, %struct.MyComplex* %src, i64 0, i32 0
  %3 = load float, float* %src.0, align 4
  %4 = call i32 @_llvm.fpga.pack.bits.i32.f32(float %3)
  %5 = zext i32 %4 to i64
  %src.1 = getelementptr %struct.MyComplex, %struct.MyComplex* %src, i64 0, i32 1
  %6 = load float, float* %src.1, align 4
  %7 = call i32 @_llvm.fpga.pack.bits.i32.f32(float %6)
  %8 = zext i32 %7 to i64
  %9 = shl i64 %8, 32
  %.partset = or i64 %9, %5
  store i64 %.partset, i64* %dst, align 512
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: alwaysinline nounwind readnone willreturn
define internal i32 @_llvm.fpga.pack.bits.i32.f32(float %A) #4 {
  %A.cast = bitcast float %A to i32
  ret i32 %A.cast
}

declare i8* @malloc(i64)

declare void @free(i8*)

declare float @apatb_MHGD_detect_accel_hw(i64*, i32, i32, i32, i64*, i64*, float, i32, i32, i32, i64*)

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_back(%struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512, %struct.MyComplex* noalias, i64* noalias readonly align 512) unnamed_addr #2 {
entry:
  call fastcc void @onebyonecpy_hls.p0struct.MyComplex(%struct.MyComplex* %0, i64* align 512 %1)
  ret void
}

declare float @MHGD_detect_accel_hw_stub(%struct.MyComplex* noalias nonnull, i32, i32, i32, %struct.MyComplex* noalias nocapture nonnull readonly, %struct.MyComplex* noalias nocapture nonnull readonly, float, i32, i32, i32, %struct.MyComplex* noalias nocapture nonnull readonly)

define float @MHGD_detect_accel_hw_stub_wrapper(i64*, i32, i32, i32, i64*, i64*, float, i32, i32, i32, i64*) #5 {
entry:
  %11 = call i8* @malloc(i64 8)
  %12 = bitcast i8* %11 to %struct.MyComplex*
  %13 = call i8* @malloc(i64 8)
  %14 = bitcast i8* %13 to %struct.MyComplex*
  %15 = call i8* @malloc(i64 8)
  %16 = bitcast i8* %15 to %struct.MyComplex*
  %17 = call i8* @malloc(i64 8)
  %18 = bitcast i8* %17 to %struct.MyComplex*
  call void @copy_out(%struct.MyComplex* %12, i64* %0, %struct.MyComplex* %14, i64* %4, %struct.MyComplex* %16, i64* %5, %struct.MyComplex* %18, i64* %10)
  %19 = call float @MHGD_detect_accel_hw_stub(%struct.MyComplex* %12, i32 %1, i32 %2, i32 %3, %struct.MyComplex* %14, %struct.MyComplex* %16, float %6, i32 %7, i32 %8, i32 %9, %struct.MyComplex* %18)
  call void @copy_in(%struct.MyComplex* %12, i64* %0, %struct.MyComplex* %14, i64* %4, %struct.MyComplex* %16, i64* %5, %struct.MyComplex* %18, i64* %10)
  call void @free(i8* %11)
  call void @free(i8* %13)
  call void @free(i8* %15)
  call void @free(i8* %17)
  ret float %19
}

attributes #0 = { noinline willreturn "fpga.wrapper.func"="wrapper" }
attributes #1 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyin" }
attributes #2 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyout" }
attributes #3 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="onebyonecpy_hls" }
attributes #4 = { alwaysinline nounwind readnone willreturn }
attributes #5 = { "fpga.wrapper.func"="stub" }

!llvm.dbg.cu = !{}
!llvm.ident = !{!0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0}
!llvm.module.flags = !{!1, !2, !3}
!blackbox_cfg = !{!4}

!0 = !{!"clang version 7.0.0 "}
!1 = !{i32 2, !"Dwarf Version", i32 4}
!2 = !{i32 2, !"Debug Info Version", i32 3}
!3 = !{i32 1, !"wchar_size", i32 4}
!4 = !{}
