// ==============================================================
// Generated by Vitis HLS v2024.2
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// ==============================================================

`timescale 1 ns / 1 ps 

(* CORE_GENERATION_INFO="MHGD_detect_accel_MHGD_detect_accel,hls_ip_2024_2,{HLS_INPUT_TYPE=c,HLS_INPUT_FLOAT=1,HLS_INPUT_FIXED=0,HLS_INPUT_PART=xcvu9p-flga2104-2-i,HLS_INPUT_CLOCK=10.000000,HLS_INPUT_ARCH=others,HLS_SYN_CLOCK=7.766000,HLS_SYN_LAT=-1,HLS_SYN_TPT=none,HLS_SYN_MEM=0,HLS_SYN_DSP=0,HLS_SYN_FF=956,HLS_SYN_LUT=2110,HLS_VERSION=2024_2}" *)

module MHGD_detect_accel (
        ap_clk,
        ap_rst,
        ap_start,
        ap_done,
        ap_idle,
        ap_ready,
        x_hat,
        Nt,
        Nr,
        mu,
        H,
        y,
        sigma2,
        mmse_init,
        lr_approx,
        iter,
        v_tb,
        ap_return
);

parameter    ap_ST_fsm_state1 = 12'd1;
parameter    ap_ST_fsm_state2 = 12'd2;
parameter    ap_ST_fsm_state3 = 12'd4;
parameter    ap_ST_fsm_state4 = 12'd8;
parameter    ap_ST_fsm_state5 = 12'd16;
parameter    ap_ST_fsm_state6 = 12'd32;
parameter    ap_ST_fsm_state7 = 12'd64;
parameter    ap_ST_fsm_state8 = 12'd128;
parameter    ap_ST_fsm_state9 = 12'd256;
parameter    ap_ST_fsm_state10 = 12'd512;
parameter    ap_ST_fsm_state11 = 12'd1024;
parameter    ap_ST_fsm_state12 = 12'd2048;

input   ap_clk;
input   ap_rst;
input   ap_start;
output   ap_done;
output   ap_idle;
output   ap_ready;
input  [63:0] x_hat;
input  [31:0] Nt;
input  [31:0] Nr;
input  [31:0] mu;
input  [63:0] H;
input  [63:0] y;
input  [31:0] sigma2;
input  [31:0] mmse_init;
input  [31:0] lr_approx;
input  [31:0] iter;
input  [63:0] v_tb;
output  [31:0] ap_return;

reg ap_done;
reg ap_idle;
reg ap_ready;

(* fsm_encoding = "none" *) reg   [11:0] ap_CS_fsm;
wire    ap_CS_fsm_state1;
reg   [31:0] iter_read_reg_123;
wire   [0:0] icmp_ln449_fu_92_p2;
reg   [0:0] icmp_ln449_reg_128;
wire   [0:0] and_ln413_fu_104_p2;
reg   [0:0] and_ln413_reg_133;
wire   [63:0] grp_fu_86_p2;
reg   [63:0] tmp_s_reg_137;
wire    ap_CS_fsm_state8;
wire   [31:0] grp_fu_83_p1;
reg   [31:0] p_acc_reg_142;
wire    ap_CS_fsm_state10;
wire   [30:0] trunc_ln519_fu_113_p1;
reg   [30:0] trunc_ln519_reg_148;
wire    ap_CS_fsm_state11;
wire   [22:0] trunc_ln519_1_fu_118_p1;
reg   [22:0] trunc_ln519_1_reg_153;
wire    grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start;
wire    grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_done;
wire    grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_idle;
wire    grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_ready;
wire   [0:0] grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_return;
reg    grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg;
wire    ap_CS_fsm_state12;
wire    ap_CS_fsm_state9;
wire   [63:0] grp_fu_86_p1;
wire   [0:0] icmp_ln434_fu_98_p2;
wire   [31:0] bitcast_ln519_fu_110_p1;
reg    ap_block_state12_on_subcall_done;
reg   [11:0] ap_NS_fsm;
reg    ap_ST_fsm_state1_blk;
wire    ap_ST_fsm_state2_blk;
wire    ap_ST_fsm_state3_blk;
wire    ap_ST_fsm_state4_blk;
wire    ap_ST_fsm_state5_blk;
wire    ap_ST_fsm_state6_blk;
wire    ap_ST_fsm_state7_blk;
wire    ap_ST_fsm_state8_blk;
wire    ap_ST_fsm_state9_blk;
wire    ap_ST_fsm_state10_blk;
wire    ap_ST_fsm_state11_blk;
reg    ap_ST_fsm_state12_blk;
wire    ap_ce_reg;

// power-on initialization
initial begin
#0 ap_CS_fsm = 12'd1;
#0 grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg = 1'b0;
end

MHGD_detect_accel_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3 grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74(
    .ap_clk(ap_clk),
    .ap_rst(ap_rst),
    .ap_start(grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start),
    .ap_done(grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_done),
    .ap_idle(grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_idle),
    .ap_ready(grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_ready),
    .iter(iter_read_reg_123),
    .bitcast_ln519(trunc_ln519_reg_148),
    .empty(trunc_ln519_1_reg_153),
    .p_acc(p_acc_reg_142),
    .icmp_ln449(icmp_ln449_reg_128),
    .ap_return(grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_return)
);

MHGD_detect_accel_fptrunc_64ns_32_2_no_dsp_1 #(
    .ID( 1 ),
    .NUM_STAGE( 2 ),
    .din0_WIDTH( 64 ),
    .dout_WIDTH( 32 ))
fptrunc_64ns_32_2_no_dsp_1_U8(
    .clk(ap_clk),
    .reset(ap_rst),
    .din0(tmp_s_reg_137),
    .ce(1'b1),
    .dout(grp_fu_83_p1)
);

MHGD_detect_accel_dexp_64ns_64ns_64_8_full_dsp_1 #(
    .ID( 1 ),
    .NUM_STAGE( 8 ),
    .din0_WIDTH( 64 ),
    .din1_WIDTH( 64 ),
    .dout_WIDTH( 64 ))
dexp_64ns_64ns_64_8_full_dsp_1_U9(
    .clk(ap_clk),
    .reset(ap_rst),
    .din0(64'd0),
    .din1(grp_fu_86_p1),
    .ce(1'b1),
    .dout(grp_fu_86_p2)
);

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_CS_fsm <= ap_ST_fsm_state1;
    end else begin
        ap_CS_fsm <= ap_NS_fsm;
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg <= 1'b0;
    end else begin
        if ((1'b1 == ap_CS_fsm_state11)) begin
            grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg <= 1'b1;
        end else if ((grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_ready == 1'b1)) begin
            grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg <= 1'b0;
        end
    end
end

always @ (posedge ap_clk) begin
    if ((1'b1 == ap_CS_fsm_state1)) begin
        and_ln413_reg_133 <= and_ln413_fu_104_p2;
        icmp_ln449_reg_128 <= icmp_ln449_fu_92_p2;
        iter_read_reg_123 <= iter;
    end
end

always @ (posedge ap_clk) begin
    if ((1'b1 == ap_CS_fsm_state10)) begin
        p_acc_reg_142 <= grp_fu_83_p1;
    end
end

always @ (posedge ap_clk) begin
    if ((1'b1 == ap_CS_fsm_state8)) begin
        tmp_s_reg_137 <= grp_fu_86_p2;
    end
end

always @ (posedge ap_clk) begin
    if ((1'b1 == ap_CS_fsm_state11)) begin
        trunc_ln519_1_reg_153 <= trunc_ln519_1_fu_118_p1;
        trunc_ln519_reg_148 <= trunc_ln519_fu_113_p1;
    end
end

assign ap_ST_fsm_state10_blk = 1'b0;

assign ap_ST_fsm_state11_blk = 1'b0;

always @ (*) begin
    if ((1'b1 == ap_block_state12_on_subcall_done)) begin
        ap_ST_fsm_state12_blk = 1'b1;
    end else begin
        ap_ST_fsm_state12_blk = 1'b0;
    end
end

always @ (*) begin
    if ((ap_start == 1'b0)) begin
        ap_ST_fsm_state1_blk = 1'b1;
    end else begin
        ap_ST_fsm_state1_blk = 1'b0;
    end
end

assign ap_ST_fsm_state2_blk = 1'b0;

assign ap_ST_fsm_state3_blk = 1'b0;

assign ap_ST_fsm_state4_blk = 1'b0;

assign ap_ST_fsm_state5_blk = 1'b0;

assign ap_ST_fsm_state6_blk = 1'b0;

assign ap_ST_fsm_state7_blk = 1'b0;

assign ap_ST_fsm_state8_blk = 1'b0;

assign ap_ST_fsm_state9_blk = 1'b0;

always @ (*) begin
    if (((1'd0 == and_ln413_reg_133) & (grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_return == 1'd1) & (icmp_ln449_reg_128 == 1'd0) & (1'b1 == ap_CS_fsm_state12) & (1'b0 == ap_block_state12_on_subcall_done))) begin
        ap_done = 1'b1;
    end else begin
        ap_done = 1'b0;
    end
end

always @ (*) begin
    if (((1'b1 == ap_CS_fsm_state1) & (ap_start == 1'b0))) begin
        ap_idle = 1'b1;
    end else begin
        ap_idle = 1'b0;
    end
end

always @ (*) begin
    if (((1'd0 == and_ln413_reg_133) & (grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_return == 1'd1) & (icmp_ln449_reg_128 == 1'd0) & (1'b1 == ap_CS_fsm_state12) & (1'b0 == ap_block_state12_on_subcall_done))) begin
        ap_ready = 1'b1;
    end else begin
        ap_ready = 1'b0;
    end
end

always @ (*) begin
    case (ap_CS_fsm)
        ap_ST_fsm_state1 : begin
            if (((1'b1 == ap_CS_fsm_state1) & (ap_start == 1'b1) & (1'd1 == and_ln413_fu_104_p2))) begin
                ap_NS_fsm = ap_ST_fsm_state12;
            end else if (((1'd0 == and_ln413_fu_104_p2) & (1'b1 == ap_CS_fsm_state1) & (ap_start == 1'b1))) begin
                ap_NS_fsm = ap_ST_fsm_state2;
            end else begin
                ap_NS_fsm = ap_ST_fsm_state1;
            end
        end
        ap_ST_fsm_state2 : begin
            ap_NS_fsm = ap_ST_fsm_state3;
        end
        ap_ST_fsm_state3 : begin
            ap_NS_fsm = ap_ST_fsm_state4;
        end
        ap_ST_fsm_state4 : begin
            ap_NS_fsm = ap_ST_fsm_state5;
        end
        ap_ST_fsm_state5 : begin
            ap_NS_fsm = ap_ST_fsm_state6;
        end
        ap_ST_fsm_state6 : begin
            ap_NS_fsm = ap_ST_fsm_state7;
        end
        ap_ST_fsm_state7 : begin
            ap_NS_fsm = ap_ST_fsm_state8;
        end
        ap_ST_fsm_state8 : begin
            ap_NS_fsm = ap_ST_fsm_state9;
        end
        ap_ST_fsm_state9 : begin
            ap_NS_fsm = ap_ST_fsm_state10;
        end
        ap_ST_fsm_state10 : begin
            ap_NS_fsm = ap_ST_fsm_state11;
        end
        ap_ST_fsm_state11 : begin
            ap_NS_fsm = ap_ST_fsm_state12;
        end
        ap_ST_fsm_state12 : begin
            if (((1'd0 == and_ln413_reg_133) & (grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_return == 1'd1) & (icmp_ln449_reg_128 == 1'd0) & (1'b1 == ap_CS_fsm_state12) & (1'b0 == ap_block_state12_on_subcall_done))) begin
                ap_NS_fsm = ap_ST_fsm_state1;
            end else begin
                ap_NS_fsm = ap_ST_fsm_state12;
            end
        end
        default : begin
            ap_NS_fsm = 'bx;
        end
    endcase
end

assign and_ln413_fu_104_p2 = (icmp_ln449_fu_92_p2 & icmp_ln434_fu_98_p2);

assign ap_CS_fsm_state1 = ap_CS_fsm[32'd0];

assign ap_CS_fsm_state10 = ap_CS_fsm[32'd9];

assign ap_CS_fsm_state11 = ap_CS_fsm[32'd10];

assign ap_CS_fsm_state12 = ap_CS_fsm[32'd11];

assign ap_CS_fsm_state8 = ap_CS_fsm[32'd7];

assign ap_CS_fsm_state9 = ap_CS_fsm[32'd8];

always @ (*) begin
    ap_block_state12_on_subcall_done = ((1'd0 == and_ln413_reg_133) & (grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_done == 1'b0));
end

assign ap_return = 32'd0;

assign bitcast_ln519_fu_110_p1 = p_acc_reg_142;

assign grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start = grp_MHGD_detect_accel_Pipeline_VITIS_LOOP_480_3_fu_74_ap_start_reg;

assign grp_fu_86_p1 = 'bx;

assign icmp_ln434_fu_98_p2 = ((mmse_init == 32'd0) ? 1'b1 : 1'b0);

assign icmp_ln449_fu_92_p2 = (($signed(Nt) > $signed(32'd0)) ? 1'b1 : 1'b0);

assign trunc_ln519_1_fu_118_p1 = bitcast_ln519_fu_110_p1[22:0];

assign trunc_ln519_fu_113_p1 = bitcast_ln519_fu_110_p1[30:0];

endmodule //MHGD_detect_accel
