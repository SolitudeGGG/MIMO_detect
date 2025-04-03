// ==============================================================
// Vitis HLS - High-Level Synthesis from C, C++ and OpenCL v2024.2 (64-bit)
// Tool Version Limit: 2024.11
// Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
// Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
// 
// ==============================================================
`timescale 1ns/1ps
module MHGD_detect_accel_control_s_axi
#(parameter
    C_S_AXI_ADDR_WIDTH = 6,
    C_S_AXI_DATA_WIDTH = 32
)(
    input  wire                          ACLK,
    input  wire                          ARESET,
    input  wire                          ACLK_EN,
    input  wire [C_S_AXI_ADDR_WIDTH-1:0] AWADDR,
    input  wire                          AWVALID,
    output wire                          AWREADY,
    input  wire [C_S_AXI_DATA_WIDTH-1:0] WDATA,
    input  wire [C_S_AXI_DATA_WIDTH/8-1:0] WSTRB,
    input  wire                          WVALID,
    output wire                          WREADY,
    output wire [1:0]                    BRESP,
    output wire                          BVALID,
    input  wire                          BREADY,
    input  wire [C_S_AXI_ADDR_WIDTH-1:0] ARADDR,
    input  wire                          ARVALID,
    output wire                          ARREADY,
    output wire [C_S_AXI_DATA_WIDTH-1:0] RDATA,
    output wire [1:0]                    RRESP,
    output wire                          RVALID,
    input  wire                          RREADY,
    output wire [63:0]                   x_hat,
    output wire [63:0]                   H,
    output wire [63:0]                   y,
    output wire [63:0]                   v_tb
);
//------------------------Address Info-------------------
// Protocol Used: ap_ctrl_none
//
// 0x00 : reserved
// 0x04 : reserved
// 0x08 : reserved
// 0x0c : reserved
// 0x10 : Data signal of x_hat
//        bit 31~0 - x_hat[31:0] (Read/Write)
// 0x14 : Data signal of x_hat
//        bit 31~0 - x_hat[63:32] (Read/Write)
// 0x18 : reserved
// 0x1c : Data signal of H
//        bit 31~0 - H[31:0] (Read/Write)
// 0x20 : Data signal of H
//        bit 31~0 - H[63:32] (Read/Write)
// 0x24 : reserved
// 0x28 : Data signal of y
//        bit 31~0 - y[31:0] (Read/Write)
// 0x2c : Data signal of y
//        bit 31~0 - y[63:32] (Read/Write)
// 0x30 : reserved
// 0x34 : Data signal of v_tb
//        bit 31~0 - v_tb[31:0] (Read/Write)
// 0x38 : Data signal of v_tb
//        bit 31~0 - v_tb[63:32] (Read/Write)
// 0x3c : reserved
// (SC = Self Clear, COR = Clear on Read, TOW = Toggle on Write, COH = Clear on Handshake)

//------------------------Parameter----------------------
localparam
    ADDR_X_HAT_DATA_0 = 6'h10,
    ADDR_X_HAT_DATA_1 = 6'h14,
    ADDR_X_HAT_CTRL   = 6'h18,
    ADDR_H_DATA_0     = 6'h1c,
    ADDR_H_DATA_1     = 6'h20,
    ADDR_H_CTRL       = 6'h24,
    ADDR_Y_DATA_0     = 6'h28,
    ADDR_Y_DATA_1     = 6'h2c,
    ADDR_Y_CTRL       = 6'h30,
    ADDR_V_TB_DATA_0  = 6'h34,
    ADDR_V_TB_DATA_1  = 6'h38,
    ADDR_V_TB_CTRL    = 6'h3c,
    WRIDLE            = 2'd0,
    WRDATA            = 2'd1,
    WRRESP            = 2'd2,
    WRRESET           = 2'd3,
    RDIDLE            = 2'd0,
    RDDATA            = 2'd1,
    RDRESET           = 2'd2,
    ADDR_BITS                = 6;

//------------------------Local signal-------------------
    reg  [1:0]                    wstate = WRRESET;
    reg  [1:0]                    wnext;
    reg  [ADDR_BITS-1:0]          waddr;
    wire [C_S_AXI_DATA_WIDTH-1:0] wmask;
    wire                          aw_hs;
    wire                          w_hs;
    reg  [1:0]                    rstate = RDRESET;
    reg  [1:0]                    rnext;
    reg  [C_S_AXI_DATA_WIDTH-1:0] rdata;
    wire                          ar_hs;
    wire [ADDR_BITS-1:0]          raddr;
    // internal registers
    reg  [63:0]                   int_x_hat = 'b0;
    reg  [63:0]                   int_H = 'b0;
    reg  [63:0]                   int_y = 'b0;
    reg  [63:0]                   int_v_tb = 'b0;

//------------------------Instantiation------------------


//------------------------AXI write fsm------------------
assign AWREADY = (wstate == WRIDLE);
assign WREADY  = (wstate == WRDATA);
assign BRESP   = 2'b00;  // OKAY
assign BVALID  = (wstate == WRRESP);
assign wmask   = { {8{WSTRB[3]}}, {8{WSTRB[2]}}, {8{WSTRB[1]}}, {8{WSTRB[0]}} };
assign aw_hs   = AWVALID & AWREADY;
assign w_hs    = WVALID & WREADY;

// wstate
always @(posedge ACLK) begin
    if (ARESET)
        wstate <= WRRESET;
    else if (ACLK_EN)
        wstate <= wnext;
end

// wnext
always @(*) begin
    case (wstate)
        WRIDLE:
            if (AWVALID)
                wnext = WRDATA;
            else
                wnext = WRIDLE;
        WRDATA:
            if (WVALID)
                wnext = WRRESP;
            else
                wnext = WRDATA;
        WRRESP:
            if (BREADY)
                wnext = WRIDLE;
            else
                wnext = WRRESP;
        default:
            wnext = WRIDLE;
    endcase
end

// waddr
always @(posedge ACLK) begin
    if (ACLK_EN) begin
        if (aw_hs)
            waddr <= {AWADDR[ADDR_BITS-1:2], {2{1'b0}}};
    end
end

//------------------------AXI read fsm-------------------
assign ARREADY = (rstate == RDIDLE);
assign RDATA   = rdata;
assign RRESP   = 2'b00;  // OKAY
assign RVALID  = (rstate == RDDATA);
assign ar_hs   = ARVALID & ARREADY;
assign raddr   = ARADDR[ADDR_BITS-1:0];

// rstate
always @(posedge ACLK) begin
    if (ARESET)
        rstate <= RDRESET;
    else if (ACLK_EN)
        rstate <= rnext;
end

// rnext
always @(*) begin
    case (rstate)
        RDIDLE:
            if (ARVALID)
                rnext = RDDATA;
            else
                rnext = RDIDLE;
        RDDATA:
            if (RREADY & RVALID)
                rnext = RDIDLE;
            else
                rnext = RDDATA;
        default:
            rnext = RDIDLE;
    endcase
end

// rdata
always @(posedge ACLK) begin
    if (ACLK_EN) begin
        if (ar_hs) begin
            rdata <= 'b0;
            case (raddr)
                ADDR_X_HAT_DATA_0: begin
                    rdata <= int_x_hat[31:0];
                end
                ADDR_X_HAT_DATA_1: begin
                    rdata <= int_x_hat[63:32];
                end
                ADDR_H_DATA_0: begin
                    rdata <= int_H[31:0];
                end
                ADDR_H_DATA_1: begin
                    rdata <= int_H[63:32];
                end
                ADDR_Y_DATA_0: begin
                    rdata <= int_y[31:0];
                end
                ADDR_Y_DATA_1: begin
                    rdata <= int_y[63:32];
                end
                ADDR_V_TB_DATA_0: begin
                    rdata <= int_v_tb[31:0];
                end
                ADDR_V_TB_DATA_1: begin
                    rdata <= int_v_tb[63:32];
                end
            endcase
        end
    end
end


//------------------------Register logic-----------------
assign x_hat = int_x_hat;
assign H     = int_H;
assign y     = int_y;
assign v_tb  = int_v_tb;
// int_x_hat[31:0]
always @(posedge ACLK) begin
    if (ARESET)
        int_x_hat[31:0] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_X_HAT_DATA_0)
            int_x_hat[31:0] <= (WDATA[31:0] & wmask) | (int_x_hat[31:0] & ~wmask);
    end
end

// int_x_hat[63:32]
always @(posedge ACLK) begin
    if (ARESET)
        int_x_hat[63:32] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_X_HAT_DATA_1)
            int_x_hat[63:32] <= (WDATA[31:0] & wmask) | (int_x_hat[63:32] & ~wmask);
    end
end

// int_H[31:0]
always @(posedge ACLK) begin
    if (ARESET)
        int_H[31:0] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_H_DATA_0)
            int_H[31:0] <= (WDATA[31:0] & wmask) | (int_H[31:0] & ~wmask);
    end
end

// int_H[63:32]
always @(posedge ACLK) begin
    if (ARESET)
        int_H[63:32] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_H_DATA_1)
            int_H[63:32] <= (WDATA[31:0] & wmask) | (int_H[63:32] & ~wmask);
    end
end

// int_y[31:0]
always @(posedge ACLK) begin
    if (ARESET)
        int_y[31:0] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_Y_DATA_0)
            int_y[31:0] <= (WDATA[31:0] & wmask) | (int_y[31:0] & ~wmask);
    end
end

// int_y[63:32]
always @(posedge ACLK) begin
    if (ARESET)
        int_y[63:32] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_Y_DATA_1)
            int_y[63:32] <= (WDATA[31:0] & wmask) | (int_y[63:32] & ~wmask);
    end
end

// int_v_tb[31:0]
always @(posedge ACLK) begin
    if (ARESET)
        int_v_tb[31:0] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_V_TB_DATA_0)
            int_v_tb[31:0] <= (WDATA[31:0] & wmask) | (int_v_tb[31:0] & ~wmask);
    end
end

// int_v_tb[63:32]
always @(posedge ACLK) begin
    if (ARESET)
        int_v_tb[63:32] <= 0;
    else if (ACLK_EN) begin
        if (w_hs && waddr == ADDR_V_TB_DATA_1)
            int_v_tb[63:32] <= (WDATA[31:0] & wmask) | (int_v_tb[63:32] & ~wmask);
    end
end


//------------------------Memory logic-------------------

endmodule
