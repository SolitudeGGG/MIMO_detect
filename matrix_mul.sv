`timescale 1ns/1ps
`default_nettype none

// Parameterized fixed-point matrix multiplication with optional scalar mode
// The design reuses a small multiplier across the inner dimension to keep
// resource usage low on wide datapaths. When SERIAL_CHUNK is asserted the
// multiplier itself is decomposed into CHUNK_W pieces so it can map cleanly
// onto narrower DSP blocks.
module matrix_mul #(
    parameter integer ROWS        = 4,   // A rows
    parameter integer COLS        = 4,   // B columns / output columns
    parameter integer INNER       = 4,   // A cols / B rows
    parameter integer INT_W       = 8,
    parameter integer FRAC_W      = 8,
    parameter bit     SERIAL_CHUNK= 0,   // enable chunked multiply to save DSPs
    parameter integer CHUNK_W     = 16,  // DSP-friendly width (e.g., 16/17/18)
    // OUT_W keeps product growth plus log2(INNER+1) headroom for accumulation
    parameter integer OUT_W       = 2*(INT_W+FRAC_W)+$clog2(INNER+1)
)(
    input  wire clk,
    input  wire rst_n,
    input  wire start,
    input  wire scalar_mode,                 // when high: element-wise scalar scale
    input  wire signed [INT_W+FRAC_W-1:0] scalar,
    input  wire signed [ROWS*INNER*(INT_W+FRAC_W)-1:0] a_flat,
    input  wire signed [INNER*COLS*(INT_W+FRAC_W)-1:0] b_flat,
    output reg  signed [ROWS*COLS*OUT_W-1:0] c_flat,
    output reg  busy,
    output reg  done
);
    localparam integer ELEM_W = INT_W + FRAC_W;
    localparam integer PROD_W = 2 * ELEM_W;
    localparam integer ACC_W  = OUT_W;
    localparam integer ROW_W  = (ROWS  > 1) ? $clog2(ROWS)  : 1;
    localparam integer COL_W  = (COLS  > 1) ? $clog2(COLS)  : 1;
    localparam integer K_W    = (INNER > 1) ? $clog2(INNER) : 1;

    // index helpers ---------------------------------------------------------
    function automatic signed [ELEM_W-1:0] a_at(
        input integer r, input integer c
    );
        integer idx;
        begin
            idx = (r*INNER + c)*ELEM_W;
            a_at = a_flat[idx +: ELEM_W];
        end
    endfunction

    function automatic signed [OUT_W-1:0] scale_down(
        input signed [ACC_W-1:0] val
    );
        begin
            scale_down = val >>> FRAC_W;
        end
    endfunction

    function automatic signed [ACC_W-1:0] extend_prod(
        input signed [PROD_W-1:0] val
    );
        begin
            extend_prod = {{(ACC_W-PROD_W){val[PROD_W-1]}}, val};
        end
    endfunction

    function automatic signed [ELEM_W-1:0] b_at(
        input integer r, input integer c
    );
        integer idx;
        begin
            idx = (r*COLS + c)*ELEM_W;
            b_at = b_flat[idx +: ELEM_W];
        end
    endfunction

    task automatic set_c(
        input integer r, input integer c, input signed [OUT_W-1:0] val
    );
        integer idx;
        begin
            idx = (r*COLS + c)*OUT_W;
            c_flat[idx +: OUT_W] <= val;
        end
    endtask

    // datapath control ------------------------------------------------------
    reg [ROW_W-1:0] row_idx;
    reg [COL_W-1:0] col_idx;
    reg [K_W-1:0]   k_idx;
    reg signed [ACC_W-1:0] acc;

    reg signed [ELEM_W-1:0] mul_a;
    reg signed [ELEM_W-1:0] mul_b;
    reg mul_start;
    wire mul_busy;
    wire mul_valid;
    wire signed [PROD_W-1:0] mul_res;
    wire                     mul_ready = !mul_busy && !mul_valid;
    wire signed [ACC_W-1:0]  mac_sum = acc + extend_prod(mul_res);

    fixed_point_mult #(
        .WIDTH   (ELEM_W),
        .CHUNK_W (CHUNK_W),
        .USE_CHUNK(SERIAL_CHUNK)
    ) u_mul (
        .clk   (clk),
        .rst_n (rst_n),
        .start (mul_start),
        .a     (mul_a),
        .b     (mul_b),
        .busy  (mul_busy),
        .valid (mul_valid),
        .prod  (mul_res)
    );

    // FSM -------------------------------------------------------------------
    integer col_limit;
    always @(*) begin
        col_limit = scalar_mode ? INNER : COLS;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy     <= 1'b0;
            done     <= 1'b0;
            row_idx  <= '0;
            col_idx  <= '0;
            k_idx    <= '0;
            acc      <= '0;
            c_flat   <= '0;
            mul_start<= 1'b0;
            mul_a    <= '0;
            mul_b    <= '0;
        end else begin
            mul_start <= 1'b0;
            done      <= 1'b0;
            if (start && !busy) begin
                busy    <= 1'b1;
                row_idx <= '0;
                col_idx <= '0;
                k_idx   <= '0;
                acc     <= '0;
            end else if (busy) begin
                // ---------------- scalar mode (element-wise scale) ----------
                if (scalar_mode) begin
                    if (mul_ready && (col_idx < col_limit)) begin
                        mul_a     <= a_at(row_idx, col_idx);
                        mul_b     <= scalar;
                        mul_start <= 1'b1;
                    end
                    if (mul_valid) begin
                        set_c(row_idx, col_idx, scale_down(extend_prod(mul_res)));
                        if (col_idx == col_limit-1) begin
                            col_idx <= '0;
                            if (row_idx == ROWS-1) begin
                                busy <= 1'b0;
                                done <= 1'b1;
                            end else begin
                                row_idx <= row_idx + 1'b1;
                            end
                        end else begin
                            col_idx <= col_idx + 1'b1;
                        end
                    end
                end else begin
                // ---------------- standard matrix multiply ------------------
                    if (mul_ready) begin
                        mul_a     <= a_at(row_idx, k_idx);
                        mul_b     <= b_at(k_idx, col_idx);
                        mul_start <= 1'b1;
                    end
                    if (mul_valid) begin
                        if (k_idx == INNER-1) begin
                            set_c(row_idx, col_idx, scale_down(mac_sum));
                            acc     <= '0;
                            k_idx   <= '0;
                            if (col_idx == COLS-1) begin
                                col_idx <= '0;
                                if (row_idx == ROWS-1) begin
                                    busy <= 1'b0;
                                    done <= 1'b1;
                                end else begin
                                    row_idx <= row_idx + 1'b1;
                                end
                            end else begin
                                col_idx <= col_idx + 1'b1;
                            end
                        end else begin
                            acc   <= mac_sum;
                            k_idx <= k_idx + 1'b1;
                        end
                    end
                end
            end
        end
    end
endmodule

// Multiplier that optionally slices one operand into CHUNK_W pieces so each
// partial product fits typical DSP widths (e.g., 18x25) and reuses the same
// hardware over multiple cycles. When USE_CHUNK is 0 it behaves as a
// single-cycle multiplier to keep latency minimal.
module fixed_point_mult #(
    parameter integer WIDTH     = 16,
    parameter integer CHUNK_W   = 16,
    parameter bit     USE_CHUNK = 0
)(
    input  wire clk,
    input  wire rst_n,
    input  wire start,
    input  wire signed [WIDTH-1:0] a,
    input  wire signed [WIDTH-1:0] b,
    output reg  busy,
    output reg  valid,
    output reg  signed [2*WIDTH-1:0] prod
);
    localparam integer NUM_CHUNKS   = (WIDTH + CHUNK_W - 1)/CHUNK_W; // ceiling(WIDTH / CHUNK_W)
    localparam integer CHUNK_IDX_W  = (NUM_CHUNKS > 1) ? $clog2(NUM_CHUNKS) : 1;
    // ACC_W leaves CHUNK_W extra guard bits so chunk shifting cannot overflow
    localparam integer ACC_W        = (2*WIDTH) + CHUNK_W;
    // widen 1-bit chunk to carry explicit sign + magnitude; otherwise match CHUNK_W
    localparam integer BPIECE_W     = (CHUNK_W == 1) ? 2 : CHUNK_W;
    localparam integer PARTIAL_W    = WIDTH + BPIECE_W;

    reg signed [WIDTH-1:0] a_reg;
    reg signed [ACC_W-1:0] acc;
    reg [CHUNK_IDX_W-1:0] chunk_idx;

    reg  [WIDTH-1:0]           b_logic_work;
    reg                        b_sign;
    wire [CHUNK_W-1:0]         b_chunk_logic = b_logic_work[CHUNK_W-1:0];
    wire                       is_last_chunk = (chunk_idx == NUM_CHUNKS-1);
    wire signed [BPIECE_W-1:0] b_piece;
    generate
        if (CHUNK_W == 1) begin : g_chunk1
            assign b_piece = is_last_chunk ? {b_sign, b_chunk_logic[0]} : {1'b0, b_chunk_logic[0]};
        end else if (CHUNK_W > 1) begin : g_chunkN
            assign b_piece = is_last_chunk ? {b_sign, b_chunk_logic[BPIECE_W-2:0]} : {1'b0, b_chunk_logic[BPIECE_W-2:0]};
        end else begin : g_chunk_err
            initial $error("CHUNK_W must be >= 1");
        end
    endgenerate

    wire signed [PARTIAL_W-1:0] partial = a_reg * b_piece;
    wire signed [ACC_W-1:0] partial_ext = {{(ACC_W-PARTIAL_W){partial[PARTIAL_W-1]}}, partial};
    wire signed [ACC_W-1:0] acc_calc    = acc + (partial_ext <<< (chunk_idx*CHUNK_W));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy      <= 1'b0;
            valid     <= 1'b0;
            prod      <= '0;
            a_reg     <= '0;
            acc       <= '0;
            chunk_idx <= '0;
            b_logic_work <= '0;
            b_sign       <= 1'b0;
        end else begin
            valid <= 1'b0;
            if (start && !busy) begin
                if (USE_CHUNK) begin
                    busy      <= 1'b1;
                    a_reg     <= a;
                    acc       <= '0;
                    chunk_idx <= '0;
                    // capture raw bits for logical right stepping; sign is latched separately in b_sign
                    b_logic_work <= $unsigned(b);
                    b_sign       <= b[WIDTH-1];
                end else begin
                    // single-cycle path leverages synthesis DSP inference
                    prod  <= a * b;
                    valid <= 1'b1;
                end
            end else if (busy && USE_CHUNK) begin
                acc <= acc_calc;
                if (chunk_idx == NUM_CHUNKS-1) begin
                    prod      <= acc_calc[2*WIDTH-1:0];
                    valid     <= 1'b1;
                    busy      <= 1'b0;
                end else begin
                    chunk_idx <= chunk_idx + 1'b1;
                    b_logic_work <= b_logic_work >> CHUNK_W;
                end
            end
        end
    end
endmodule

`default_nettype wire
