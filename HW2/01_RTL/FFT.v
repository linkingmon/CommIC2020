
// multi path delay commutator FFT (top module)
module FFT(
    clk,
    rst_n,
    in_valid,
    din_r,
    din_i,
    out_valid,
    dout_r,
    dout_i
);

parameter n_points = 32;
parameter n_levels = $clog2(n_points);

parameter S_IDLE = 0;
parameter S_START = 1;
parameter S_BUFFER = 2;
parameter S_DONE = 3;

input clk;
input rst_n;
input in_valid;
input signed [11:0] din_r;
input signed [11:0] din_i;
output out_valid;
output signed [15:0] dout_r;
output signed [15:0] dout_i;

// registers
reg [n_levels:0] counter_r, counter_w;
reg signed [11:0] din_r_r;
reg signed [11:0] din_i_r;
reg [1:0] state_r, state_w;

// commutator wires
wire signed [15:0] commut_in_top_r [0:n_levels-1];
wire signed [15:0] commut_in_top_i [0:n_levels-1];
wire signed [15:0] commut_out_top_r [0:n_levels-1];
wire signed [15:0] commut_out_top_i [0:n_levels-1];

wire signed [15:0] commut_in_bottom_r [0:n_levels-1];
wire signed [15:0] commut_in_bottom_i [0:n_levels-1];
wire signed [15:0] commut_out_bottom_r [0:n_levels-1];
wire signed [15:0] commut_out_bottom_i [0:n_levels-1];

// butterfly wires
wire signed [15:0] butterfly_in_top_r [0:n_levels-1];
wire signed [15:0] butterfly_in_top_i [0:n_levels-1];

wire signed [15:0] butterfly_out_bottom_r [0:n_levels-1];
wire signed [15:0] butterfly_out_bottom_i [0:n_levels-1];

// mult wires
wire signed [15:0] twiddle_mult_out_r [0:n_levels-2];
wire signed [15:0] twiddle_mult_out_i [0:n_levels-2];

// output buffer wires
wire signed [15:0] buff_out1_r, buff_out1_i, buff_out2_r, buff_out2_i;

// control wires
wire swap [0:n_levels-1];
wire buffer_valid;

assign out_valid = (state_r == S_DONE) ? 1 : 0;
assign buffer_valid = (state_r == S_BUFFER) ? 1 : 0;

// commutator(in1_r, in1_i, in2_r, in2_i, out1_r, out1_i, out2_r, out2_i, swap)
commutator commutator_level4({{4{din_r_r[11]}},din_r_r}, {{4{din_i_r[11]}},din_i_r}, 16'b0, 16'b0, 
    commut_out_top_r[4], commut_out_top_i[4], commut_out_bottom_r[4], commut_out_bottom_i[4], swap[4]);
commutator commutator_level3(commut_in_top_r[3], commut_in_top_i[3], commut_in_bottom_r[3], commut_in_bottom_i[3], 
    commut_out_top_r[3], commut_out_top_i[3], commut_out_bottom_r[3], commut_out_bottom_i[3], swap[3]);
commutator commutator_level2(commut_in_top_r[2], commut_in_top_i[2], commut_in_bottom_r[2], commut_in_bottom_i[2], 
    commut_out_top_r[2], commut_out_top_i[2], commut_out_bottom_r[2], commut_out_bottom_i[2], swap[2]);
commutator commutator_level1(commut_in_top_r[1], commut_in_top_i[1], commut_in_bottom_r[1], commut_in_bottom_i[1], 
    commut_out_top_r[1], commut_out_top_i[1], commut_out_bottom_r[1], commut_out_bottom_i[1], swap[1]);
commutator commutator_level0(commut_in_top_r[0], commut_in_top_i[0], commut_in_bottom_r[0], commut_in_bottom_i[0], 
    commut_out_top_r[0], commut_out_top_i[0], commut_out_bottom_r[0], commut_out_bottom_i[0], swap[0]);

// delay_line(clk, rst_n, in_r, in_i, out_r, out_i);
delay_line #(.delay_nums(16)) delay_top_level4(clk, rst_n, commut_out_top_r[4], commut_out_top_i[4], butterfly_in_top_r[4], butterfly_in_top_i[4]);
delay_line #(.delay_nums(8)) delay_top_level3(clk, rst_n, commut_out_top_r[3], commut_out_top_i[3], butterfly_in_top_r[3], butterfly_in_top_i[3]);
delay_line #(.delay_nums(4)) delay_top_level2(clk, rst_n, commut_out_top_r[2], commut_out_top_i[2], butterfly_in_top_r[2], butterfly_in_top_i[2]);
delay_line #(.delay_nums(2)) delay_top_level1(clk, rst_n, commut_out_top_r[1], commut_out_top_i[1], butterfly_in_top_r[1], butterfly_in_top_i[1]);
delay_line #(.delay_nums(1)) delay_top_level0(clk, rst_n, commut_out_top_r[0], commut_out_top_i[0], butterfly_in_top_r[0], butterfly_in_top_i[0]);

delay_line #(.delay_nums(8)) delay_bottom_level3(clk, rst_n, twiddle_mult_out_r[3], twiddle_mult_out_i[3], commut_in_bottom_r[3], commut_in_bottom_i[3]);
delay_line #(.delay_nums(4)) delay_bottom_level2(clk, rst_n, twiddle_mult_out_r[2], twiddle_mult_out_i[2], commut_in_bottom_r[2], commut_in_bottom_i[2]);
delay_line #(.delay_nums(2)) delay_bottom_level1(clk, rst_n, twiddle_mult_out_r[1], twiddle_mult_out_i[1], commut_in_bottom_r[1], commut_in_bottom_i[1]);
delay_line #(.delay_nums(1)) delay_bottom_level0(clk, rst_n, twiddle_mult_out_r[0], twiddle_mult_out_i[0], commut_in_bottom_r[0], commut_in_bottom_i[0]);

// module butterfly(in1_r, in1_i, in2_r, in2_i, out1_r, out1_i, out2_r, out2_i)
butterfly butterfly_level4(butterfly_in_top_r[4], butterfly_in_top_i[4], commut_out_bottom_r[4], commut_out_bottom_i[4],
    commut_in_top_r[3], commut_in_top_i[3], butterfly_out_bottom_r[4], butterfly_out_bottom_i[4]);
butterfly butterfly_level3(butterfly_in_top_r[3], butterfly_in_top_i[3], commut_out_bottom_r[3], commut_out_bottom_i[3],
    commut_in_top_r[2], commut_in_top_i[2], butterfly_out_bottom_r[3], butterfly_out_bottom_i[3]);
butterfly butterfly_level2(butterfly_in_top_r[2], butterfly_in_top_i[2], commut_out_bottom_r[2], commut_out_bottom_i[2],
    commut_in_top_r[1], commut_in_top_i[1], butterfly_out_bottom_r[2], butterfly_out_bottom_i[2]);
butterfly butterfly_level1(butterfly_in_top_r[1], butterfly_in_top_i[1], commut_out_bottom_r[1], commut_out_bottom_i[1],
    commut_in_top_r[0], commut_in_top_i[0], butterfly_out_bottom_r[1], butterfly_out_bottom_i[1]);
butterfly butterfly_level0(butterfly_in_top_r[0], butterfly_in_top_i[0], commut_out_bottom_r[0], commut_out_bottom_i[0],
    buff_out1_r, buff_out1_i, buff_out2_r,buff_out2_i);

// module output_buffer(clk, rst_n in1_r, in1_i, in2_r, in2_i, i_valid, out_r, out_i)
output_buffer output_buffer_0(clk, rst_n, buff_out1_r, buff_out1_i, buff_out2_r, buff_out2_i, buffer_valid, dout_r, dout_i);

// module twidle_mult(in_point, in_r,　in_i,　out_r,　out_i);
twidle_mult3 twidle_mult3_0(counter_r[3:0], butterfly_out_bottom_r[4], butterfly_out_bottom_i[4], twiddle_mult_out_r[3], twiddle_mult_out_i[3]);
twidle_mult2 twidle_mult2_0(counter_r[2:0], butterfly_out_bottom_r[3], butterfly_out_bottom_i[3], twiddle_mult_out_r[2], twiddle_mult_out_i[2]);
twidle_mult1 twidle_mult1_0(counter_r[1:0], butterfly_out_bottom_r[2], butterfly_out_bottom_i[2], twiddle_mult_out_r[1], twiddle_mult_out_i[1]);
twidle_mult0 twidle_mult0_0(counter_r[0], butterfly_out_bottom_r[1], butterfly_out_bottom_i[1], twiddle_mult_out_r[0], twiddle_mult_out_i[0]);

assign swap[4] = counter_r[4]; 
assign swap[3] = counter_r[3]; 
assign swap[2] = counter_r[2]; 
assign swap[1] = counter_r[1]; 
assign swap[0] = counter_r[0]; 

always@(*) begin
    state_w = state_r;
    counter_w = counter_r;
    if(state_r != S_IDLE)
        counter_w = counter_r + 1;
    if(in_valid == 1)
        state_w = S_START;
    if(counter_r == 29)
        state_w = S_BUFFER;
    if(counter_r == 40)
        state_w = S_DONE;
end

always@(posedge clk or negedge rst_n) begin
    if(!rst_n) begin
        counter_r <= 0;
        state_r <= S_IDLE;
    end
    else begin
        counter_r <= counter_w;
        din_r_r <= din_r;
        din_i_r <= din_i;
        state_r <= state_w;
    end
end

endmodule

// delay line module
module delay_line(
    clk,
    rst_n,
    in_r,
    in_i,
    out_r,
    out_i
);

parameter data_bits = 16;
parameter delay_nums = 4;

input clk, rst_n;
input signed [data_bits-1:0] in_r, in_i;
output signed [data_bits-1:0] out_r, out_i;

reg signed [data_bits-1:0] buffers_r [0:delay_nums-1], buffers_i [0:delay_nums-1];
integer i;

assign out_r = buffers_r[delay_nums-1];
assign out_i = buffers_i[delay_nums-1];

always@(posedge clk or negedge rst_n) begin
    if(!rst_n) begin
        for(i=0;i<delay_nums;i=i+1)
            buffers_r[i] <= 0;
        for(i=0;i<delay_nums;i=i+1)
            buffers_i[i] <= 0;
    end
    else begin
        buffers_r[0] <= in_r;
        buffers_i[0] <= in_i;
        for(i=1;i<delay_nums;i=i+1)
            buffers_r[i] <= buffers_r[i-1];
        for(i=1;i<delay_nums;i=i+1)
            buffers_i[i] <= buffers_i[i-1];
    end
end

endmodule

// commutator module
module commutator(
    in1_r,
    in1_i,
    in2_r,
    in2_i,
    out1_r,
    out1_i,
    out2_r,
    out2_i,
    swap
);

parameter data_bits = 16;

input swap;
input signed [data_bits-1:0] in1_r, in1_i, in2_r, in2_i;
output signed [data_bits-1:0] out1_r, out1_i, out2_r, out2_i;

assign out1_r = (swap == 1) ? in2_r : in1_r;
assign out1_i = (swap == 1) ? in2_i : in1_i;
assign out2_r = (swap == 1) ? in1_r : in2_r;
assign out2_i = (swap == 1) ? in1_i : in2_i;

endmodule

// radix 2 butterfly
module butterfly(
    in1_r,
    in1_i,
    in2_r,
    in2_i,
    out1_r,
    out1_i,
    out2_r,
    out2_i
);

parameter data_bits = 16;

input signed [data_bits-1:0] in1_r, in1_i, in2_r, in2_i;
output signed [data_bits-1:0] out1_r, out1_i, out2_r, out2_i;

assign out1_r = in1_r + in2_r;
assign out1_i = in1_i + in2_i;
assign out2_r = in1_r - in2_r;
assign out2_i = in1_i - in2_i;

endmodule

// twidle multiplier level 1 (degree: 0, 90)
module twidle_mult0(
    in_point,
    in_r,
    in_i,
    out_r,
    out_i
);

parameter data_bits = 16;


input in_point;
input signed [data_bits-1:0] in_r, in_i;
output reg signed [data_bits-1:0] out_r, out_i;

always@(*) begin
    case(in_point)
        0: // 0
        begin
            out_r = in_r;
            out_i = in_i;
        end
        1: // 90
        begin
            out_r = in_i;
            out_i = -in_r;
        end
    endcase
end
endmodule

// twidle multiplier level 2 (degree: 0, 45, 90, 135)
// 1/2 + 1/8 + 1/16 = 0.6875
// 1/2 + 1/8 + 1/16 + 1/128 = 0.6953
module twidle_mult1(
    in_point,
    in_r,
    in_i,
    out_r,
    out_i
);

parameter data_bits = 16;

input [1:0] in_point;
input signed [data_bits-1:0] in_r, in_i;
output reg signed [data_bits-1:0] out_r, out_i;


always@(*) begin
    case(in_point)
        0: // 0
        begin   
            out_r = in_r;
            out_i = in_i;
        end
        1: // 45
        begin
            out_r =  in_r - (in_r >>> 2) - (in_r >>> 4) + (in_r >>> 6) + (in_r >>> 8)
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) + (in_i >>> 8) ;
            out_i = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) - (in_r >>> 8) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) + (in_i >>> 8);
        end
        2: // 90
        begin
            out_r = in_i;
            out_i = -in_r;
        end
        3: // 135
        begin
            out_r = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) - (in_r >>> 8) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) + (in_i >>> 8);
            out_i = - (in_r >>> 1) - (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) - (in_r >>> 8)
                - in_i + (in_i >>> 2) + (in_i >>> 4) - (in_i >>> 6) - (in_i >>> 8);
        end
    endcase
end

endmodule

// twidle multiplier level 3 (degree: 0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5)
// rotate 21.8014 degree, 
// [ 5 -2 ] [a] = [a']
// [ 2  5 ] [b] = [b']
// 1/8 + 1/32 = 0.1563
// 1/8 + 1/32 + 1/64 = 0.1719
// a' = (4 + 1) * a * (1/8 + 1/32) + (-2) * b * (1/8 + 1/32)
module twidle_mult2(
    in_point,
    in_r,
    in_i,
    out_r,
    out_i
);

parameter data_bits = 16;

input [2:0] in_point;
input signed [data_bits-1:0] in_r, in_i;
output reg signed [data_bits-1:0] out_r, out_i;


always@(*) begin
    case(in_point)
        0: // 0
        begin
            out_r = in_r;
            out_i = in_i;
        end
        1: // 22.5
        begin
            out_r = in_r - (in_r >>> 4) - (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
            out_i = - (in_r >>> 1) + (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
        end
        2: // 45
        begin
            out_r =  in_r - (in_r >>> 2) - (in_r >>> 4) + (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
            out_i = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
        end
        3: // 67.5
        begin
            out_r = (in_r >>> 1) - (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
            out_i = -in_r + (in_r >>> 4) + (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
        end
        4: // 90
        begin
            out_r = in_i;
            out_i = -in_r;
        end
        5: // 112.5
        begin
            out_r = -(in_r >>> 1) + (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
            out_i = -in_r + (in_r >>> 4) + (in_r >>> 6) - (in_i >>> 1) + (in_i >>> 3);
        end
        6: // 135
        begin
            out_r = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
            out_i = - (in_r >>> 1) - (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                - in_i + (in_i >>> 2) + (in_i >>> 4) - (in_i >>> 6) ;
        end
        7: // 157.5
        begin
            out_r = -in_r + (in_r >>> 4) + (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
            out_i = - (in_r >>> 1) + (in_r >>> 3) - in_i +  (in_i >>> 4) + (in_i >>> 6);
        end
    endcase
end

endmodule

module twidle_mult3(
    in_point,
    in_r,
    in_i,
    out_r,
    out_i
);

parameter data_bits = 16;

input [3:0] in_point;
input signed [data_bits-1:0] in_r, in_i;
output reg signed [data_bits-1:0] out_r, out_i;


always@(*) begin
    case(in_point)
        0: // 0
        begin
            out_r = in_r;
            out_i = in_i;
        end
        1: // 11.25
        begin
            out_r = in_r - (in_r >>> 6) + (in_i >>> 2) - (in_i >>> 4);
            out_i = - (in_r >>> 2) + (in_r >>> 4) + in_i - (in_i >>> 6);
        end
        2: // 22.5
        begin
            out_r = in_r - (in_r >>> 4) - (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
            out_i = - (in_r >>> 1) + (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
        end
        3: // 33.75
        begin
            out_r = in_r - (in_r >>> 2) + (in_r >>> 4) + (in_r >>> 6) + (in_i >>> 1) + (in_i >>> 4) - (in_i >>> 6);
            out_i = - (in_r >>> 1) - (in_r >>> 4) + (in_r >>> 6) + in_i - (in_i >>> 2) + (in_i >>> 4) + (in_i >>> 6);
        end
        4: // 45
        begin
            out_r =  in_r - (in_r >>> 2) - (in_r >>> 4) + (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
            out_i = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
        end
        5: // 56.25
        begin
            out_r = (in_r >>> 1) + (in_r >>> 4) - (in_r >>> 6) + in_i - (in_i >>> 2) + (in_i >>> 4) ;
            out_i = - in_r + (in_r >>> 2) - (in_r >>> 4) - (in_r >>> 6) + (in_i >>> 1) + (in_i >>> 4) ;
        end
        6: // 67.5
        begin
            out_r = (in_r >>> 1) - (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
            out_i = -in_r + (in_r >>> 4) + (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
        end
        7: // 78.25
        begin
            out_r =   (in_r >>> 2) - (in_r >>> 4) + in_i - (in_i >>> 5);
            out_i = - in_r + (in_r >>> 5) + (in_i >>> 2) - (in_i >>> 4);
        end
        8: // 90
        begin
            out_r = in_i;
            out_i = -in_r;
        end
        9: // 101.25
        begin
            out_r = - (in_r >>> 2) + (in_r >>> 4) + in_i - (in_i >>> 5);
            out_i = - in_r + (in_r >>> 5) - (in_i >>> 2) + (in_i >>> 4);
        end
        10: // 112.5
        begin
            out_r = -(in_r >>> 1) + (in_r >>> 3) + in_i - (in_i >>> 4) - (in_i >>> 6);
            out_i = -in_r + (in_r >>> 4) + (in_r >>> 6) - (in_i >>> 1) + (in_i >>> 3);
        end
        11: // 123.75
        begin
            out_r = - (in_r >>> 1) - (in_r >>> 4) + (in_r >>> 6) + in_i - (in_i >>> 2) + (in_i >>> 4) + (in_i >>> 6);
            out_i = - in_r + (in_r >>> 2) - (in_r >>> 4) - (in_r >>> 6) - (in_i >>> 1) - (in_i >>> 4) + (in_i >>> 6);
        end
        12: // 135
        begin
            out_r = - in_r + (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                + in_i - (in_i >>> 2) - (in_i >>> 4) + (in_i >>> 6) ;
            out_i = - (in_r >>> 1) - (in_r >>> 2) + (in_r >>> 4) - (in_r >>> 6) 
                - in_i + (in_i >>> 2) + (in_i >>> 4) - (in_i >>> 6) ;
        end
        13: // 146.25
        begin
            out_r = - in_r + (in_r >>> 2) - (in_r >>> 4) - (in_r >>> 6) + (in_i >>> 1) + (in_i >>> 4) - (in_i >>> 6);
            out_i = - (in_r >>> 1) - (in_r >>> 4) + (in_r >>> 6) - in_i + (in_i >>> 2) - (in_i >>> 4) - (in_i >>> 6);
        end
        14: // 157.5
        begin
            out_r = -in_r + (in_r >>> 4) + (in_r >>> 6) + (in_i >>> 1) - (in_i >>> 3);
            out_i = - (in_r >>> 1) + (in_r >>> 3) - in_i +  (in_i >>> 4) + (in_i >>> 6);
        end
        15: // 168.75
        begin
            out_r = - in_r + (in_r >>> 5) + (in_i >>> 2) - (in_i >>> 4);
            out_i = - (in_r >>> 2) + (in_r >>> 4) - in_i + (in_i >>> 5);
        end
    endcase
end

endmodule

// output buffering
module output_buffer(
    clk,
    rst_n,
    in1_r,
    in1_i, 
    in2_r, 
    in2_i, 
    i_valid, 
    out_r, 
    out_i
);

parameter data_bits = 16;
parameter n_points = 32;

input clk, rst_n, i_valid;
input signed [data_bits-1:0] in1_r, in1_i, in2_r, in2_i;
output signed [data_bits-1:0] out_r, out_i;

parameter S_IDLE = 0;
parameter S_START = 1;
parameter S_DONE = 2;

integer i;

reg [data_bits-1:0] buffer_r_r [0:n_points-1], buffer_r_w [0:n_points-1];
reg [data_bits-1:0] buffer_i_r [0:n_points-1], buffer_i_w [0:n_points-1];

reg [$clog2(n_points)-1:0] counter_r, counter_w;
reg [$clog2(n_points)-1:0] counter2_r, counter2_w;
reg [$clog2(n_points)-1:0] counter3_r, counter3_w;
reg [$clog2(n_points)-1:0] reverse_w, reverse2_w;
reg [1:0] state_r, state_w;

assign out_r = buffer_r_r[counter3_r];
assign out_i = buffer_i_r[counter3_r];

always@(*) begin
    state_w = state_r;
    counter_w = counter_r;
    counter2_w = counter2_r;
    counter3_w = counter3_r;
    reverse_w = {counter_r[0],counter_r[1],counter_r[2],counter_r[3],counter_r[4]};
    reverse2_w = {counter2_r[0],counter2_r[1],counter2_r[2],counter2_r[3],counter2_r[4]};
    for(i=0;i<n_points;i=i+1) begin
        buffer_r_w[i] = buffer_r_r[i];
        buffer_i_w[i] = buffer_i_r[i];
    end
    if(i_valid == 1)
        state_w = S_START;
    if(state_r == S_START) begin
        counter_w = counter_r + 2;
        counter2_w = counter2_r + 2;
    end
    if(counter2_r == 31)
        state_w = S_DONE;
    if(state_r ==  S_START) begin
        buffer_r_w[reverse_w] = in1_r;
        buffer_i_w[reverse_w] = in1_i;
        buffer_r_w[reverse2_w] = in2_r;
        buffer_i_w[reverse2_w] = in2_i;
    end
    if(state_r != S_IDLE)
        counter3_w = counter3_r + 1;

end

always@(posedge clk or negedge rst_n) begin
    if(!rst_n) begin
        state_r <= S_IDLE;
        counter_r <= 0;
        counter2_r <= 1;
        counter3_r <= 0-10;
    end
    else begin
        state_r <= state_w;
        counter_r <= counter_w;
        counter2_r <= counter2_w;
        counter3_r <= counter3_w;
        for(i=0;i<n_points;i=i+1) begin
            buffer_r_r[i] <= buffer_r_w[i];
            buffer_i_r[i] <= buffer_i_w[i];
        end
    end
end

endmodule