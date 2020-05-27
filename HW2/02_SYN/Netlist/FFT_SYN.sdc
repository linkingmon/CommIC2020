###################################################################

# Created by write_sdc on Wed May 27 01:12:42 2020

###################################################################
set sdc_version 2.0

set_units -time ns -resistance kOhm -capacitance pF -voltage V -current mA
set_wire_load_mode top
set_max_area 0
set_load -pin_load 0.05 [get_ports out_valid]
set_load -pin_load 0.05 [get_ports {dout_r[15]}]
set_load -pin_load 0.05 [get_ports {dout_r[14]}]
set_load -pin_load 0.05 [get_ports {dout_r[13]}]
set_load -pin_load 0.05 [get_ports {dout_r[12]}]
set_load -pin_load 0.05 [get_ports {dout_r[11]}]
set_load -pin_load 0.05 [get_ports {dout_r[10]}]
set_load -pin_load 0.05 [get_ports {dout_r[9]}]
set_load -pin_load 0.05 [get_ports {dout_r[8]}]
set_load -pin_load 0.05 [get_ports {dout_r[7]}]
set_load -pin_load 0.05 [get_ports {dout_r[6]}]
set_load -pin_load 0.05 [get_ports {dout_r[5]}]
set_load -pin_load 0.05 [get_ports {dout_r[4]}]
set_load -pin_load 0.05 [get_ports {dout_r[3]}]
set_load -pin_load 0.05 [get_ports {dout_r[2]}]
set_load -pin_load 0.05 [get_ports {dout_r[1]}]
set_load -pin_load 0.05 [get_ports {dout_r[0]}]
set_load -pin_load 0.05 [get_ports {dout_i[15]}]
set_load -pin_load 0.05 [get_ports {dout_i[14]}]
set_load -pin_load 0.05 [get_ports {dout_i[13]}]
set_load -pin_load 0.05 [get_ports {dout_i[12]}]
set_load -pin_load 0.05 [get_ports {dout_i[11]}]
set_load -pin_load 0.05 [get_ports {dout_i[10]}]
set_load -pin_load 0.05 [get_ports {dout_i[9]}]
set_load -pin_load 0.05 [get_ports {dout_i[8]}]
set_load -pin_load 0.05 [get_ports {dout_i[7]}]
set_load -pin_load 0.05 [get_ports {dout_i[6]}]
set_load -pin_load 0.05 [get_ports {dout_i[5]}]
set_load -pin_load 0.05 [get_ports {dout_i[4]}]
set_load -pin_load 0.05 [get_ports {dout_i[3]}]
set_load -pin_load 0.05 [get_ports {dout_i[2]}]
set_load -pin_load 0.05 [get_ports {dout_i[1]}]
set_load -pin_load 0.05 [get_ports {dout_i[0]}]
set_ideal_network -no_propagate  [get_ports clk]
create_clock [get_ports clk]  -period 10  -waveform {0 5}
set_clock_uncertainty 0.1  [get_clocks clk]
set_input_delay -clock clk  5  [get_ports clk]
set_input_delay -clock clk  5  [get_ports rst_n]
set_input_delay -clock clk  5  [get_ports in_valid]
set_input_delay -clock clk  5  [get_ports {din_r[11]}]
set_input_delay -clock clk  5  [get_ports {din_r[10]}]
set_input_delay -clock clk  5  [get_ports {din_r[9]}]
set_input_delay -clock clk  5  [get_ports {din_r[8]}]
set_input_delay -clock clk  5  [get_ports {din_r[7]}]
set_input_delay -clock clk  5  [get_ports {din_r[6]}]
set_input_delay -clock clk  5  [get_ports {din_r[5]}]
set_input_delay -clock clk  5  [get_ports {din_r[4]}]
set_input_delay -clock clk  5  [get_ports {din_r[3]}]
set_input_delay -clock clk  5  [get_ports {din_r[2]}]
set_input_delay -clock clk  5  [get_ports {din_r[1]}]
set_input_delay -clock clk  5  [get_ports {din_r[0]}]
set_input_delay -clock clk  5  [get_ports {din_i[11]}]
set_input_delay -clock clk  5  [get_ports {din_i[10]}]
set_input_delay -clock clk  5  [get_ports {din_i[9]}]
set_input_delay -clock clk  5  [get_ports {din_i[8]}]
set_input_delay -clock clk  5  [get_ports {din_i[7]}]
set_input_delay -clock clk  5  [get_ports {din_i[6]}]
set_input_delay -clock clk  5  [get_ports {din_i[5]}]
set_input_delay -clock clk  5  [get_ports {din_i[4]}]
set_input_delay -clock clk  5  [get_ports {din_i[3]}]
set_input_delay -clock clk  5  [get_ports {din_i[2]}]
set_input_delay -clock clk  5  [get_ports {din_i[1]}]
set_input_delay -clock clk  5  [get_ports {din_i[0]}]
set_output_delay -clock clk  5  [get_ports out_valid]
set_output_delay -clock clk  5  [get_ports {dout_r[15]}]
set_output_delay -clock clk  5  [get_ports {dout_r[14]}]
set_output_delay -clock clk  5  [get_ports {dout_r[13]}]
set_output_delay -clock clk  5  [get_ports {dout_r[12]}]
set_output_delay -clock clk  5  [get_ports {dout_r[11]}]
set_output_delay -clock clk  5  [get_ports {dout_r[10]}]
set_output_delay -clock clk  5  [get_ports {dout_r[9]}]
set_output_delay -clock clk  5  [get_ports {dout_r[8]}]
set_output_delay -clock clk  5  [get_ports {dout_r[7]}]
set_output_delay -clock clk  5  [get_ports {dout_r[6]}]
set_output_delay -clock clk  5  [get_ports {dout_r[5]}]
set_output_delay -clock clk  5  [get_ports {dout_r[4]}]
set_output_delay -clock clk  5  [get_ports {dout_r[3]}]
set_output_delay -clock clk  5  [get_ports {dout_r[2]}]
set_output_delay -clock clk  5  [get_ports {dout_r[1]}]
set_output_delay -clock clk  5  [get_ports {dout_r[0]}]
set_output_delay -clock clk  5  [get_ports {dout_i[15]}]
set_output_delay -clock clk  5  [get_ports {dout_i[14]}]
set_output_delay -clock clk  5  [get_ports {dout_i[13]}]
set_output_delay -clock clk  5  [get_ports {dout_i[12]}]
set_output_delay -clock clk  5  [get_ports {dout_i[11]}]
set_output_delay -clock clk  5  [get_ports {dout_i[10]}]
set_output_delay -clock clk  5  [get_ports {dout_i[9]}]
set_output_delay -clock clk  5  [get_ports {dout_i[8]}]
set_output_delay -clock clk  5  [get_ports {dout_i[7]}]
set_output_delay -clock clk  5  [get_ports {dout_i[6]}]
set_output_delay -clock clk  5  [get_ports {dout_i[5]}]
set_output_delay -clock clk  5  [get_ports {dout_i[4]}]
set_output_delay -clock clk  5  [get_ports {dout_i[3]}]
set_output_delay -clock clk  5  [get_ports {dout_i[2]}]
set_output_delay -clock clk  5  [get_ports {dout_i[1]}]
set_output_delay -clock clk  5  [get_ports {dout_i[0]}]
set_drive 1  [get_ports clk]
set_drive 1  [get_ports rst_n]
set_drive 1  [get_ports in_valid]
set_drive 1  [get_ports {din_r[11]}]
set_drive 1  [get_ports {din_r[10]}]
set_drive 1  [get_ports {din_r[9]}]
set_drive 1  [get_ports {din_r[8]}]
set_drive 1  [get_ports {din_r[7]}]
set_drive 1  [get_ports {din_r[6]}]
set_drive 1  [get_ports {din_r[5]}]
set_drive 1  [get_ports {din_r[4]}]
set_drive 1  [get_ports {din_r[3]}]
set_drive 1  [get_ports {din_r[2]}]
set_drive 1  [get_ports {din_r[1]}]
set_drive 1  [get_ports {din_r[0]}]
set_drive 1  [get_ports {din_i[11]}]
set_drive 1  [get_ports {din_i[10]}]
set_drive 1  [get_ports {din_i[9]}]
set_drive 1  [get_ports {din_i[8]}]
set_drive 1  [get_ports {din_i[7]}]
set_drive 1  [get_ports {din_i[6]}]
set_drive 1  [get_ports {din_i[5]}]
set_drive 1  [get_ports {din_i[4]}]
set_drive 1  [get_ports {din_i[3]}]
set_drive 1  [get_ports {din_i[2]}]
set_drive 1  [get_ports {din_i[1]}]
set_drive 1  [get_ports {din_i[0]}]