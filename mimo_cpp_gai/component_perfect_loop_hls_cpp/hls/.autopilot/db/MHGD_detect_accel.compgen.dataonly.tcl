# This script segment is generated automatically by AutoPilot

set axilite_register_dict [dict create]
set port_control {
x_hat { 
	dir I
	width 64
	depth 1
	mode ap_none
	offset 16
	offset_end 27
}
H { 
	dir I
	width 64
	depth 1
	mode ap_none
	offset 28
	offset_end 39
}
y { 
	dir I
	width 64
	depth 1
	mode ap_none
	offset 40
	offset_end 51
}
v_tb { 
	dir I
	width 64
	depth 1
	mode ap_none
	offset 52
	offset_end 63
}
}
dict set axilite_register_dict control $port_control


