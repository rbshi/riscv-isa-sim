# Instruction of E-LSTM Co-processor

Function: 0
RS1: size of `x` in block
RS2: size of `h` in block
RD: -
Behavior:  load the size of `x`, unit: block (ROCC interface is 64 bits and each input block contains 4 scalars); load the hidden size (in block)

Function: 1
RS1: threshold for hidden state puring (h_th)
RS2: sizxe of time step (size_ts)
RD: 
Behavior: load the threshold in hidden state puring and the size of time step 

Function: 2
RS1: pointer of `W`
RS2: pointer of `b`
RD: 
Behavior: load the pointer of `W` and pointer of `b`

Function: 3
RS1: pointer of `U` column address
RS2: `N_grp`
RD: 
Behavior: load the pointer of `U` column address; load `N_grp`

Function: 4
RS1: pointer of `x`
RS2: pointer of `h`
RD: 
Behavior: load pointers of input `x` and output `h`

Function: 5
RS1: pointer of `u`
RS2: 
RD: 
Behavior: 

# Hardware

Buffer0: `x`
Buffer1: `b`
Buffer2: address offset of `U` columns
Buffer3: `Rw` result of Wx_t
Buffer4: `Ru` result of Uh_t
Buffer5: `f`
Buffer6: `i`
Buffer7: `c`
Buffer8: `h`


# Work-flow

1. Load network parameters
* `x` size (in block) and `h` size
* `weight W` size (after eSELL compression, in block) and time step
* address of compressed `W`
* address offset of `U` columns -> Buffer2
* store `b` -> Buffer1
* `N_grp`

2. Load the input and output address
* address of `x`
* address of `h`





