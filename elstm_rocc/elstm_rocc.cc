#include "rocc.h"
#include "mmu.h"
#include <cstring>
#include "half.hpp"
#include <math.h>
using namespace half_float;

#define NPE 2
#define CHKH 4
#define CHKW 4
#define BLKH 2

class elstm_rocc_t : public rocc_t
{
 public:
  const char* name() { return "elstm_rocc"; }

  reg_t custom0(rocc_insn_t insn, reg_t xs1, reg_t xs2)
  {
    
    switch (insn.funct)
    {
      case 0:
        size_x = xs1; size_h = xs2;
        break;
      case 1:
        h_th_uint = xs1; size_ts = xs2;
        break;
      case 2:
        p_w = xs1; p_b = xs2;
        break;
      case 3:
        p_ucolumn = xs1; ngrp = xs2;
        break;
      case 4:
        p_x = xs1; p_h = xs2;
        break;
      case 5:
        p_u = xs1; dram_trace_flag = xs2;
        compute();
        break;

      default:
        illegal_instruction();
        break;
    }
  }

  elstm_rocc_t()
  {
    printf("Construted the ELSTM ROCC co-processor.\n");
  }

 private:

  // control register
  reg_t size_x, size_h, size_ts;
  reg_t p_w, p_b, p_ucolumn, p_u;
  reg_t ngrp;
  reg_t p_x, p_h;
  reg_t h_th_uint;
  float *h_th_float;
  reg_t dram_trace_flag;

  uint clock;
  uint w_clock = 0;
  uint u_clock = 0;  

  // on-chip buffer
  float *buf_x, *buf_rw, *buf_ru, *buf_c, *buf_h;
  half *buf_b;
  uint16_t *buf_ucolumn;

  uint16_t h_total =0;
  uint16_t h_puring =0;
  uint16_t h_total_batch =0;
  uint16_t h_puring_batch =0;  


  // unpack the 64-bit ROCC interface and obtain the certain data
  uint32_t unpack32_rocc(reg_t rocc_data, uint n){
    return static_cast<uint32_t>(rocc_data>>(n*32) & 0xFFFFFFFF);
  }

  uint16_t unpack16_rocc(reg_t rocc_data, uint n){
    return static_cast<uint16_t>(rocc_data>>(n*16) & 0xFFFF);
  }

  uint8_t unpack8_rocc(reg_t rocc_data, uint n){
    return static_cast<uint8_t>(rocc_data>>(n*8) & 0xFF);
  }

  // initialize the storage
  // void init_network(reg_t size_x, reg_t size_h, reg_t size_w, reg_t size_u, reg_t p_b, reg_t p_ucolumn)
  void init_network()
  {
    // intialize the on-chip buffer
    buf_x = (float*) malloc(NPE*ngrp*size_x*sizeof(float));
    memset(buf_x, 0, NPE*ngrp*size_x*sizeof(float)); 
    buf_b = (half*) malloc(size_h*4*sizeof(half));     
    memset(buf_b, 0, size_h*4*sizeof(half));
    buf_ucolumn = (uint16_t*) malloc(size_h/4*sizeof(uint16_t));
    memset(buf_ucolumn, 0, size_h/4*sizeof(uint16_t));
    buf_h = (float*) malloc(size_h*sizeof(float));
    memset(buf_h, 0, size_h*sizeof(float));
    buf_c = (float*) malloc(size_h*sizeof(float));
    memset(buf_c, 0, size_h*sizeof(float));    
    buf_rw = (float*) malloc(NPE*ngrp*size_h*4*sizeof(float));
    memset(buf_rw, 0, NPE*ngrp*size_h*4*sizeof(float));
    buf_ru = (float*) malloc(size_h*4*sizeof(float));
    memset(buf_ru, 0, size_h*4*sizeof(float));
  }

  void load_ucolumn(){
    reg_t rocc_data;
    reg_t idx_ucolumn = 0;
    for (uint ii=0; ii<(uint((size_h/CHKW)/4)+1); ii++){
      clock++;
      rocc_data = p->get_mmu()->load_uint64(p_ucolumn+ii*sizeof(reg_t));
      uint cnt_rocc_word = 0;
      while ((idx_ucolumn < size_h) && (cnt_rocc_word <4)) {
        *(buf_ucolumn+idx_ucolumn) = unpack16_rocc(rocc_data, cnt_rocc_word);
        idx_ucolumn++;
        cnt_rocc_word++;
      }
    }
  }

  // load the bias
  void load_bias(){
    reg_t rocc_data;
    // number of bias is equal to size*4 and occupies size rocc words
    for (uint ii=0; ii<size_h; ii++){
      rocc_data = p->get_mmu()->load_uint64(p_b+ii*sizeof(reg_t));
      clock++;
      for (uint jj=0; jj<4; jj++){
        uint16_t rocc_unpack = unpack16_rocc(rocc_data,jj);
        half* tmp = reinterpret_cast<half*>(&rocc_unpack);
        *(buf_b+ii*4+jj) = *tmp;
      }
    }
  }

  // load input
  void load_x(uint idx_batch){
    uint size_batch = ngrp * NPE * size_x;
    reg_t rocc_data;
    // printf("idx_batch=%d\n", idx_batch);
    for (uint ii=0; ii<size_batch; ii++){
      if (ii%2==0){
        clock++;
        rocc_data = p->get_mmu()->load_uint64(p_x+(idx_batch*size_batch/2+(ii/2))*sizeof(reg_t));
        if (dram_trace_flag)
          printf("0,RD,0x%x\n",p_x+(idx_batch*size_batch/2+(ii/2))*sizeof(reg_t));
      }
      // printf("%x\n", p_x+idx_batch*size_batch/2+(ii/2)*sizeof(reg_t));
      // printf("rocc_data=%lx\n", rocc_data);
      uint32_t rocc_unpack = unpack32_rocc(rocc_data,ii%2);
      // printf("rocc_unpack=%x\n", rocc_unpack);
      float *tmp = reinterpret_cast<float*>(&rocc_unpack);
      // printf("cast=%f\n", *tmp);
      *(buf_x+ii) = *tmp;
    }
  }

  // store output (h)
  void store_h(){

    // for (uint ii=0; ii<size_h; ii+=2){
    //   uint32_t* low32 = reinterpret_cast<uint32_t*>((buf_h+ii));
    //   uint32_t* high32 = reinterpret_cast<uint32_t*>((buf_h+ii+1));
    //   reg_t tmp = (*high32 << 32) | *low32;      
    //   p->get_mmu()->store_uint64(p_h + ii/2*sizeof(reg_t), tmp);
    // }
    for (uint ii=0; ii<size_h; ii++){
      uint32_t *tmp = reinterpret_cast<uint32_t*>((buf_h+ii));
      p->get_mmu()->store_uint32(p_h + ii*sizeof(uint32_t), *tmp);
    }
  }

  void compute(){
    h_total =0;
    h_puring =0;

    clock = 0;
    w_clock = 0;
    u_clock = 0;

    uint size_batch = ngrp * NPE;
    uint idx_iter = 0;
    h_th_float = reinterpret_cast<float*>(&h_th_uint);
    // printf("h_th_float=%f\n", *h_th_float);
    init_network();
    load_bias();
    load_ucolumn();
    while (idx_iter < size_ts){
      load_x(idx_iter/size_batch);    
      batch_cal((size_ts-idx_iter)>size_batch? size_batch : (size_ts-idx_iter));
      idx_iter += size_batch;

      // printf("idx_iter = %d\n", idx_iter);
      // for (int ii=0; ii<size_h*4*NPE*ngrp; ii++){
      //   printf("buf_rw[%d]=%f\n", ii, *(buf_rw+ii));
      // }

      // NOTE: reset the memory
      memset(buf_rw, 0, NPE*ngrp*size_h*4*sizeof(float));
    }
    // add the tail u_clock
    clock += u_clock;

    // printf("time consumption = %d clocks\n", clock);
    printf("%d\n", clock);

    store_h();

    // for (int ii=0; ii<size_h; ii++){
    //   printf("buf_h[%d]=%f\n", ii, *(buf_h+ii));
    // }

    // printf("h_sparsity=%f\n", float(h_puring)/float(h_total));
    // printf("h_puring=%d, h_total=%d\n", h_puring, h_total);
  }

  // calcualte a batch
  void batch_cal(uint size_h_iter){

    // the matrix height is 4x the size_h
    uint h_blk = size_h*4 / (CHKH*BLKH);
    uint w_chk_wx = size_x / CHKW;
    uint w_chk_uh = size_h / CHKW;

    reg_t w_ptr;
    float* x_ptr;
    float* res_ptr;

    float f_t, i_t, c_t, o_t;

    // weight pointer = p_w
    w_ptr = p_w;
    x_ptr = buf_x;
    res_ptr = buf_rw;

    // printf("wptr=%x, xptr=%x\n", w_ptr, x_ptr);
    
    // compute all block columns for Wx (SpMV)
    for (uint idx_chkcol=0; idx_chkcol<w_chk_wx; idx_chkcol++){
      // obtain the input x vector which is needed by the block column computation in SpMV 
      float* in_ptr = x_ptr+idx_chkcol*CHKW;
      // compute block column and obtain the return value of dynamic w_ptr
      w_ptr = blkcol_cal(in_ptr, h_blk, res_ptr, w_ptr, NPE*ngrp, &w_clock);
    }

    // accumulate the time consumption
    uint wu = (w_clock*(ngrp-1) > u_clock) ? w_clock*ngrp : (u_clock - w_clock*(ngrp-1))+w_clock*ngrp;
    // printf("wu_clock=%d\n", wu);
    clock = clock + wu;
    // clear 
    w_clock = 0;
    u_clock = 0;    

    // compute SpM-SpV in Uh
    // initialize the address for the Uh of next time step
    w_ptr = p_u;
    x_ptr = buf_h;
    res_ptr = buf_ru;

    // printf("wptr=%x, xptr=%x\n", w_ptr, x_ptr);
    
    for (uint idx_batch=0; idx_batch<size_h_iter; idx_batch++){
        // reset the w_ptr
        w_ptr = p_u;      
      // compute one Uh time step
      for (uint idx_chkcol=0; idx_chkcol<w_chk_uh; idx_chkcol++){
        float* in_ptr = x_ptr+idx_chkcol*CHKW;

        //TODO: judge whether to compuate the block column with the values in buf_h;
        if (blk_puring(in_ptr)){
          w_ptr = blkcol_cal(in_ptr, h_blk, res_ptr, w_ptr, 1, &u_clock);
        } else{
          // skip the weight block column
          w_ptr += *(buf_ucolumn+idx_chkcol) * sizeof(reg_t);
        }
      }

      // for (int ii=0; ii<size_h*4; ii++){
      //   printf("buf_ru[%d]=%f\n", ii, *(buf_ru+ii));
      // }

      // vector operation to generate the new hidden state value
      float* buf_rw_batch = buf_rw + (idx_batch * size_h * 4);

      for (uint ii=0; ii<size_h; ii++){
        f_t = sigmod(*(buf_rw_batch+ii)+*(buf_ru+ii)+*(buf_b+ii));
        i_t = sigmod(*(buf_rw_batch+size_h+ii)+*(buf_ru+size_h+ii)+*(buf_b+size_h+ii));
        c_t = f_t * (*(buf_c+ii)) + i_t * tanh(*(buf_rw_batch+2*size_h+ii)+*(buf_ru+2*size_h+ii)+*(buf_b+2*size_h+ii));
        // store the c_t
        *(buf_c+ii) = c_t;
        o_t = sigmod(*(buf_rw_batch+3*size_h+ii)+*(buf_ru+3*size_h+ii)+*(buf_b+3*size_h+ii));
        // store the h_t
        *(buf_h+ii) = o_t * tanh(c_t);
      }

      // for (int ii=0; ii<size_h; ii++){
      //   printf("buf_c[%d]=%f\n", ii, *(buf_c+ii));
      // }      

      // reset the result buffer for Uh
      memset(buf_ru, 0, size_h*4*sizeof(float));
    }

    // printf("h_total_batch=%d, h_puring_batch=%d\n", h_total_batch, h_puring_batch);
    h_total_batch = 0;
    h_puring_batch = 0;

  }
  


  // process a block column in SpMV (including batch processing)
  reg_t blkcol_cal(float* in_ptr, uint h_blk, float* res_ptr, reg_t w_ptr, uint size_batch, uint* wu_clock){

    for (uint idx_blk=0; idx_blk<h_blk; idx_blk++){

      // obtain the intermediate results needed by the current block
      float* vec_res = res_ptr+idx_blk*CHKH*BLKH;
      // memcpy(vec_res, res_ptr+idx_blk*CHKH*BLKH, sizeof(float)*CHKH*BLKH);

      // process one block
      reg_t rocc_data = p->get_mmu()->load_uint64(w_ptr);
      *wu_clock =  *wu_clock + 1;
      if (dram_trace_flag)
        printf("0,RD,0x%x\n", w_ptr);
      // increase the pointer of weight load
      w_ptr +=sizeof(reg_t);
      uint32_t blk_header[2];
      blk_header[0] = unpack32_rocc(rocc_data,0);
      blk_header[1] = unpack32_rocc(rocc_data,1);

      for (uint idx_chk=0; idx_chk<BLKH; idx_chk++){

        // process one chunk
        uint8_t idxrow[CHKH];
        uint8_t eidxcol[CHKH];
        uint8_t chkw = blk_header[idx_chk] & 0x07;

        // printf("%d\n", chkw);        

        // obtain the col and row index
        for (uint ii=0; ii<CHKH; ii++){
          eidxcol[ii] = blk_header[idx_chk] >> (3+3*ii) & 0x07;
          idxrow[ii] = blk_header[idx_chk] >> (3+12+3*ii) & 0x07;
        }

        // compute
        for (uint ii=0; ii<chkw; ii++){
          // obtain the weight values
          // printf("%x\n", w_ptr);
          rocc_data = p->get_mmu()->load_uint64(w_ptr);
          *wu_clock =  *wu_clock + 1;
          if (dram_trace_flag)
            printf("0,RD,0x%x\n", w_ptr);          
          w_ptr +=sizeof(reg_t);

          for (uint jj=0; jj<CHKH; jj++){
            // printf("%lx\n", rocc_data);
            uint16_t rocc_unpack = unpack16_rocc(rocc_data,jj);
            // printf("%x\n", rocc_unpack);
            half* cur_w = reinterpret_cast<half*>(&rocc_unpack);
            // std::cout << *cur_w << std::endl;
            // batch processing
            for (uint idx_batch=0; idx_batch<size_batch; idx_batch++){
              *(vec_res+idx_batch*size_h*4+idxrow[jj]) = *(vec_res+idx_batch*size_h*4+idxrow[jj]) + (*cur_w) * dec_x(in_ptr+idx_batch*size_x, eidxcol[jj], ii);
            }
          }
        }
      }
    }
    return w_ptr;   
  }

  float dec_x(float* in_ptr, uint8_t eidx, uint8_t idxcol){
    // define the decoder book
    uint decbook[4][8]={0};
    decbook[0][4] = 1;
    decbook[0][5] = 1;
    decbook[0][6] = 2;    
    decbook[0][7] = 3;    
    decbook[1][0] = 1;
    decbook[1][1] = 1;    
    decbook[1][2] = 2;
    decbook[1][4] = 2;
    decbook[1][3] = 3;
    decbook[1][5] = 3;
    decbook[1][6] = 3;        
    decbook[2][0] = 2;
    decbook[2][1] = 3;
    decbook[2][2] = 3;
    decbook[2][4] = 3;
    decbook[3][0] = 3;

    // printf("eidx=%d,idxcol=%d\n", eidx, idxcol);

    // printf("decx=%f\n", *(in_ptr+decbook[idxcol][eidx]));

    return *(in_ptr+decbook[idxcol][eidx]);
  }

  // sigmode function
  inline float sigmod(float x){
    return 1/(1+exp(-x));
  }

  // 
  inline bool blk_puring(float* x){
    h_total++;
    h_total_batch++;
    float h_th = *h_th_float;
    // printf("th=%f\n", h_th);
    // printf("%f,%f,%f,%f\n", fabs(*x), fabs(*(x+1)), fabs(*(x+2)), fabs(*(x+3)));
    if (((fabs(*x) + fabs(*(x+1)) + fabs(*(x+2)) + fabs(*(x+3))) >= h_th) || (fabs(*x) >= h_th/2) || (fabs(*(x+1)) >= h_th/2) || (fabs(*(x+2)) >= h_th/2) || (fabs(*(x+3)) >= h_th/2)){
      // printf("TRUE\n");
      return true;
    }
    // printf("FALSE\n");
    // printf("%f,%f,%f,%f\n", fabs(*x), fabs(*(x+1)), fabs(*(x+2)), fabs(*(x+3)));
    h_puring++;
    h_puring_batch++;
    return false;
  }

};

REGISTER_EXTENSION(elstm_rocc, []() { return new elstm_rocc_t; })
