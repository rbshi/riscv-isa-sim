#include "rocc.h"
#include "mmu.h"
#include <cstring>
#include "half.hpp"
using half_float::half;

class gemv_rocc_t : public rocc_t
{
 public:
  const char* name() { return "gemv_rocc"; }

  reg_t custom0(rocc_insn_t insn, reg_t xs1, reg_t xs2)
  {
    
    // printf("instruction OPCODE=%d\n", insn.opcode);
    // printf("instruction RD=%d\n", insn.rd);
    // printf("instruction XS2=%d\n", insn.xs2);
    // printf("instruction XS1=%d\n", insn.xs1);
    // printf("instruction XD=%d\n", insn.xd);            
    // printf("instruction RS1=%d\n", insn.rs1);
    // printf("instruction RS2=%d\n", insn.rs2);
    // printf("instruction FUNCT=%d\n", insn.funct);            

    // printf("instruction xs1=%d\n", xs1);
    // printf("instruction xs2=%d\n", xs2);    


    switch (insn.funct)
    {
      case 0:
        mat_w = xs1; mat_h = xs2;
        break;
      case 1:
        mat_addr = xs1; vec_addr = xs2;
        break;
      case 2:
        res_addr = xs1;
        gemv(mat_w, mat_h, mat_addr, vec_addr, res_addr);
        break;
      default:
        illegal_instruction();
        break;
    }

  }

  gemv_rocc_t()
  {
    printf("Construted the GEMV ROCC co-processor.\n");
    half a(3.4), b(5);
    half c = a * b;
    c += 3;
    if(c > a)
        std::cout << c << std::endl;
  }

 private:
  reg_t mat_w, mat_h;
  reg_t mat_addr, vec_addr;
  reg_t res_addr;


  void gemv(reg_t mat_w, reg_t mat_h, reg_t mat_addr, reg_t vec_addr, reg_t res_addr){
    reg_t tmp;
    for (uint ii=0; ii<mat_h; ii++){
      tmp = 0;
      for (uint jj=0; jj<mat_w; jj++){
        tmp = tmp + p->get_mmu()->load_uint64(mat_addr+(ii*mat_w+jj)*sizeof(reg_t)) * p->get_mmu()->load_uint64(vec_addr+jj*sizeof(reg_t));
    }
    p->get_mmu()->store_uint64(res_addr+ii*sizeof(reg_t), tmp);
  }
}
};

REGISTER_EXTENSION(gemv_rocc, []() { return new gemv_rocc_t; })
