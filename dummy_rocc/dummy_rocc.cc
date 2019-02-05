#include "rocc.h"
#include "mmu.h"
#include <cstring>

class dummy_rocc_t : public rocc_t
{
 public:
  const char* name() { return "dummy_rocc"; }

  reg_t custom0(rocc_insn_t insn, reg_t xs1, reg_t xs2)
  {
    reg_t prev_acc;

    printf("instruction OPCODE=%d\n", insn.opcode);
    printf("instruction RD=%d\n", insn.rd);
    printf("instruction XS2=%d\n", insn.xs2);
    printf("instruction XS1=%d\n", insn.xs1);
    printf("instruction XD=%d\n", insn.xd);            
    printf("instruction RS1=%d\n", insn.rs1);
    printf("instruction RS2=%d\n", insn.rs2);
    printf("instruction FUNCT=%d\n", insn.funct);            

    printf("instruction xs1=%d\n", xs1);
    printf("instruction xs2=%d\n", xs2);    


    if (xs2 >= num_acc)
      illegal_instruction();

    switch (insn.funct)
    {
      case 0: // acc <- xs1
        acc[xs2] = xs1;
        break;
      case 1: // xd <- acc (the only real work is the return statement below)
        prev_acc = acc[xs2];
        break;
      case 2: // acc[rs2] <- Mem[xs1]
        acc[xs2] = p->get_mmu()->load_uint64(xs1);
        break;
      case 3: // acc[rs2] <- accX + xs1
        acc[xs2] += xs1;
        break;
      default:
        illegal_instruction();
        break;
    }

    return prev_acc; // in all cases, xd <- previous value of acc[rs2]
  }

  dummy_rocc_t()
  {
    memset(acc, 0, sizeof(acc));
  }

 private:
  static const int num_acc = 4;
  reg_t acc[num_acc];
};

REGISTER_EXTENSION(dummy_rocc, []() { return new dummy_rocc_t; })
