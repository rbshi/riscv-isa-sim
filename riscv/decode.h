#ifndef _RISCV_DECODE_H
#define _RISCV_DECODE_H

#define __STDC_LIMIT_MACROS
#include <stdint.h>
typedef int int128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));

#define support_64bit 1
typedef int64_t sreg_t;
typedef uint64_t reg_t;
typedef uint64_t freg_t;

const int OPCODE_BITS = 7;
const int JTYPE_OPCODE_BITS = 5;

const int GPR_BITS = 8*sizeof(reg_t);
const int GPRID_BITS = 5;
const int NGPR = 1 << GPRID_BITS;

const int FPR_BITS = 64;
const int FPRID_BITS = 5;
const int NFPR = 1 << FPRID_BITS;

const int IMM_BITS = 12;
const int TARGET_BITS = 27;
const int SHAMT_BITS = 6;
const int FUNCT_BITS = 3;
const int FFUNCT_BITS = 5;
const int BIGIMM_BITS = 20;
const int BRANCH_ALIGN_BITS = 1;
const int JUMP_ALIGN_BITS = 1;

#define SR_ET    0x0000000000000001ULL
#define SR_PS    0x0000000000000004ULL
#define SR_S     0x0000000000000008ULL
#define SR_EF    0x0000000000000010ULL
#define SR_UX    0x0000000000000020ULL
#define SR_SX    0x0000000000000040ULL
#define SR_IM    0x000000000000FF00ULL
#define SR_ZERO  ~(SR_ET | SR_PS | SR_S | SR_EF | SR_UX | SR_SX | SR_IM)
#define SR_IM_SHIFT 8
#define TIMER_IRQ 7

#define FP_RD_NE  0
#define FP_RD_0   1
#define FP_RD_DN  2
#define FP_RD_UP  3
#define FP_RD_NMM 4
#define FSR_RD_SHIFT 10
#define FSR_RD   (0x7 << FSR_RD_SHIFT)

#define FPEXC_NX 0x01
#define FPEXC_UF 0x02
#define FPEXC_OF 0x04
#define FPEXC_DZ 0x02
#define FPEXC_NV 0x10

#define FSR_CEXC_SHIFT 5
#define FSR_NVC  (FPEXC_NV << FSR_CEXC_SHIFT)
#define FSR_OFC  (FPEXC_OF << FSR_CEXC_SHIFT)
#define FSR_UFC  (FPEXC_UF << FSR_CEXC_SHIFT)
#define FSR_DZC  (FPEXC_DZ << FSR_CEXC_SHIFT)
#define FSR_NXC  (FPEXC_NX << FSR_CEXC_SHIFT)
#define FSR_CEXC (FSR_NVC | FSR_OFC | FSR_UFC | FSR_DZC | FSR_NXC)

#define FSR_AEXC_SHIFT 0
#define FSR_NVA  (FPEXC_NV << FSR_AEXC_SHIFT)
#define FSR_OFA  (FPEXC_OF << FSR_AEXC_SHIFT)
#define FSR_UFA  (FPEXC_UF << FSR_AEXC_SHIFT)
#define FSR_DZA  (FPEXC_DZ << FSR_AEXC_SHIFT)
#define FSR_NXA  (FPEXC_NX << FSR_AEXC_SHIFT)
#define FSR_AEXC (FSR_NVA | FSR_OFA | FSR_UFA | FSR_DZA | FSR_NXA)

#define FSR_ZERO ~(FSR_RD | FSR_AEXC | FSR_CEXC)

// note: bit fields are in little-endian order
struct itype_t
{
  unsigned imm : IMM_BITS;
  unsigned funct : FUNCT_BITS;
  unsigned rb : GPRID_BITS;
  unsigned ra : GPRID_BITS;
  unsigned opcode : OPCODE_BITS;
};

struct jtype_t
{
  unsigned target : TARGET_BITS;
  unsigned jump_opcode : JTYPE_OPCODE_BITS;
};

struct rtype_t
{
  unsigned rc : GPRID_BITS;
  unsigned shamt : SHAMT_BITS;
  unsigned unused : 1;
  unsigned funct : FUNCT_BITS;
  unsigned rb : GPRID_BITS;
  unsigned ra : GPRID_BITS;
  unsigned opcode : OPCODE_BITS;
};

struct btype_t
{
  unsigned bigimm : BIGIMM_BITS;
  unsigned rt : GPRID_BITS;
  unsigned opcode : OPCODE_BITS;
};

struct ftype_t
{
  unsigned rc : FPRID_BITS;
  unsigned rd : FPRID_BITS;
  unsigned ffunct : FFUNCT_BITS;
  unsigned rb : FPRID_BITS;
  unsigned ra : FPRID_BITS;
  unsigned opcode : OPCODE_BITS;
};

union insn_t
{
  itype_t itype;
  jtype_t jtype;
  rtype_t rtype;
  btype_t btype;
  ftype_t ftype;
  uint32_t bits;
};

// helpful macros, etc
#define RA R[insn.rtype.ra]
#define RB R[insn.rtype.rb]
#define RC R[insn.rtype.rc]
#define FRA FR[insn.ftype.ra]
#define FRB FR[insn.ftype.rb]
#define FRC FR[insn.ftype.rc]
#define FRD FR[insn.ftype.rd]
#define BIGIMM insn.btype.bigimm
#define IMM insn.itype.imm
#define SIMM ((int32_t)((uint32_t)insn.itype.imm<<(32-IMM_BITS))>>(32-IMM_BITS))
#define SHAMT insn.rtype.shamt
#define TARGET insn.jtype.target
#define BRANCH_TARGET (pc + (SIMM << BRANCH_ALIGN_BITS))
#define JUMP_TARGET ((pc & ~((1<<(TARGET_BITS+JUMP_ALIGN_BITS))-1)) + (TARGET << JUMP_ALIGN_BITS))

#define require_supervisor if(!(sr & SR_S)) throw trap_privileged_instruction
#define require64 if(gprlen != 64) throw trap_illegal_instruction
#define require_fp if(!(sr & SR_EF)) throw trap_fp_disabled
#define cmp_trunc(reg) (reg_t(reg) << (64-gprlen))
#define set_fp_exceptions ({ set_fsr((fsr & ~FSR_CEXC) | \
                                (softfloat_exceptionFlags << FSR_AEXC_SHIFT) | \
                                (softfloat_exceptionFlags << FSR_CEXC_SHIFT)); \
                             softfloat_exceptionFlags = 0; })

static inline sreg_t sext32(int32_t arg)
{
  return arg;
}

#endif
