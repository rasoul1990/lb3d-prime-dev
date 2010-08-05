#ifndef MCMP_PRIME_H
#define MCMP_PRIME_H
//##############################################################################
//
// lb3d_prime.h
//
//  - Header file for lb3d_prime.
//

#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef CLK_TCK
// Old versions of gcc use CLK_TCK.
#define CLK_TCK CLOCKS_PER_SEC
#endif

#include "flags.h"
#include "bc_flags.h"
#include "ic_flags.h"
#include "process.h"
#include "lattice.h"
#include "forward_declarations.h"
#include "params.h"

// D3Q19
//
//             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//             |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
int vx[19] = { 0,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0};
int vy[19] = { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
int vz[19] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};
//             C  W  E  N  S  T  B  N  N  S  S  T  T  B  B  T  T  B  B
//                                  W  E  W  E  W  E  W  E  N  S  N  S

#define C  0
#define W  1
#define E  2
#define N  3
#define S  4
#define T  5
#define B  6
#define NW 7
#define NE 8
#define SW 9
#define SE 10
#define TW 11
#define TE 12
#define BW 13
#define BE 14
#define TN 15
#define TS 16
#define BN 17
#define BS 18

#define EPS .0000000000001
#define PI  3.1415926535

#define HRULE0 "- - - - - - - - - - - - - - - - - - - - " \
               "- - - - - - - - - - - - - - - - - - - - "

#define HRULE1 "----------------------------------------" \
               "----------------------------------------"

#define HRULE2 "========================================" \
               "========================================"

#if DO_NOT_STORE_SOLIDS
// TODO: Incorporate version that omits storage of interior solids (which
// are not involved in flow).
//#include "min_nodes/compute.c"
//#include "min_nodes/stream.c"
//#include "min_nodes/bcs.c"
//#include "min_nodes/collide.c"
//#include "min_nodes/lbio.c"
//#include "min_nodes/latman.c"
#else /* !( DO_NOT_STORE_SOLIDS) */
#include "lattice.c"
#include "process.c"
#include "compute.c"
#include "stream.c"
#include "collide.c"
#include "bcs.c"
#include "lbio.c"
#include "latman.c"
#endif /* DO_NOT_STORE_SOLIDS */

#if 0
void report_flags()
{
  struct report_struct flags;
  report_open(          &flags, "./out/flags");
  report_integer_entry( &flags, "VERBOSITY_LEVEL", VERBOSITY_LEVEL, "");
  report_integer_entry( &flags, "SAY_HI", SAY_HI, "");
  report_integer_entry( &flags, "NUM_FLUID_COMPONENTS", NUM_FLUID_COMPONENTS, "");
  report_integer_entry( &flags, "INAMURO_SIGMA_COMPONENT", INAMURO_SIGMA_COMPONENT, "");
  report_integer_entry( &flags, "ZHANG_AND_CHEN_ENERGY_TRANSPORT", ZHANG_AND_CHEN_ENERGY_TRANSPORT, "");
  report_integer_entry( &flags, "POROUS_MEDIA", POROUS_MEDIA, "");
  report_integer_entry( &flags, "STORE_UEQ", STORE_UEQ, "");
  report_integer_entry( &flags, "DO_NOT_STORE_SOLIDS", DO_NOT_STORE_SOLIDS, "");
  report_integer_entry( &flags, "NON_LOCAL_FORCES", NON_LOCAL_FORCES, "");
  report_integer_entry( &flags, "MANAGE_BODY_FORCE", MANAGE_BODY_FORCE, "");
  report_integer_entry( &flags, "STORE_BTC", STORE_BTC, "");
  report_integer_entry( &flags, "DETERMINE_FLOW_DIRECTION", DETERMINE_FLOW_DIRECTION, "");
  report_integer_entry( &flags, "BC_SOLID_NODE", BC_SOLID_NODE, "");
  report_integer_entry( &flags, "BC_FLUID_NODE", BC_FLUID_NODE, "");
  report_integer_entry( &flags, "BC_SLIP_NODE", BC_SLIP_NODE, "");
  report_integer_entry( &flags, "BC_FILM_NODE", BC_FILM_NODE, "");
  report_integer_entry( &flags, "BC_PRESSURE_N_IN", BC_PRESSURE_N_IN, "");
  report_integer_entry( &flags, "BC_PRESSURE_S_IN", BC_PRESSURE_S_IN, "");
  report_integer_entry( &flags, "BC_PRESSURE_E_IN", BC_PRESSURE_E_IN, "");
  report_integer_entry( &flags, "BC_PRESSURE_W_IN", BC_PRESSURE_W_IN, "");
  report_integer_entry( &flags, "BC_PRESSURE_N_OUT", BC_PRESSURE_N_OUT, "");
  report_integer_entry( &flags, "BC_PRESSURE_S_OUT", BC_PRESSURE_S_OUT, "");
  report_integer_entry( &flags, "BC_PRESSURE_E_OUT", BC_PRESSURE_E_OUT, "");
  report_integer_entry( &flags, "BC_PRESSURE_W_OUT", BC_PRESSURE_W_OUT, "");
  report_integer_entry( &flags, "BC_VELOCITY_N_IN", BC_VELOCITY_N_IN, "");
  report_integer_entry( &flags, "BC_VELOCITY_S_IN", BC_VELOCITY_S_IN, "");
  report_integer_entry( &flags, "BC_VELOCITY_E_IN", BC_VELOCITY_E_IN, "");
  report_integer_entry( &flags, "BC_VELOCITY_W_IN", BC_VELOCITY_W_IN, "");
  report_integer_entry( &flags, "BC_VELOCITY_N_OUT", BC_VELOCITY_N_OUT, "");
  report_integer_entry( &flags, "BC_VELOCITY_S_OUT", BC_VELOCITY_S_OUT, "");
  report_integer_entry( &flags, "BC_VELOCITY_E_OUT", BC_VELOCITY_E_OUT, "");
  report_integer_entry( &flags, "BC_VELOCITY_W_OUT", BC_VELOCITY_W_OUT, "");
  report_integer_entry( &flags, "WRITE_MACRO_VAR_DAT_FILES", WRITE_MACRO_VAR_DAT_FILES, "");
  report_integer_entry( &flags, "WRITE_RHO_AND_U_TO_TXT", WRITE_RHO_AND_U_TO_TXT, "");
  report_integer_entry( &flags, "WRITE_PDF_DAT_FILES", WRITE_PDF_DAT_FILES, "");
  report_integer_entry( &flags, "WRITE_PDF_TO_TXT", WRITE_PDF_TO_TXT, "");
  report_integer_entry( &flags, "INACTIVE_NODE", INACTIVE_NODE, "");
  report_integer_entry( &flags, "PUKE_NEGATIVE_DENSITIES", PUKE_NEGATIVE_DENSITIES, "");
  report_integer_entry( &flags, "SOLID_COLOR_IS_CHECKERBOARD", SOLID_COLOR_IS_CHECKERBOARD, "");
  report_integer_entry( &flags, "SOLID_COLOR_IS_BLACK", SOLID_COLOR_IS_BLACK, "");
  report_integer_entry( &flags, "DELAY", DELAY, "");
  report_integer_entry( &flags, "END_GRAV", END_GRAV, "");
  report_integer_entry( &flags, "MARK_ORIGIN_FOR_REFERENCE", MARK_ORIGIN_FOR_REFERENCE, "");
  report_integer_entry( &flags, "PERTURBATIONS", PERTURBATIONS, "");
  report_integer_entry( &flags, "IC_UNIFORM_RHO_A", IC_UNIFORM_RHO_A, "");
  report_integer_entry( &flags, "IC_UNIFORM_RHO_B", IC_UNIFORM_RHO_B, "");
  report_integer_entry( &flags, "IC_UNIFORM_RHO_IN", IC_UNIFORM_RHO_IN, "");
  report_integer_entry( &flags, "IC_BUBBLE", IC_BUBBLE, "");
  report_integer_entry( &flags, "IC_DIAGONAL", IC_DIAGONAL, "");
  report_integer_entry( &flags, "IC_2X2_CHECKERS", IC_2X2_CHECKERS, "");
  report_integer_entry( &flags, "IC_STATIC", IC_STATIC, "");
  report_integer_entry( &flags, "IC_RECTANGLE", IC_RECTANGLE, "");
  report_integer_entry( &flags, "IC_DOT", IC_DOT, "");
  report_integer_entry( &flags, "IC_WOLF_GLADROW_DIFFUSION", IC_WOLF_GLADROW_DIFFUSION, "");
  report_integer_entry( &flags, "IC_YIN_YANG", IC_YIN_YANG, "");
  report_integer_entry( &flags, "IC_HYDROSTATIC", IC_HYDROSTATIC, "");
  report_close(         &flags);
}
#endif

#endif /* MCMP_PRIME_H */
