//##############################################################################
//
// lb3d_prime.c
//
//  - Lattice Boltzmann
//
//  - D3Q19
//

#include "lb3d_prime.h"

int main( int argc, char **argv)
{
  int n;
  double k;
  int time, frames;

  struct lattice_struct *lattice;

  setbuf( stdout, (char*)NULL); // Don't buffer screen output.

#if VERBOSITY_LEVEL > 0
  printf("%s %d >> lb3d_prime.c: main() -- Hello.\n", __FILE__, __LINE__);
#endif /* VERBOSITY_LEVEL > 0 */

  construct_lattice( &lattice, argc, argv);

  init_problem( lattice);

  output_frame( lattice);
  
  for( frames = 0, time=1; time<=lattice->NumTimeSteps; time++)
  {
    lattice->time = time;

    stream( lattice);

    if( !lattice->param.GZL && (!lattice->param.AllBoundaryPeriodic)) bcs( lattice);

    compute_macro_vars( lattice);

#if SAVE_MEMO
//feq was calculated in collide_save(lattice)
#else    
    compute_feq( lattice);
#endif

    collide( lattice);

    if( lattice->param.GZL && (!lattice->param.AllBoundaryPeriodic)) bcs( lattice);
//compute_feq( lattice);

    if( !(time%lattice->param.FrameRate)) { output_frame( lattice); frames++;}

  } /* for( time=1; time<=lattice->NumTimeSteps; time++) */

  destruct_lattice( lattice);

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf("%s %d >> lb3d_prime.c: main() -- Terminating normally.\n",
      __FILE__, __LINE__);
  printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */

  return 0;

} /* int main( int argc, char **argv) */
