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

  process_tic( lattice);

  init_problem( lattice);

  output_frame( lattice);

  process_toc( lattice);
  display_etime( lattice);
  process_tic( lattice);

  for( frames = 0, time=1; time<=lattice->NumTimeSteps; time++)
  {
    lattice->time = time;

    stream( lattice);

    if( do_post_streaming_bcs(lattice)) { bcs( lattice);}

    compute_macro_vars( lattice);

#if SAVE_MEMO
    // feq was calculated in collide_save(lattice)
#else
    compute_feq( lattice);
#endif

    collide( lattice);

    if( do_post_collision_bcs(lattice)) { bcs( lattice);}

    //compute_feq( lattice);

    if( !(time%lattice->param.FrameRate))
    {
      output_frame( lattice); frames++;
      process_toc( lattice);
      display_etime( lattice);
    }

  } /* for( time=1; time<=lattice->NumTimeSteps; time++) */

  process_barrier();
  process_toc( lattice);
  display_etime( lattice);

  destruct_lattice( lattice);

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf("%s %d >> lb3d_prime.c: main() -- Terminating normally.\n",
      __FILE__, __LINE__);
  printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */

  return 0;

} /* int main( int argc, char **argv) */
