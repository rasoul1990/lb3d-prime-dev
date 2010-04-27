//##############################################################################
//
// get_r_and_dp.c
//
//  - Lattice Boltzmann
//
//  - Read in frames from a previously complete LB run with a bubble.
//
//  - Compute the radius of the bubble and the pressure difference
//    inside versus outside the bubble.
//

#include "lb3d_prime.h"

double input_frame( lattice_ptr lattice, char *dirname);
double get_radius_bubble( lattice_ptr lattice, double rho_cut);
double get_radius_drop( lattice_ptr lattice, double rho_cut);
double get_dp( lattice_ptr lattice);
void read_dat( 
       lattice_ptr lattice,
       double *a, 
       int     stride,
       char   *filename);

int main( int argc, char **argv)
{
  int n;
  double k;
  int time, frames;

  struct lattice_struct *lattice;

  char dirname_out[1024];

  setbuf( stdout, (char*)NULL); // Don't buffer screen output.

#if VERBOSITY_LEVEL > 0
  printf("%s %d >> get_r_and_dp.c: main() -- Hello.\n", __FILE__, __LINE__);
#endif /* VERBOSITY_LEVEL > 0 */

  construct_lattice( &lattice, argc, argv);

  if( argc == 2)
  {
    sprintf( dirname_out, "out_%s", argv[1]);
  }
  else
  {
    sprintf( dirname_out, "out");
  }

  lattice->time = 0;

  for( frames = 0; frames<=lattice->param.NumFrames; frames++)
  {
    lattice->time = frames * lattice->param.FrameRate;

    input_frame( lattice, dirname_out);

  //printf( "  Bubble: radius = %f\n", get_radius_bubble( lattice, 487.));
  //printf( "  Drop:   radius = %f\n", get_radius_drop( lattice, 487.));
    printf( "  Bubble: radius = %f\n", get_radius_bubble( lattice, 356.6148 ));
    printf( "  Drop:   radius = %f\n", get_radius_drop(   lattice, 356.6148 ));
    
    get_dp( lattice);

  } /* for( time=1; time<=lattice->NumTimeSteps; time++) */

  destruct_lattice( lattice);

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf("%s %d >> get_r_and_dp.c: main() -- Terminating normally.\n",
      __FILE__, __LINE__);
  printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */

  return 0;

} /* int main( int argc, char **argv) */

double input_frame( lattice_ptr lattice, char *dirname)
{
  char   filename[1024];
  FILE   *o, *o_u, *o_rho, *o_ux, *o_uy, *o_ueq, *o_ueq_x, *o_ueq_y;
  int    *node_ptr;
  int    i, n;
  double *macro_vars_ptr;
  double *ueq;
  int    frame;
#if WRITE_RHO_AND_U_TO_TXT
  int    i, j;
#endif /* WRITE_RHO_AND_U_TO_TXT */
  double min_u[3], max_u[3],  ave_u[3];
  double min_rho, max_rho,   ave_rho;
  double rho_ratio, u_x_ratio, u_y_ratio;
  int    subs;
  int    ni, nj, nk;
  int    njdiv2, nkdiv2;

  frame = (int)((double)lattice->time/(double)lattice->param.FrameRate);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./%s/rho_%dx%dx%d_frame%04d_subs%02d.dat", 
        dirname, ni, nj, nk, frame, subs);
    read_dat( 
      lattice,
      &(lattice->macro_vars[subs][0].rho),
      /*stride*/ 4,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */
    compute_max_rho( lattice, &max_rho, subs);
    compute_min_rho( lattice, &min_rho, subs);
    compute_ave_rho( lattice, &ave_rho, subs);
    printf("%s %d >> min_rho = %f, max_rho = %f, ave_rho = %f\n",
      __FILE__,__LINE__, min_rho, max_rho, ave_rho);

#if 0

#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_x_%dx%dx%d_frame%04d_subs%02d.dat", 
        ni, nj, nk, frame, subs);
    read_dat( 
      lattice,
      &(lattice->macro_vars[subs][0].u[0]),
      /*stride*/ 4,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_y_%dx%dx%d_frame%04d_subs%02d.dat", 
        ni, nj, nk, frame, subs);
    read_dat( 
      lattice,
      &(lattice->macro_vars[subs][0].u[1]),
      /*stride*/ 4,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_z_%dx%dx%d_frame%04d_subs%02d.dat", 
        ni, nj, nk, frame, subs);
    read_dat( 
      lattice,
      &(lattice->macro_vars[subs][0].u[2]),
      /*stride*/ 4,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */

    compute_max_u( lattice, max_u, subs);
    compute_min_u( lattice, min_u, subs);
    compute_ave_u( lattice, ave_u, subs);
    printf("%s %d >> min_u = [ %f %f %f], max_u = [ %f %f %f], ave_u = [ %f %f %f]\n",
      __FILE__,__LINE__,
    min_u[0], min_u[1], min_u[2],
    max_u[0], max_u[1], max_u[2], 
    ave_u[0], ave_u[1], ave_u[2] );
#endif

    njdiv2 = (int)floor((double)nj/2.);
    nkdiv2 = (int)floor((double)nk/2.);
    printf("%s %d >> slice_rho = [", __FILE__,__LINE__);
    for( i=0; i<ni; i++)
    {
      printf(" %f ", 
        lattice->macro_vars[subs][XYZ2N(i,njdiv2,nkdiv2,ni,nj)].rho);
    }
    printf("]\n");


  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

} /* double input_frame( lattice_ptr lattice, char *dirname) */

// double get_radius_bubble( lattice_ptr lattice, double rho_cut)
//##############################################################################
// Assumes sphereical region.
// Computes r = ((4/3)*volume/pi)^(1/3)
double get_radius_bubble( lattice_ptr lattice, double rho_cut)
{
  int i, j, k;
  int n, ni, nj, nk;

  int volume = 0;

  for( n=0; n<lattice->NumNodes; n++)
  {
    if( lattice->macro_vars[0][n].rho <= rho_cut)
    {
      volume++;
    }
  }

  return pow((3./4.)*(double)volume/3.141593,1./3.);

} /* double get_radius_bubble( lattice_ptr lattice) */
 
// double get_radius_drop( lattice_ptr lattice, double rho_cut)
//##############################################################################
// Assumes sphereical region.
// Computes r = ((4/3)*volume/pi)^(1/3)
double get_radius_drop( lattice_ptr lattice, double rho_cut)
{
  int i, j, k;
  int n, ni, nj, nk;

  int volume = 0;

  for( n=0; n<lattice->NumNodes; n++)
  {
    if( lattice->macro_vars[0][n].rho >= rho_cut)
    {
      volume++;
    }
  }

  return pow((3./4.)*(double)volume/3.141593,1./3.);

} /* double get_radius_drop( lattice_ptr lattice) */
    
double get_dp( lattice_ptr lattice)
{
} /* double get_dp( lattice_ptr lattice) */

void read_dat( 
       lattice_ptr lattice,
       double *a, 
       int     stride,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  double *Xpix1;
  FILE *infile;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(double);

  if( !(Xpix1 = ( double*)malloc( size)))
  {
    printf("%s %d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__);
    process_exit(1);
  }

  if (!(infile = fopen(filename,"r")))
  {
    printf("%s %d >> Can't open file %s\n", __FILE__, __LINE__, filename);
    process_exit(-1);
  }

  fread( (double*)Xpix1, 1, size, infile);

  fclose(infile);

  for( n=0; n<ni*nj*nk; n++)
  {
    a[stride*n] = Xpix1[n];
  }

  printf("%s %d >> write_dat() -- Read file \"%s\".\n", 
      __FILE__, __LINE__, filename);

} /* void write_dat( lattice_ptr lattice, double *a,  ... */


