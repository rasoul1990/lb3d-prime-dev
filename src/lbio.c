//##############################################################################
//
// lbio.c
//
//  - Lattice Boltzmann I/O routines.
//
//  - Mainly, dump the data to files that can be read by Matlab.
//
//  - Also, output routines for facilitating debugging.
//
//  - Should have a routine that dumps a matlab script?
//

// Some compilers, e.g., VC++, don't have the usual round() function
// in their math library.  Alternatively, ROUND can be defined as
// ceil or floor or some other rounding function.  It is used in
// the below routines for converting the real number valued of
// quantities at a lattice node into integer RGB values for writing
// to BMP files.
#define ROUND floor

//void output_frame( lattice_ptr lattice)
//##############################################################################
//
// O U T P U T   F R A M E
//
void output_frame( lattice_ptr lattice)
{
  double s, u[3];
  double nu;
  double L;
  int    n;

#if VERBOSITY_LEVEL > 0
      printf("\n");
      printf( "========================================"
              "========================================\n");
      printf("Begin file I/O at time = %d, frame = %d.\n",
          lattice->time, lattice->time/lattice->param.FrameRate);
      printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */
      dump_frame_summary( lattice);

#if WRITE_MACRO_VAR_DAT_FILES || WRITE_MACRO_VAR_RAW_FILES  || WRITE_PLOT_FILE
      dump_macro_vars( lattice, lattice->time);
#endif /* WRITE_MACRO_VAR_DAT_FILES || WRITE_MACRO_VAR_RAW_FILES */
#if 0
#if WRITE_PDF_DAT_FILES
      dump_pdf( lattice, lattice->time);
#endif /* WRITE_PDF_DAT_FILES */

      slice( lattice);
#endif
#if VERBOSITY_LEVEL > 0
      printf("\n");
      printf("File I/O done.\n");
      printf("--\n");
#endif /* VERBOSITY_LEVEL > 0 */

#if 0
  nu = (1./3.)*(lattice->param.tau[0] - .5);
  L = lattice->param.length_scale;

  compute_ave_u( lattice, u, 0);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("subs 0: Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("subs 0: Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("subs 0: Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );

#if NUM_FLUID_COMPONENTS == 2
  compute_ave_u( lattice, u, 1);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("subs 1: Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("subs 1: Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("subs 1: Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );
#endif /* NUM_FLUID_COMPONENTS == 2 */

#if STORE_UEQ
  compute_ave_ueq( lattice, u);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("eq:     Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("eq:     Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("eq:     Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );
#endif /* STORE_UEQ */
#endif

} /* void output_frame( lattice_ptr lattice) */

// void dump_frame_info( struct lattice_struct *lattice)
//##############################################################################
//
// D U M P   F R A M E   I N F O
//
void dump_frame_summary( struct lattice_struct *lattice)
{
  char   filename[1024];
  FILE   *o;
  double min_u[3], max_u[3],  ave_u[3];
  double min_rho, max_rho,   ave_rho;
  double rho_ratio, u_x_ratio, u_y_ratio, u_z_ratio;
  int    subs;

 for( subs = 0; subs < NUM_FLUID_COMPONENTS; subs++)
 {
  if( 1 || is_on_root_proc( lattice))
  {
    if( get_num_procs( lattice) > 1)
    {
      sprintf( filename, "./out/frames%dx%dx%d_subs%02d_proc%04d.dat",
        get_g_LX( lattice),
        get_g_LY( lattice),
        get_g_LZ( lattice), subs, get_proc_id( lattice));
    }
    else
    {
      sprintf( filename, "./out/frames%dx%dx%d_subs%02d.dat",
        get_g_LX( lattice),
        get_g_LY( lattice),
        get_g_LZ( lattice), subs );
    }

  // On the first timestep, make sure we start with a new file.
  if( lattice->time==0)
  {
      if( !( o = fopen(filename,"w+")))
      {
        printf("ERROR: fopen(\"%s\",\"w+\") = NULL.  Bye, bye!\n", filename);
        process_exit(1);
      }
      else
      {
        // Put a header on the file.
        fprintf( o, "\n");
        fprintf( o,
        "        time "
        "    min_u[x] "
        "    min_u[y] "
        "    min_u[z] "
        "    max_u[x] "
        "    max_u[y] "
        "    max_u[z] "
        "    ave_u[x] "
        "    ave_u[y] "
        "    ave_u[z] "
        "  max/ave[x] "
        "  max/ave[y] "
        "  max/ave[z] "
        "     min_rho "
        "     max_rho "
        "     ave_rho "
        "     max/ave "
        "\n");fprintf( o,
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "\n");
        fclose(o);
      }
  }
  if( !( o = fopen(filename,"a+")))
  {
    printf("ERROR: fopen(\"%s\",\"a+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }
  }

  compute_min_u( lattice, min_u, subs);
  compute_max_u( lattice, max_u, subs);
  compute_ave_u( lattice, ave_u, subs);
  compute_min_rho( lattice, &min_rho, subs);
  compute_max_rho( lattice, &max_rho, subs);
  compute_ave_rho( lattice, &ave_rho, subs);

  if( 1|| is_on_root_proc( lattice))
  {
    rho_ratio = ( ave_rho  != 0.) ? ( max_rho /ave_rho ):( 1.);
    u_x_ratio = ( ave_u[0] != 0.) ? ( max_u[0]/ave_u[0]):( 1.);
    u_y_ratio = ( ave_u[1] != 0.) ? ( max_u[1]/ave_u[1]):( 1.);
    u_z_ratio = ( ave_u[2] != 0.) ? ( max_u[2]/ave_u[2]):( 1.);

    fprintf( o,
      "%12d "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f\n",
      lattice->time,
      min_u[0], min_u[1], min_u[2],
      max_u[0], max_u[1], max_u[2],
      ave_u[0], ave_u[1], ave_u[2],
      (fabs(u_x_ratio)<=9999.)?(u_x_ratio):(1./0.),
      (fabs(u_y_ratio)<=9999.)?(u_y_ratio):(1./0.),
      (fabs(u_z_ratio)<=9999.)?(u_z_ratio):(1./0.),
      min_rho,
      max_rho,
      ave_rho,
      (fabs(rho_ratio)<=9999.)?(rho_ratio):(1./0.) );

  fclose(o);

#if VERBOSITY_LEVEL > 0
    printf("%s %d %04d >> dump_frame_info() -- "
      "Wrote/appended to file \"%s\"\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
#endif /* VERBOSITY_LEVEL > 0 */

#if VERBOSITY_LEVEL > 0
    printf("%s %d %04d >> dump_frame_info() -- "
      "frame = %d/%d = %d\n",
      __FILE__, __LINE__, get_proc_id(lattice),
      lattice->time,
      lattice->param.FrameRate,
      (int)((double)lattice->time/(double)lattice->param.FrameRate));
#endif /* VERBOSITY_LEVEL > 0 */
  }
 }

} /* void dump_frame_info( struct lattice_struct *lattice) */

// void dump_macro_vars( struct lattice_struct *lattice)
//##############################################################################
//
// D U M P   M A C R O S C O P I C
//
//  - Output the macro_vars variables to files.
//
void dump_macro_vars( struct lattice_struct *lattice, int time)
{
  char   filename[1024],filename2[1024], fn3[1024];
  FILE   *o, *o_u, *o_rho, *o_ux, *o_uy, *o_ueq, *o_ueq_x, *o_ueq_y, *sat;
  int    *node_ptr;
  int    i, j, n, count, count1, count2, basek;
  double *macro_vars_ptr;
  double *ueq;
  int    frame;
  double min_u[3], max_u[3],  ave_u[3];
  double min_rho, max_rho,   ave_rho;
  double rho_ratio, u_x_ratio, u_y_ratio;
  int    subs;
  int    ni, nj, nk;
  int    j_slice, k_slice;

  frame = (int)((double)lattice->time/(double)lattice->param.FrameRate);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    compute_max_rho( lattice, &max_rho, subs);
    compute_min_rho( lattice, &min_rho, subs);
    compute_ave_rho( lattice, &ave_rho, subs);

   if( is_on_root_proc( lattice))
    {
    printf("%s %d %04d >> "
      "min_rho = %20.17f, max_rho = %20.17f, ave_rho = %20.17f\n",
      __FILE__,__LINE__, get_proc_id(lattice),
      min_rho, max_rho, ave_rho);
   } /* if( is_on_root_proc( lattice)) */

#if WRITE_MACRO_VAR_RAW_FILES
    sprintf( filename, "./out/rho_%dx%dx%d_frame%04d_subs%02d_proc%04d.raw",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_raw(
      lattice,
      &(lattice->macro_vars[subs][0].rho),
      /*stride*/ 4,
      max_rho,
      /*min_rho*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_RAW_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/rho_%dx%dx%d_frame%04d_subs%02d_proc%04d.dat",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_dat(
      lattice,
      &(lattice->macro_vars[subs][0].rho),
      /*stride*/ 4,
      max_rho,
      /*min_rho*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */

    compute_max_u( lattice, max_u, subs);
    compute_min_u( lattice, min_u, subs);
    compute_ave_u( lattice, ave_u, subs);


   if( is_on_root_proc( lattice))
   {
    printf("%s %d %04d >> "
      "min_u = [ %f %f %f], max_u = [ %f %f %f], ave_u = [ %f %f %f]\n",
      __FILE__,__LINE__, get_proc_id(lattice),
    min_u[0], min_u[1], min_u[2],
    max_u[0], max_u[1], max_u[2],
    ave_u[0], ave_u[1], ave_u[2] );
   } /* if( is_on_root_proc( lattice)) */

 /*
   j_slice = (int)floor((double)nj/2.);

    k_slice = 0;
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ", get_uz( lattice, subs, XYZ2N(i,j_slice,k_slice,ni,nj)));
    }
    printf("]\n");

    k_slice = 1;
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ", get_uz( lattice, subs, XYZ2N(i,j_slice,k_slice,ni,nj)));
    }
    printf("]\n");

    k_slice = (int)floor((double)nk/4.);
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ",
        lattice->macro_vars[subs][XYZ2N(i,j_slice,k_slice,ni,nj)].u[2]);
    }
    printf("]\n");

    k_slice = (int)floor((double)nk/2.);
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ",
        lattice->macro_vars[subs][XYZ2N(i,j_slice,k_slice,ni,nj)].u[2]);
    }
    printf("]\n");

    k_slice = (int)floor(3.*(double)nk/4.);
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ",
        lattice->macro_vars[subs][XYZ2N(i,j_slice,k_slice,ni,nj)].u[2]);
    }
    printf("]\n");


    j_slice = (int)floor((double)nj/2.);
    k_slice = get_LZ( lattice) - 1;
    printf("%s %d %04d >> k_slice(%d) = [",
      __FILE__,__LINE__, get_proc_id(lattice), k_slice+get_g_SZ(lattice));
    for( i=0; i<ni; i++)
    {
      printf(" %f ",
        lattice->macro_vars[subs][XYZ2N(i,j_slice,k_slice,ni,nj)].u[2]);
    }
    printf("]\n");
*/

#if WRITE_MACRO_VAR_RAW_FILES
    sprintf( filename, "./out/u_%dx%dx%d_frame%04d_subs%02d_proc%04d.raw",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_raw_u(
      lattice,
      &(lattice->macro_vars[subs][0].u[0]),
      &(lattice->macro_vars[subs][0].u[1]),
      &(lattice->macro_vars[subs][0].u[2]),
      /*stride*/ 4,
      max_u[0],
      /*min_u[0]*/0.,
      max_u[1],
      /*min_u[1]*/0.,
      max_u[2],
      /*min_u[2]*/0.,
      filename);
#endif
#if WRITE_MACRO_VAR_RAW_FILES
    sprintf( filename, "./out/u_x_%dx%dx%d_frame%04d_subs%02d_proc%04d.raw",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_raw(
      lattice,
      &(lattice->macro_vars[subs][0].u[0]),
      /*stride*/ 4,
      max_u[0],
      /*min_u[0]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_RAW_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_x_%dx%dx%d_frame%04d_subs%02d_proc%04d.dat",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_dat(
      lattice,
      &(lattice->macro_vars[subs][0].u[0]),
      /*stride*/ 4,
      max_u[0],
      /*min_u[0]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */
#if WRITE_MACRO_VAR_RAW_FILES
    sprintf( filename, "./out/u_y_%dx%dx%d_frame%04d_subs%02d_proc%04d.raw",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_raw(
      lattice,
      &(lattice->macro_vars[subs][0].u[1]),
      /*stride*/ 4,
      max_u[1],
      /*min_u[1]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_RAW_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_y_%dx%dx%d_frame%04d_subs%02d_proc%04d.dat",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_dat(
      lattice,
      &(lattice->macro_vars[subs][0].u[1]),
      /*stride*/ 4,
      max_u[1],
      /*min_u[1]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */
#if WRITE_MACRO_VAR_RAW_FILES
    sprintf( filename, "./out/u_z_%dx%dx%d_frame%04d_subs%02d_proc%04d.raw",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_raw(
      lattice,
      &(lattice->macro_vars[subs][0].u[2]),
      /*stride*/ 4,
      max_u[2],
      /*min_u[2]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_RAW_FILES) */
#if WRITE_MACRO_VAR_DAT_FILES
    sprintf( filename, "./out/u_z_%dx%dx%d_frame%04d_subs%02d_proc%04d.dat",
        ni, nj, nk, frame, subs, get_proc_id( lattice));
    write_dat(
      lattice,
      &(lattice->macro_vars[subs][0].u[2]),
      /*stride*/ 4,
      max_u[2],
      /*min_u[2]*/0.,
      filename);
#endif /* (WRITE_MACRO_VAR_DAT_FILES) */

//------------------------------------------------------------------------------------------
//
//  Write Vmag files for post processing large domains that have memory issues
//
//  edit: peter
//
#if WRITE_VMAG

int eye;
float vx, vy, vz, vxsq, vysq, vzsq;
float mag=0;
double vmagmax=0, vmagmin=0;
float vmag[ni*nj*nk];
double vmagscale[ni*nj*nk];

  for(eye=0; eye<(ni*nj*nk);eye++)
  {
    vx=lattice->macro_vars[subs][eye].u[0];
    vxsq = vx*vx;
    vy=lattice->macro_vars[subs][eye].u[1];
    vysq = vy*vy;
    vz=lattice->macro_vars[subs][eye].u[2];
    vzsq = vz*vz;
    mag = sqrt(vxsq+vysq+vzsq);
    vmag[eye] = mag;
    if(vmag[eye]>vmagmax)
    {
      vmagmax = vmag[eye];
    }
    if(vmag[eye]<vmagmin)
    {
      vmagmin = vmag[eye];
    }
  }

char *ffilename[1024];
FILE *out;
  sprintf( ffilename,"./out/vmag_proc%04d.txt", get_proc_id( lattice));
  out = fopen( ffilename, "w");
    for(eye=0; eye<(ni*nj*nk);eye++)
    {
      fprintf(out,"%f\n",vmag[eye]);
    //  vmagscale[eye] = (255.*vmag[eye] / vmagmax);
    //  printf("vmag [%d] = %1.10f,   vmagscale[%d] = %3.5f \n", eye, vmag[eye], eye, vmagscale[eye]);
    }
  fclose(out);
//printf("vmagmax = %f, vmagmin = %f \n", vmagmax, vmagmin);
//printf("velocity [%d] = %1.10f, Squared = %1.10f \n", eye, lattice->macro_vars[subs][eye].u[0], (lattice->macro_vars[subs][eye].u[0]*lattice->macro_vars[subs][eye].u[0]));
#endif /*WRITE_VMAG*/
//-------------------------------------------------------------------------------------------

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if WRITE_PLOT_FILE
    if(NUM_FLUID_COMPONENTS == 2)
    {
    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    sprintf( filename2, "./out/R2_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 2 component\n");
    write_plt(
      lattice,
      &(lattice->macro_vars[0][0].rho),
      &(lattice->macro_vars[1][0].rho),
      filename, filename2);

    //get saturation
   basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif
        sprintf( fn3, "./out/saturation_proc%04d.dat",
        get_proc_id( lattice));
  if(frame ==0)
  {
  if(  !( sat = fopen(fn3, "w+")))
   {
    printf(" %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      fn3);
    process_exit(1);
   }
  }
  if(frame !=0)
  {
  if(  !( sat = fopen(fn3, "a+")))
   {
    printf(" %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      fn3);
    process_exit(1);
   }
  }
  count =0;
  count1 =0;
  count2 =0;

 for( n=0; n<ni*nj*nk; n++)
  {

     if( (N2Z(n,ni,nj,nk)+basek >lattice->param.z1 &&
          N2Z(n,ni,nj,nk)+basek <lattice->param.z2    ) )
     {
       count++;
      if ( lattice->macro_vars[0][n].rho> (lattice->param.rho_A[0]+lattice->param.rho_A[1])/2.) count1 ++;
      if ( lattice->solids[0][n].is_solid) count2++;
     }
  }
  fprintf( sat, " %7d, %7d, %7d, %7d\n",   count, count2, count-count2, count1  );


    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_uvw.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 2 UVZ component\n");
    write_plt_uvw(
      lattice,
      &(lattice->macro_vars[0][0].rho),
      &(lattice->macro_vars[1][0].rho),
      filename);
 fclose(sat);
    }

    else
    {
    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 1 component\n");
    write_plt_single(
      lattice,
      &(lattice->macro_vars[0][0].rho),
      filename);
    }
#endif

} /* void dump_macro_vars( struct lattice_struct *lattice, int time) */

void read_solids( lattice_ptr lattice, char *filename)
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *fd;
#if 1
  int i, j, k, p;
#endif

  for( n=0; n<lattice->NumNodes; n++)
  {
    lattice->solids[0][n].is_solid = 0;
  }

  fd = fopen( filename, "r+");
  if( !fd)
  {
    printf("%s %d %04d >> ERROR: Can't open file \"%s\". (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice), filename);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't malloc image buffer. (Exiting!)\n", __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  printf("%s %d %04d >> Reading %d bytes from file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), size, filename);

  //size_read = read( fd, raw, size);
  size_read = fread( raw, 1, size, fd);

  if( size_read != size)
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't read image data: read = %d. (Exiting!)\n",
          __FILE__, __LINE__, get_proc_id(lattice), size_read);
    process_exit(1);
  }

  fclose( fd);

  g_n = get_g_StartNode( lattice); // Global node index.
  for( n=0; n<lattice->NumNodes; n++)
  {
    for( subs=0; subs<(NUM_FLUID_COMPONENTS); subs++)
    {
      lattice->solids[subs][n].is_solid = raw[g_n];
    }
    //printf("%s %d %04d >> solids[%d] = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice), n, (int)lattice->solids[0][n].is_solid);
    g_n++;
  }

  free(raw);

#if 1
  if( /* Domain not too big. (Tweak to suit.) */
      ( get_LX(lattice) <= 12
      &&
        get_LY(lattice) <= 12
      &&
        get_LZ(lattice) <= 12
      )
    )
  {
#if PARALLEL
  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);
    if( p == get_proc_id( lattice))
    {
      printf("%s %d %04d >> Solids:\n", __FILE__, __LINE__, p);
#endif
      for( j=0; j<get_LY( lattice); j++)
      {
        for( k=0; k<get_LZ( lattice); k++)
        {
          for( i=0; i<get_LX( lattice); i++)
          {
            if( is_solid( lattice,
                          XYZ2N( i, j, k,
                                 get_LX( lattice),
                                 get_LY( lattice)) ) )
            {
              printf("#");
            }
            else
            {
              printf(" ");
            }
          }
          printf(" ");
        }
        printf("\n");
      }
#if PARALLEL
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
#endif
  }
#endif

} /* read_solids( lattice_ptr lattice, char *filename) */


/*void read_solids_from_plt( lattice_ptr lattice, char *filename)
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *in;
#if 1
  int i, j, k, p;
#endif


  for( n=0; n<lattice->NumNodes; n++)
  {
    lattice->solids[0][n].is_solid = 0;
  }


   if( !( in = fopen( filename, "r")))
  {
    printf("%s %d %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      infile);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't malloc image buffer. (Exiting!)\n", __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }
 for (i =1 ; i<=7; i++)
 {
   skip_label( in);
 }

  for(k =1; k<=lattice->param.LZ; k++)
  {
  for(j =1; j<=lattice->param.LY; j++)
  {
  for(i =1; i<=lattice->param.LX; i++)
  {
   fscanf( in, "%d", &( raw[n])         );
  }
  }
  }

  fclose( in);

  for( n=0; n<lattice->NumNodes; n++)
  {
    for( subs=0; subs<(NUM_FLUID_COMPONENTS); subs++)
    {
      lattice->solids[subs][n].is_solid = raw[g_n]*255;
    }
    //printf("%s %d %04d >> solids[%d] = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice), n, (int)lattice->solids[0][n].is_solid);
    g_n++;
  }

  free(raw);

#if 1
  if( // Domain not too big. (Tweak to suit.)
      ( get_LX(lattice) <= 12
      &&
        get_LY(lattice) <= 12
      &&
        get_LZ(lattice) <= 12
      )
    )
  {
#if PARALLEL
  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);
    if( p == get_proc_id( lattice))
    {
      printf("%s %d %04d >> Solids:\n", __FILE__, __LINE__, p);
#endif
      for( j=0; j<get_LY( lattice); j++)
      {
        for( k=0; k<get_LZ( lattice); k++)
        {
          for( i=0; i<get_LX( lattice); i++)
          {
            if( is_solid( lattice,
                          XYZ2N( i, j, k,
                                 get_LX( lattice),
                                 get_LY( lattice)) ) )
            {
              printf("#");
            }
            else
            {
              printf(" ");
            }
          }
          printf(" ");
        }
        printf("\n");
      }
#if PARALLEL
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
#endif
  }
#endif
}
*/

void write_raw(
       lattice_ptr lattice,
       double *a,
       int     stride,
       double  a_max,
       double  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  FILE *infile;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
    if( !lattice->solids[0][n].is_solid)
    {
//      Xpix1[n] = (unsigned char)ROUND(255.*(a[stride*n]-a_min)/(a_max-a_min));
      Xpix1[n] = (unsigned char)ROUND(1.+254.*(a[stride*n]-a_min)/(a_max-a_min));
    }
    else
    {
      Xpix1[n] = (unsigned char)0;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} /* void write_raw( lattice_ptr lattice, double *a,  ... */


void write_raw_u(
       lattice_ptr lattice,
       double *ux,
       double *uy,
       double *uz,
       int     stride,
       double  ux_max,
       double  ux_min,
       double  uy_max,
       double  uy_min,
       double  uz_max,
       double  uz_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  FILE *infile;
  double a_max = sqrt( ux_max*ux_max + uy_max*uy_max + uz_max*uz_max);
  double a_min = sqrt( ux_min*ux_min + uy_min*uy_min + uz_min*uz_min);
  double a;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
    if( !lattice->solids[0][n].is_solid)
    {
//      Xpix1[n] = (unsigned char)ROUND(255.*(a[stride*n]-a_min)/(a_max-a_min));
      a = sqrt( ux[stride*n]*ux[stride*n]
              + uy[stride*n]*uy[stride*n]
              + uz[stride*n]*uz[stride*n]);
      Xpix1[n] = (unsigned char)ROUND(1.+254.*(a-a_min)/(a_max-a_min));
    }
    else
    {
      Xpix1[n] = (unsigned char)0;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} /* void write_raw( lattice_ptr lattice, double *a,  ... */


/*void write_raw1(
       lattice_ptr lattice,
       double *a,
       int     stride,
       double  a_max,
       double  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  double half;
  FILE *infile;
//  half = (a_max+ a_min)/2.;
    half = lattice->param.rho_A[0]-0.5;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
//       R      G       B
//14  88  0  0
//15  114  0  1
//16  140  0  1
//17  166  0  1        decane
//18  192  0  1
//19  218  0  1
//20  244  0  1
//
//33  255  0  2
//34  255  0  78
//35  255  0  171
//36  127  0  254
//37  0  0  254
//38  0  0  254      water
//
//64  0  0  254
//65  0  0  178
//66  88  36  101
//67  213  198  24
//68  254  254  0       solid

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
    if( !lattice->solids[0][n].is_solid && a[stride*n]<half)
    {
      Xpix1[n] = (unsigned char)ROUND(64.-29.*(a[stride*n]-a_min)/(half-a_min));
//According to Utah color Table, please refer to Utah Color Table!  35--64
    }
    else if(!lattice->solids[0][n].is_solid && a[stride*n]>half)
    {
      Xpix1[n] = (unsigned char)ROUND(35.-15.*(a[stride*n]-half)/(a_max-half));
//According to Utah color Table, please refer to Utah Color Table!  20---35
    }

    else
    {
      Xpix1[n] = (unsigned char)255.;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} // void write_raw( lattice_ptr lattice, double *a,  ... */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void write_dat(
       lattice_ptr lattice,
       double *a,
       int     stride,
       double  a_max,
       double  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  double*Xpix1;
  FILE *infile ;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(double);

  if(!(Xpix1 = ( double*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  for( n=0; n<ni*nj*nk; n++)
  {
    if( !lattice->solids[0][n].is_solid)
    {
      Xpix1[n] = a[stride*n];
    }
    else
    {
      Xpix1[n] = 0.;
    }
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }

  fwrite( (double*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_dat() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_dat( lattice_ptr lattice, double *a,  ... */

#if (NUM_FLUID_COMPONENTS==2)
void write_plt(
       lattice_ptr lattice,
       double *a, double *b,
       char   *filename, char *filename2)
{
  int size;
  int n, ni, nj, nk, count, basek;
  FILE *infile, *infile2, *FL1, *FL2, *FL3, *FL4;
  double v_in, v_out, v_inupr, v_outupr;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif

  size = ni*nj*nk*sizeof(double);


  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  if (!(infile2 = fopen(filename2,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }


  //  fprintf( infile, "variables = x, y, z, rho1, rho2, solid\n"  );
  fprintf( infile, "variables =rho1\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );
  fprintf( infile2, "variables =rho2\n"  );
  fprintf( infile2, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );

  count = 0;

  if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  {

  if (!(FL1 = fopen("./out/in_vels.dat","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  if (!(FL3 = fopen("./out/zupper.plt","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  //ptinf out the 2D plt file of upper inlet surface
  fprintf(FL3,"variables = i, j, rho1, u, v, w\n");
  fprintf(FL3,"zone i=%d, j=%d,  f=point\n", ni,  nj);
   v_in = 0.;
   v_inupr = 0.;
  }
  // endif get_proc_id(lattice) == get_num_procs(lattice)-1)
  if(get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }

  if (!(FL4 = fopen("./out/zbottom.plt","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL4,"variables = i, j, rho1, u, v, w\n");
  fprintf(FL4,"zone i=%d, j=%d,  f=point\n", ni,  nj);

   v_out = 0.;
   v_outupr = 0.;
  }
  //endif (get_proc_id(lattice) == 0

  for( n=0; n<ni*nj*nk; n++)
  {
   if( lattice->solids[0][n].is_solid)
     {
      a[4*n] = -1.;
      lattice->macro_vars[1][n].rho = -1.;
     }

     if(  get_proc_id(lattice) == get_num_procs(lattice)-1 &&
N2Z(n,ni,nj,nk)+basek == (get_num_procs(lattice) *lattice->param.LZ-1) )
     {
       v_in= v_in + a[4*n+3];
       v_inupr= v_inupr +lattice->ueq[n].u[2];

    fprintf(FL3,"%4d %4d %7.4f  %2.9f  %2.9f %2.9f \n",N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), a[4*n], a[4*n+1], a[4*n+2], a[4*n+3] ); //MS changed format from 8.5 to 2.9
     }
     if(get_proc_id(lattice) == 0 && N2Z(n,ni,nj,nk)+basek == 1)
     {
       v_out= v_out + a[4*n+3];
             v_outupr= v_outupr +lattice->ueq[n].u[2];
    fprintf(FL4,"%4d %4d %7.4f  %2.9f  %2.9f %2.9f \n",N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), a[4*n], a[4*n+1], a[4*n+2], a[4*n+3] ); //MS changed format from 8.5 to 2.9

     }
     if( (N2Z(n,ni,nj,nk)+basek >lattice->param.z1 &&
    N2Z(n,ni,nj,nk)+basek <lattice->param.z2
           ) && (a[4*n] > (lattice->param.rho_A[0]+lattice->param.rho_A[1])/2.) )
     {count++;}
//  fprintf( infile, "%4d  %4d %4d  %7.4f  %7.4f  %2d\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk),
//  a[/*stride*/ 4*n], lattice->macro_vars[1][n].rho , lattice->solids[0][n].is_solid/255  );
  fprintf( infile, " %7.4f \n",   a[/*stride*/ 4*n]  );
  fprintf( infile2, " %7.4f \n",   b[/*stride*/ 4*n]  );
  }
  if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  fprintf(FL1,"position=%4d",(get_num_procs(lattice) *lattice->param.LZ-2 ));
  fprintf(FL1, "integration in =%12.7f\n", v_in );
  fprintf(FL1, "integration in upr=%12.7f\n", v_inupr );
  fclose(FL1);
  fclose(FL3);
  }
  if(get_proc_id(lattice) == 0 )
  {
  fprintf(FL2, "integration out=%12.7f\n", v_out );
  fprintf(FL2, "integration out upr =%12.7f\n", v_outupr );
  fclose(FL2);
  fclose(FL4);
  }

  fclose(infile);
  fclose(infile2);

  printf("%s %d %04d >> write_tecplot_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
  printf("%s %d >> Non-Wetting fluid nodes: %d \n",
      __FILE__, __LINE__, count);

} /* void write_plt( lattice_ptr lattice, double *a,  ... */
#endif

void write_plt_uvw(
       lattice_ptr lattice,
       double *a, double *b,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  FILE *infile;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(double);

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf( infile, "variables = x, y, z, u, v, w\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );


 for( n=0; n<ni*nj*nk; n++)
  {
   if( lattice->solids[0][n].is_solid)
     {
      a[4*n] = -1.;
      lattice->macro_vars[1][n].rho = -1.;
     }

//   if(NUM_FLUID_COMPONENT ==2)
   {
  fprintf( infile, "%4d  %4d %4d  %7.4f  %7.4f  %7.4f\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk),
    b[/*stride*/ 4*n+1], b[/*stride*/ 4*n+2], b[/*stride*/ 4*n+3]  );
   }

  }

  fclose(infile);

  printf("%s %d %04d >> write_tecplot_uvw_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_plt_uvw( lattice_ptr lattice, double *a,  ... */


void write_plt_single(
       lattice_ptr lattice,
       double *a,
       char   *filename)
{
  int size;
  int n, ni, nj, nk, basek;
  double v_in, v_out;
  FILE *infile, *FL1, *FL2;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
  basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif

  size = ni*nj*nk*sizeof(double);

  if(lattice->time==0 && get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  if (!(FL1 = fopen("./out/in_vels.dat","w+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL1, "time  integration in \n" );

  }//endif   if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  if(lattice->time!=0 && get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  if (!(FL1 = fopen("./out/in_vels.dat","a+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  }//endif   if(get_proc_id(lattice) == get_num_procs(lattice)-1)

  if (lattice->time==0 && get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","w+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL2, "time  integration in \n" );

  }
//if   if(get_proc_id(lattice) == 0
  if (lattice->time!=0 && get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","a+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  }
//if   if(get_proc_id(lattice) == 0

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf( infile, "variables = x, y, z, rho1, u, v, w, solid\n"  );
//  fprintf( infile, "variables = rho1, w, solid\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );

  v_in = 0.;
  v_out = 0.;
 for( n=0; n<ni*nj*nk; n++)
  {
   if( lattice->solids[0][n].is_solid)
     {
      a[4*n] = -1.;
     }
   else
     {
       if(  (get_proc_id(lattice) == get_num_procs(lattice)-1) &&
           N2Z(n,ni,nj,nk)+basek == get_num_procs(lattice) *lattice->param.LZ-1) v_in= v_in + a[4*n+3];
       if(  (get_proc_id(lattice) == 0) && N2Z(n,ni,nj,nk)+basek ==1) v_out= v_out + a[4*n+3];
     }
  fprintf( infile, "%4d  %4d %4d  %7.4f  %2.9f %2.9f %2.9f %2d\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk)+basek,
  a[/*stride*/ 4*n], a[/*stride*/ 4*n+1], a[/*stride*/ 4*n+2], a[/*stride*/ 4*n+3], lattice->solids[0][n].is_solid/255  ); //MS changed format from 8.5 to 2.9

//  fprintf( infile, " %7.4f  %7.4f  %2d\n",
//  a[/*stride*/ 4*n],  a[/*stride*/ 4*n+3], lattice->solids[0][n].is_solid/255  );
  }

if(get_proc_id(lattice) == get_num_procs(lattice)-1) //Works only if Parallel (MS)

  {
  fprintf(FL1, "%8d  %12.7f\n", lattice->time, v_in );
  fclose(FL1);
  }

  if(get_proc_id(lattice) == 0)
  {
  fprintf(FL2, "%8d  %12.7f\n", lattice->time, v_out );
  fclose(FL2);
  }

  fclose(infile);

  printf("%s %d %04d >> write_tecplot_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_plt( lattice_ptr lattice, double *a,  ... */
