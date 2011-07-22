//##############################################################################
//
// compute.c
//
//  - Routines for computing on the lattice:
//
//    - compute_rho_and_u
//    - compute_feq
//    - compute_big_u
//    - compute_gforce
//    - compute_fluid_fluid_force
//    - etc...
//

// void compute_macro_vars( struct lattice_struct *lattice)
//##############################################################################
//
// C O M P U T E   M A C R O   V A R S
//
//  - Compute macroscopic variables.

void rho_send_recv_begin( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n, a;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

//rev Huang
//To send "rho" and "solid" to neighbour blocks
  n =0;
  k = get_LZ(lattice)-1;
// for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 for( a=0; a<NUM_FLUID_COMPONENTS; a++)
   //changed by Shadab on Danny's suggestion!!
  {
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
#if SAVE_MEMO

#else

      lattice->process.z_pos_rho_to_send[n] =
        lattice->macro_vars[a][ XYZ2N( i , j , k, ni, nj)].rho;
      lattice->process.z_neg_rho_to_send[n] =
        lattice->macro_vars[a][ XYZ2N( i , j , 0, ni, nj)].rho;
      n++;
#endif
    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */
  } //for  subs=0; subs<NUM_FLUID_COMPONENTS; subs++

     // S E N D   I N   P O S I T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Isend( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Isend(
       /*void *buf*/          lattice->process.z_pos_rho_to_send,
       /*int count*/          2*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_DOUBLE,
       /*int dest*/         ( get_proc_id(lattice)
                            + get_num_procs(lattice)+1)
                            % get_num_procs(lattice),
       /*int tag*/            2,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.send_req_2)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Isend( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // R E C V   F R O M   N E G A T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Irecv( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Irecv(
       /*void *buf*/          lattice->process.z_pos_rho_to_recv,
       /*int count*/          2*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_DOUBLE,
       /*int src*/          ( get_proc_id(lattice)
                            + get_num_procs(lattice)-1)
                            % get_num_procs(lattice),
       /*int tag*/            2,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.recv_req_2)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Irecv( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // S E N D   I N   N E G A T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Isend( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Isend(
       /*void *buf*/          lattice->process.z_neg_rho_to_send,
       /*int count*/          2*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_DOUBLE,
       /*int dest*/         ( get_proc_id(lattice)
                            + get_num_procs(lattice)-1)
                            % get_num_procs(lattice),
       /*int tag*/            3,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.send_req_3)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Isend( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // R E C V   F R O M   P O S I T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Irecv( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Irecv(
       /*void *buf*/          lattice->process.z_neg_rho_to_recv,
       /*int count*/          2*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_DOUBLE,
       /*int src*/          ( get_proc_id(lattice)
                            + get_num_procs(lattice)+1)
                            % get_num_procs(lattice),
       /*int tag*/            3,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.recv_req_3)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Irecv( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }



#endif
} /* void rho_send_recv_begin( lattice_ptr lattice, const int subs) */


void rho_send_recv_end( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int ip, in;
  int jp, jn;
  int mpierr;

  mpierr = MPI_Wait(
  /* mpi_request *req */&(lattice->process.send_req_2),
  /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* mpi_request *req */&(lattice->process.recv_req_2),
  /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* mpi_request *req */&(lattice->process.send_req_3),
  /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* mpi_request *req */&(lattice->process.recv_req_3),
  /* mpi_status *stat */&(lattice->process.mpi_status));

#endif
} /* void rho_send_recv_end( lattice_ptr lattice, const int subs) */
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void solid_send_recv_begin( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

//rev Huang
//To send "rho" and "solid" to neighbour blocks
  n =0;
  k = get_LZ(lattice)-1;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
#if SAVE_MEMO

#else
      lattice->process.z_pos_solid_to_send[n] =
        lattice->solids[subs][ XYZ2N( i , j , k, ni, nj)].is_solid;
      lattice->process.z_neg_solid_to_send[n] =
        lattice->solids[subs][ XYZ2N( i , j , 0, ni, nj)].is_solid;
      n++;
#endif
    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */


     // SOLID S E N D   I N   P O S I T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Isend( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Isend(
       /*void *buf*/          lattice->process.z_pos_solid_to_send,
       /*int count*/          1*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_INT,
       /*int dest*/         ( get_proc_id(lattice)
                            + get_num_procs(lattice)+1)
                            % get_num_procs(lattice),
       /*int tag*/            4,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.send_req_4)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Isend( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // SOLID R E C V   F R O M   N E G A T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Irecv( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Irecv(
       /*void *buf*/          lattice->process.z_pos_solid_to_recv,
       /*int count*/          1*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_INT,
       /*int src*/          ( get_proc_id(lattice)
                            + get_num_procs(lattice)-1)
                            % get_num_procs(lattice),
       /*int tag*/            4,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.recv_req_4)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Irecv( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // SOLID S E N D   I N   N E G A T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Isend( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Isend(
       /*void *buf*/          lattice->process.z_neg_solid_to_send,
       /*int count*/          1*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_INT,
       /*int dest*/         ( get_proc_id(lattice)
                            + get_num_procs(lattice)-1)
                            % get_num_procs(lattice),
       /*int tag*/            5,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.send_req_5)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Isend( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }
     // SOLID R E C V   F R O M   P O S I T I V E   D I R E C T I O N
     //#########################################################################
#if VERBOSITY_LEVEL > 1
     printf( "%s %d %04d >> "
       "MPI_Irecv( %04d)"
       "\n",
       __FILE__,__LINE__,get_proc_id(lattice),
       (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
     mpierr =
       MPI_Irecv(
       /*void *buf*/          lattice->process.z_neg_solid_to_recv,
       /*int count*/          1*(get_LX(lattice)*get_LY(lattice)),
       /*MPI_Datatype dtype*/ MPI_INT,
       /*int src*/          ( get_proc_id(lattice)
                            + get_num_procs(lattice)+1)
                            % get_num_procs(lattice),
       /*int tag*/            5,
       /*MPI_Comm comm*/      MPI_COMM_WORLD,
       /*MPI_Request *req*/   &(lattice->process.recv_req_5)
       );
     if( mpierr != MPI_SUCCESS)
     {
       printf( "%s %d %04d >> "
         "ERROR: %d <-- MPI_Irecv( %04d)"
         "\n",
         __FILE__,__LINE__,get_proc_id(lattice),
         mpierr,
         (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
       process_finalize();
       process_exit(1);
     }



#endif
} /* void rho_send_recv_begin( lattice_ptr lattice, const int subs) */


void solid_send_recv_end( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int ip, in;
  int jp, jn;
  int mpierr;

  mpierr = MPI_Wait(
  /* MPI_Request *req */&(lattice->process.send_req_4),
  /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* MPI_Request *req */&(lattice->process.recv_req_4),
  /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* MPI_Request *req */&(lattice->process.send_req_5),
  /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
  /* MPI_Request *req */&(lattice->process.recv_req_5),
  /* MPI_Status *stat */&(lattice->process.mpi_status));


#endif
} /* void rho_send_recv_end( lattice_ptr lattice, const int subs) */
//#####################################################################
//#####################################################################


#if INAMURO_SIGMA_COMPONENT
// TODO:
#else /* !( INAMURO_SIGMA_COMPONENT) */
void compute_macro_vars( struct lattice_struct *lattice)
{
  int a, n, k, i,j, ni, nj, nk;

  double *rho[ NUM_FLUID_COMPONENTS],
         *u_x[ NUM_FLUID_COMPONENTS],
         *u_y[ NUM_FLUID_COMPONENTS],
         *u_z[ NUM_FLUID_COMPONENTS];

  double ux_sum, uy_sum, uz_sum;

  double *ueq;

  double *ftemp;

  int    is_solid;

  int    subs;

  double tau0,
         tau1;

  int mpierr;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  for(subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {

    for( n=0; n<lattice->NumNodes; n++)
    {
      rho[subs] = &( lattice->macro_vars[subs][n].rho);
      u_x[subs] = &( lattice->macro_vars[subs][n].u[0]);
      u_y[subs] = &( lattice->macro_vars[subs][n].u[1]);
      u_z[subs] = &( lattice->macro_vars[subs][n].u[2]);

      ftemp     =   lattice->pdf[subs][n].ftemp;
      is_solid  = ( lattice->solids[subs][n].is_solid);

      if(!(is_solid))
      {
        *rho[subs] = 0.;
        *u_x[subs] = 0.;
        *u_y[subs] = 0.;
        *u_z[subs] = 0.;

        for(a=0; a<Q; a++)
        {
          (*rho[subs]) +=       (ftemp[a]);
          (*u_x[subs]) += vx[a]*(ftemp[a]);
          (*u_y[subs]) += vy[a]*(ftemp[a]);
          (*u_z[subs]) += vz[a]*(ftemp[a]);
        }

      } /*if(!(is_solid))*/
    } /*for( n=0; n<lattice->NumNodes; n++) */
  } /* for(subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

  if( NUM_FLUID_COMPONENTS==2)
  {
    tau0 = lattice->param.tau[0];
    tau1 = lattice->param.tau[1];

    for( n=0; n<lattice->NumNodes; n++)
    {
      for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
      {
        rho[subs]   = &( lattice->macro_vars[subs][n].rho);
        u_x[subs]   = &( lattice->macro_vars[subs][n].u[0]);
        u_y[subs]   = &( lattice->macro_vars[subs][n].u[1]);
        u_z[subs]   = &( lattice->macro_vars[subs][n].u[2]);
      }

#if STORE_UEQ
      ueq = lattice->ueq[n].u;
#endif /* STORE_UEQ */

      is_solid  = ( lattice->solids[0][n].is_solid);

      if( !( is_solid))
      {
        ux_sum =  *u_x[0]/tau0 + *u_x[1]/tau1;
        uy_sum =  *u_y[0]/tau0 + *u_y[1]/tau1;
        uz_sum =  *u_z[0]/tau0 + *u_z[1]/tau1;

#if STORE_UEQ
        if( *rho[0] + *rho[1] != 0.)
        {
          //        ueq[0] = ( ux_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          //        ueq[1] = ( uy_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          //        ueq[2] = ( uz_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( ux_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( uy_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( uz_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
        }
        else
        {
          *ueq++ = 0.;
          *ueq++ = 0.;
          *ueq++ = 0.;
        }
#endif /* STORE_UEQ */


        if( ux_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_x[0] = *u_x[0] / *rho[0]; }
          else {             *u_x[0] = 0.; }
          if( *rho[1] != 0.) { *u_x[1] = *u_x[1] / *rho[1]; }
          else {             *u_x[1] = 0.; }
        }
        else { *u_x[0] = 0.; *u_x[1] = 0.; }

        if( uy_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_y[0] = *u_y[0] / *rho[0]; }
          else {             *u_y[0] = 0.; }
          if( *rho[1] != 0.) { *u_y[1] = *u_y[1] / *rho[1]; }
          else {             *u_y[1] = 0.; }
        }
        else { *u_y[0] = 0.; *u_y[1] = 0.; }

        if( uz_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_z[0] = *u_z[0] / *rho[0]; }
          else {             *u_z[0] = 0.; }
          if( *rho[1] != 0.) { *u_z[1] = *u_z[1] / *rho[1]; }
          else {             *u_z[1] = 0.; }
        }
        else { *u_z[0] = 0.; *u_z[1] = 0.; }

      } /* if( !( is_solid)) */

    } /* for( n=0; n<lattice->NumNodes; n++) */

  } /* if( NUM_FLUID_COMPONENTS==2) */

  else if( NUM_FLUID_COMPONENTS == 1)
  {
    for( n=0; n<lattice->NumNodes; n++)
    {
      rho[0]      = &( lattice->macro_vars[0][n].rho);
      u_x[0]      = &( lattice->macro_vars[0][n].u[0]);
      u_y[0]      = &( lattice->macro_vars[0][n].u[1]);
      u_z[0]      = &( lattice->macro_vars[0][n].u[2]);
      is_solid    =  ( lattice->solids[0][n].is_solid);


      if( !( is_solid) )
      {
        if( *rho[0] != 0. && *u_x[0] != 0.)
        {
          *u_x[0] = *u_x[0] / *rho[0];
        }
        else
        {
          *u_x[0] = 0.;
        }

        if( *rho[0] != 0. && *u_y[0] != 0.)
        {
          *u_y[0] = *u_y[0] / *rho[0];
        }
        else
        {
          *u_y[0] = 0.;
        }

        if( *rho[0] != 0. && *u_z[0] != 0.)
        {
          *u_z[0] = *u_z[0] / *rho[0];
        }
        else
        {
          *u_z[0] = 0.;
        }

      } /* if( !( is_solid)) */

    } /* for( n=0; n<lattice->NumNodes; n++) */
  }
  else
  {
    printf(
        "compute_macro_vars() -- "
        "Unhandled case "
        "NUM_FLUID_COMPONENTS = %d . "
        "Exiting!\n",
        NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if PARALLEL

  rho_send_recv_begin(  lattice, 0);		//these calls deal with BOTH subs 0 and 1
  rho_send_recv_end(  lattice, 0);

  solid_send_recv_begin(  lattice, 0);		//these calls deal with subs 0, which is all that is required
  solid_send_recv_end(  lattice, 0);

#endif
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if NON_LOCAL_FORCES
  switch( NUM_FLUID_COMPONENTS)
  {
    case 1:
      {
        compute_phase_force( lattice, 0);
        compute_single_fluid_solid_force( lattice, 0);
        break;
      }

    case 2:
      {
        compute_fluid_fluid_force( lattice);
        compute_double_fluid_solid_force( lattice);
        break;
      }

    default:
      {
        printf("%s %d >> ERROR: Unhandled case %d\n",
            __FILE__,__LINE__,
            NUM_FLUID_COMPONENTS);
        process_exit(1);
        break;
      }
  }
#endif /* NON_LOCAL_FORCES */


} /* void compute_macro_vars( struct lattice_struct *lattice) */
#endif /* INAMURO_SIGMA_COMPONENT */

// void compute_feq( struct lattice_struct *lattice)
//##############################################################################
//
// C O M P U T E   F E Q
//
//  - Compute equilibrium distribution function, feq.
//
#if SAVE_MEMO

#else
void compute_feq( struct lattice_struct *lattice)
{
  int n, a;

  double *feq;
  int    subs;

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {

    for( n=0; n<lattice->NumNodes; n++)
    {
      feq       =    lattice->pdf[subs][n].feq;
      compute_single_feq(  lattice, n, subs, feq);

    }/*for( n=0; n<lattice->NumNodes; n++)*/
  }/*for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); ... */

} /* void compute_feq( struct lattice_struct *lattice) */
#endif

void compute_single_feq( struct lattice_struct *lattice, int n, int subs, double *feq)
{
  int  a;

  double tau;

  double W0,   W1,   W2;
  double ux,   uy,  uz,
         uxsq, uysq, uzsq, usq;
  double c; /*sound velocity*/
  double udotx;

  double rho, u_x, u_y, u_z;


#if STORE_UEQ
  double *ueq;
#endif /* STORE_UEQ */
  int    is_solid;



  /*Refers to Klaus Lglburger Thesis*/
  W0 = 1./3.;
  W1 = 1./18.;
  W2 = 1./36.;

      tau       =    lattice->param.tau[subs];
      rho       =    lattice->macro_vars[subs][n].rho;
      is_solid  =    lattice->solids[subs][n].is_solid;

//      if( /*fcountone*/  1 || !(is_solid))
      if( !lattice->time || !(is_solid))
      {
#if STORE_UEQ
        // Start with the composite macroscopic velocities.
        u_x       =    lattice->ueq[n].u[0];
        u_y       =    lattice->ueq[n].u[1];
        u_z       =    lattice->ueq[n].u[2];
        //  u_x       = *ueq;     *ueq++;
        //  u_y       = *ueq;     *ueq++;
        //  u_z       = *ueq;     *ueq++;
#else /* !( STORE_UEQ) */
        // Start with the individual component velocities.
        u_x       =    lattice->macro_vars[subs][n].u[0];
        u_y       =    lattice->macro_vars[subs][n].u[1];
        u_z       =    lattice->macro_vars[subs][n].u[2];
#endif /* STORE_UEQ */

        u_x       = u_x
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[0]
          +tau*lattice->force[ subs][n].sforce[0]
#endif /* NON_LOCAL_FORCES */
          +tau*lattice->param.gforce[ subs][0]
          ;
        u_y       = u_y
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[1]
          +tau*lattice->force[ subs][n].sforce[1]
#endif /* NON_LOCAL_FORCES */
          +tau*lattice->param.gforce[ subs][1]
          ;
        u_z       = u_z
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[2]
          +tau*lattice->force[ subs][n].sforce[2]
#endif /* NON_LOCAL_FORCES */
          + tau*lattice->param.gforce[ subs][2]
          ;

        usq = u_x*u_x + u_y*u_y + u_z*u_z;

        feq[0] = W0 * rho*(1. - 1.5*usq);

        for( a=1; a<=6; a++)
        {
          udotx = ((double)vx[a]*u_x+(double)vy[a]*u_y+(double)vz[a]*u_z);
          feq[a] = W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq);
        }

        for( a=7; a<Q; a++)
        {
          udotx = ((double)vx[a]*u_x+(double)vy[a]*u_y+(double)vz[a]*u_z);
          feq[a]  = W2*rho*(1. + 3.*udotx + 4.5*udotx*udotx - 1.5*usq);
        }

      }/*if(!is_solid)*/
      else
      {
        for( a=0; a<Q; a++)
        {
          feq[a] = 0.;
        }
      }

} /* void compute_single_feq( struct lattice_struct *lattice) */


#if NON_LOCAL_FORCES
void compute_phase_force( lattice_ptr lattice, int subs)
{
  double ***psi; //psi[LZ][LY][LX];
  double psi_temp;

  double rho;

  double *force;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk;

//printf("%s %d >> compute_phase_force() -- Hi!\n",
//    __FILE__,__LINE__);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  psi =    ( double***)malloc( nk*sizeof(double**));
  for( k=0; k<nk; k++)
  {
    psi[k] = ( double**)malloc( nj*sizeof(double*));
  }
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      psi[k][j] = ( double*)malloc( ni*sizeof(double));
    }
  }

  // Initialize psi.
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        psi[k][j][i] = 0.;
      }
    }
  }

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = XYZ2N( i, j, k, ni, nj);
        rho = lattice->macro_vars[subs][n].rho;

        if( rho!=0 && !( lattice->solids[subs][n].is_solid))
        {
          psi[k][j][i] = 4.*exp(-200./(rho));
        }
        else
        {
          psi[k][j][i] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = XYZ2N( i, j, k, ni, nj);
        force = lattice->force[subs][n].force;

        force[0] = 0.;
        force[1] = 0.;
        force[2] = 0.;

        if( !( lattice->solids[subs][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if( !( lattice->solids[subs][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              force[0] += WM*vx[a]*psi[ ka][ ja][ ia];
              force[1] += WM*vy[a]*psi[ ka][ ja][ ia];
              force[2] += WM*vz[a]*psi[ ka][ ja][ ia];
            }
          }

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if( !( lattice->solids[subs][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              force[0] += WD*vx[a]*psi[ ka][ ja][ ia];
              force[1] += WD*vy[a]*psi[ ka][ ja][ ia];
              force[2] += WD*vz[a]*psi[ ka][ ja][ ia];
            }
          }

          force[0] = -lattice->param.big_V0*psi[k][j][i]*( force[0]);
          force[1] = -lattice->param.big_V0*psi[k][j][i]*( force[1]);
          force[2] = -lattice->param.big_V0*psi[k][j][i]*( force[2]);

        } /* if( !( lattice->solids[subs][ j*ni+i].is_solid )) */

        else
        {
            *( force  ) = 0.;
            *( force+1) = 0.;
            *( force+2) = 0.;
        }

      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

 for( k=0; k<nk; k++)
 {
   for( j=0; j<nj; j++)
   {
     free( psi[k][j]);
   }
 }
 for( k=0; k<nk; k++)
 {
   free( psi[k]);
 }
 free( psi);

//printf("%s %d >> compute_phase_force() -- Bi!\n",
//    __FILE__,__LINE__);

} /* void compute_phase_force( lattice_ptr lattice) */

void compute_fluid_fluid_force( lattice_ptr lattice)
{
  double psi_temp;

  double rho, rho0, rho1;

  double *force[2];

  double sumx[2], sumy[2], sumz[2];

  int    a, subs;
  int    id;

  int    i,  j,  k, k1, k2,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n, na,
         ni, nj, nk;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

#if PARALLEL
        id = get_proc_id(lattice);
//  printf("%d\n",id);

  { k1 =0;   k2 = nk;}
  if(id ==0 && (!lattice->param.AllBoundaryPeriodic))
  { k1 =3;   k2 = nk;}
  if((id ==get_num_procs(lattice)-1) && (!lattice->param.AllBoundaryPeriodic))
  { k1 =0;   k2 = nk-2;}
//  printf("k1=%d, k2=%d, #########\n", k1, k2);
#else
  if(lattice->param.AllBoundaryPeriodic)
  { k1 =0;   k2 = nk;}
  else
  { k1 =3;   k2 = nk-2;}

#endif


  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
// for fluid 1 injection into porous media filled with fluid 2, rev_Huang
#if PARALLEL
//     rho_send_recv_begin( lattice, subs);
#endif
    for( k=k1; k<k2; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = XYZ2N( i, j, k, ni, nj);
          force[subs] = lattice->force[subs][n].force;

          force[subs][0] = 0.;
          force[subs][1] = 0.;
          force[subs][2] = 0.;

          if( !( lattice->solids[subs][ n].is_solid ))
          {
            for( a=1; a<=6; a++)
            {
              ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
              ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );

#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = XYZ2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
               {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( !( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_neg_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
// lattice->process.z_neg_solid_to_recv[n];
         }

        }
        else
        {
              if( !( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_pos_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
         }
        }

#else
        ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
        na = XYZ2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
              {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
              }

#endif
      }//a=0,... 6

            for( a=7; a<Q; a++)
            {
              ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
              ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = XYZ2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
               {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( !( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_neg_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
// lattice->process.z_neg_solid_to_recv[n];
         }

        }
        else
        {
              if( !( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_pos_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
         }
        }


#else
              ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
              na = XYZ2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
              {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
              }
#endif

            }//a=7....18

          } /* if( !( lattice->solids[subs][ n].is_solid )) */
        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */
    } /* for( k=0; k<nk; k++) */

#if PARALLEL
//    rho_send_recv_end( lattice, subs);
#endif
  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

    for( k=0; k<nk; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = XYZ2N( i, j, k, ni, nj);

          if( !( lattice->solids[0][ n].is_solid ))
          {
            force[0] = lattice->force[0][n].force;
            force[1] = lattice->force[1][n].force;
            rho0  = lattice->macro_vars[0][n].rho;
            rho1  = lattice->macro_vars[1][n].rho;

            sumx[/*subs*/0] = force[/*subs*/0][/*x*/0];
            sumx[/*subs*/1] = force[/*subs*/1][/*x*/0];

//            force[0][0] = -lattice->param.big_V0*rho0*sumx[1];
//            force[1][0] = -lattice->param.big_V0*rho1*sumx[0];
            force[0][0] = -lattice->param.big_V0*sumx[1];
            force[1][0] = -lattice->param.big_V0*sumx[0];

            sumy[/*subs*/0] = force[/*subs*/0][/*y*/1];
            sumy[/*subs*/1] = force[/*subs*/1][/*y*/1];

//            force[0][1] = -lattice->param.big_V0*rho0*sumy[1];
//            force[1][1] = -lattice->param.big_V0*rho1*sumy[0];
            force[0][1] = -lattice->param.big_V0*sumy[1];
            force[1][1] = -lattice->param.big_V0*sumy[0];

            sumz[/*subs*/0] = force[/*subs*/0][/*z*/2];
            sumz[/*subs*/1] = force[/*subs*/1][/*z*/2];

//            force[0][2] = -lattice->param.big_V0*rho0*sumz[1];
//            force[1][2] = -lattice->param.big_V0*rho1*sumz[0];
            force[0][2] = -lattice->param.big_V0*sumz[1];
            force[1][2] = -lattice->param.big_V0*sumz[0];

#if 0
            if( force[0][0]+ force[1][0]+ force[0][1]+ force[1][1]+ force[0][2]+ force[1][2] != 0.)
            {
printf("%s %d >> %f %f, %f %f, %f %f\n",
  __FILE__,__LINE__,
    force[0][0],
    force[1][0],
    force[0][1],
    force[1][1],
    force[0][2],
    force[1][2]
    );
            }
#endif

          } /* if( !( lattice->solids[subs][ n].is_solid )) */
        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */
    } /* for( k=0; k<nk; k++) */


} /* void compute_fluid_fluid_force( lattice_ptr lattice) */


void compute_single_fluid_solid_force( lattice_ptr lattice, int subs)
{
  double ***psi; //psi[ NUM_FLUID_COMPONENTS][LX][LY];
  double psi_temp;

  double rho;

  double *sforce;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk;

//printf("%s %d >> compute_single_fluid_solid_force() -- Hi!\n",
//    __FILE__,__LINE__);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  psi =    ( double***)malloc( nk*sizeof(double**));
  for( k=0; k<nk; k++)
  {
    psi[k] = ( double**)malloc( nj*sizeof(double*));
  }
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      psi[k][j] = ( double*)malloc( ni*sizeof(double));
    }
  }

  // Initialize psi.
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        psi[k][j][i] = 0.;
      }
    }
  }

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = XYZ2N( i, j, k, ni, nj);
        rho = lattice->macro_vars[subs][n].rho;

        if( rho!=0 && !( lattice->solids[subs][n].is_solid))
        {
          psi[k][j][i] = 4.*exp(-200./(rho));
        }
        else
        {
          psi[k][j][i] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = XYZ2N( i, j, k, ni, nj);
        sforce = lattice->force[subs][n].sforce;

        sforce[0] = 0.;
        sforce[1] = 0.;
        sforce[2] = 0.;

        if( !( lattice->solids[subs][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[subs][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sforce[0] += WM*vx[a];
              sforce[1] += WM*vy[a];
              sforce[2] += WM*vz[a];
            }
          }

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[subs][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sforce[0] += WD*vx[a];
              sforce[1] += WD*vy[a];
              sforce[2] += WD*vz[a];
            }
          }

          sforce[0] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[0]);
          sforce[1] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[1]);
          sforce[2] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[2]);

        } /* if( !( lattice->solids[subs][ j*ni+i].is_solid )) */

        else
        {
            *( sforce  ) = 0.;
            *( sforce+1) = 0.;
            *( sforce+2) = 0.;
        }

      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

 for( k=0; k<nk; k++)
 {
   for( j=0; j<nj; j++)
   {
     free( psi[k][j]);
   }
 }
 for( k=0; k<nk; k++)
 {
   free( psi[k]);
 }
 free( psi);

//printf("%s %d >> compute_single_fluid_solid_force() -- Bi!\n",
//    __FILE__,__LINE__);

} /* void compute_single_fluid_solid_sforce( lattice_ptr lattice, int subs) */

void compute_double_fluid_solid_force( lattice_ptr lattice)
{
  double rho;

  double *sforce[ /*NUM_FLUID_COMPONENTS*/ 2];

  double sum_x, sum_y, sum_z;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk, na;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
#if PARALLEL
//     solid_send_recv_begin( lattice, 0);
#endif
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = XYZ2N( i, j, k, ni, nj);
        sforce[0] = lattice->force[0][n].sforce;
        sforce[1] = lattice->force[1][n].sforce;

        sum_x = 0.;
        sum_y = 0.;
        sum_z = 0.;

        if( !( lattice->solids[0][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = XYZ2N( ia, ja, ka, ni, nj);
              if( ( lattice->solids[0][ na].is_solid ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( ( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
         }

        }
        else
        {
              if( ( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
         }
        }

#else
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[0][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sum_x += WM*vx[a];
              sum_y += WM*vy[a];
              sum_z += WM*vz[a];
            }
#endif

          } // a = 0....6

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = XYZ2N( ia, ja, ka, ni, nj);
              if( ( lattice->solids[0][ na].is_solid ))
               {
                sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( ( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
         }

        }
        else
        {
              if( ( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
          sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
         }
        }

#else
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[0][
                  XYZ2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sum_x += WD*vx[a];
              sum_y += WD*vy[a];
              sum_z += WD*vz[a];
            }
#endif
          } //a = 7,....18

#if 0
    if( lattice->param.big_V0_solid[0]+lattice->param.big_V0_solid[1] != 0.)
    {
       printf("%s %d >> %f %f\n",
         __FILE__,__LINE__,
         lattice->param.big_V0_solid[0],
         lattice->param.big_V0_solid[1] );
    }
#endif

//          sforce[0][0] = -lattice->param.big_V0_solid[0]*( sum_x) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
//          sforce[0][1] = -lattice->param.big_V0_solid[0]*( sum_y) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
//          sforce[0][2] = -lattice->param.big_V0_solid[0]*( sum_z) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
          sforce[0][0] = -lattice->param.big_V0_solid[0]*( sum_x)  *lattice->param.rhow;
          sforce[0][1] = -lattice->param.big_V0_solid[0]*( sum_y)  *lattice->param.rhow;
          sforce[0][2] = -lattice->param.big_V0_solid[0]*( sum_z)  *lattice->param.rhow;

//          sforce[1][0] = -lattice->param.big_V0_solid[1]*( sum_x) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
//          sforce[1][1] = -lattice->param.big_V0_solid[1]*( sum_y) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
//          sforce[1][2] = -lattice->param.big_V0_solid[1]*( sum_z) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
          sforce[1][0] = -lattice->param.big_V0_solid[1]*( sum_x)  *lattice->param.rhow;
          sforce[1][1] = -lattice->param.big_V0_solid[1]*( sum_y)  *lattice->param.rhow;
          sforce[1][2] = -lattice->param.big_V0_solid[1]*( sum_z)  *lattice->param.rhow;

        } /* if( !( lattice->solids[0][ j*ni+i].is_solid )) */

        else
        {
          sforce[0][0] = 0.;
          sforce[0][1] = 0.;
          sforce[0][2] = 0.;

          sforce[1][0] = 0.;
          sforce[1][1] = 0.;
          sforce[1][2] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

#if PARALLEL
//  solid_send_recv_end( lattice, 0);
#endif

} /* void compute_double_fluid_solid_force( lattice_ptr lattice) */

#endif /* NON_LOCAL_FORCES */

//void compute_max_rho( lattice_ptr lattice, double *max_rho, int subs)
void compute_max_rho( lattice_ptr lattice, double *max_rho, int subs)
{
  int n;
  double rho;

  *max_rho = 0.;

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho > *max_rho) { *max_rho = rho;}
    }
  }

  process_allreduce_double_max( lattice, max_rho);

} /* void compute_max_rho( lattice_ptr lattice, double *max_rho, int subs) */

//void compute_min_rho( lattice_ptr lattice, double *min_rho, int subs)
void compute_min_rho( lattice_ptr lattice, double *min_rho, int subs)
{
  int n;
  double rho;

  compute_max_rho( lattice, min_rho, subs);

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho < *min_rho) { *min_rho = rho;}
    }
  }

  process_allreduce_double_min( lattice, min_rho);

} /* void compute_min_rho( lattice_ptr lattice, double *min_rho, int subs) */

//void compute_ave_rho( lattice_ptr lattice, double *ave_rho, int subs)
void compute_ave_rho( lattice_ptr lattice, double *ave_rho, int subs)
{
  int n;
  int nn; // Number of non-solid nodes.
  *ave_rho = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes(lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_rho += fabs( get_rho( lattice, subs, n));
      nn++;
    }
  }

  process_allreduce_double_sum( lattice, ave_rho);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0) { *ave_rho = (*ave_rho) / nn;}

} /* void compute_ave_rho( lattice_ptr lattice, double *ave_rho, int subs) */

//void compute_max_u( lattice_ptr lattice, double *max_u, int subs)
void compute_max_u( lattice_ptr lattice, double *max_u, int subs)
{
  int n;
  double *max_ux = max_u+0;
  double *max_uy = max_u+1;
  double *max_uz = max_u+2;
  double u;
  *max_ux = 0.;
  *max_uy = 0.;
  *max_uz = 0.;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      u = fabs( get_ux( lattice, subs, n));
      if( u > *(max_ux)) { *max_ux = u;}

      u = fabs( get_uy( lattice, subs, n));
      if( u > *(max_uy)) { *max_uy= u;}

      u = fabs( get_uz( lattice, subs, n));
      if( u > *(max_uz)) { *max_uz= u;}
    }
  }

  process_allreduce_double_max( lattice, max_ux);
  process_allreduce_double_max( lattice, max_uy);
  process_allreduce_double_max( lattice, max_uz);

} /* void compute_max_u( lattice_ptr lattice, double *max_u, int subs) */

//void compute_min_u( lattice_ptr lattice, double *min_u, int subs)
void compute_min_u( lattice_ptr lattice, double *min_u, int subs)
{
  int n;
  double *min_ux = min_u+0;
  double *min_uy = min_u+1;
  double *min_uz = min_u+2;
  double u;

  compute_max_u( lattice, min_u,   subs);

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      u = fabs( get_ux( lattice, subs, n));
      if( u < *(min_ux)) { *min_ux = u;}

      u = fabs( get_uy( lattice, subs, n));
      if( u < *(min_uy)) { *min_uy = u;}

      u = fabs( get_uz( lattice, subs, n));
      if( u < *(min_uz)) { *min_uz = u;}
    }
  }

  process_allreduce_double_min( lattice, min_ux);
  process_allreduce_double_min( lattice, min_uy);
  process_allreduce_double_min( lattice, min_uz);

} /* void compute_min_u( lattice_ptr lattice, double *min_u) */

//void compute_ave_u( lattice_ptr lattice, double *ave_u, int subs)
void compute_ave_u( lattice_ptr lattice, double *ave_u, int subs)
{
  int n;
  int nn; // Number of non-solid nodes.
  double *ave_ux = ave_u+0;
  double *ave_uy = ave_u+1;
  double *ave_uz = ave_u+2;
  *ave_ux = 0.;
  *ave_uy = 0.;
  *ave_uz = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_ux += get_ux( lattice, subs, n);
      *ave_uy += get_uy( lattice, subs, n);
      *ave_uz += get_uz( lattice, subs, n);
      nn++;
    }
  }

  process_allreduce_double_sum( lattice, ave_ux);
  process_allreduce_double_sum( lattice, ave_uy);
  process_allreduce_double_sum( lattice, ave_uz);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_ux = (*ave_ux)/nn;
    *ave_uy = (*ave_uy)/nn;
    *ave_uz = (*ave_uz)/nn;
  }

} /* void compute_ave_u( lattice_ptr lattice, double *ave_u, int subs) */

void compute_flux( lattice_ptr lattice, double *flux, int subs)
{
  int n, nn;
  double rho, u_x, u_y, u_z;
  *(flux+0) = 0.;
  *(flux+1) = 0.;
  *(flux+2) = 0.;
  *(flux+3) = 0.;
  nn = 0;
  for( n=0; n<lattice->NumNodes; n++)
  {
    rho = lattice->macro_vars[subs][n].rho;
    u_x = lattice->macro_vars[subs][n].u[0];
    u_y = lattice->macro_vars[subs][n].u[1];
    u_z = lattice->macro_vars[subs][n].u[2];

    if( is_not_solid(lattice,subs))
    {
      *(flux+0) += rho*sqrt(u_x*u_x+u_y*u_y+u_z*u_z);
      *(flux+1) += rho*u_x;
      *(flux+2) += rho*u_y;
      *(flux+3) += rho*u_z;
      nn++;
    }
  }
  if( nn != 0)
  {
    *(flux+0) = (*(flux+0))/nn;
    *(flux+1) = (*(flux+1))/nn;
    *(flux+2) = (*(flux+2))/nn;
    *(flux+3) = (*(flux+3))/nn;
  }

  process_allreduce_double_sum( lattice, flux+0);
  process_allreduce_double_sum( lattice, flux+1);
  process_allreduce_double_sum( lattice, flux+2);
  process_allreduce_double_sum( lattice, flux+3);

  flux[0] /= get_num_procs(lattice);
  flux[1] /= get_num_procs(lattice);
  flux[2] /= get_num_procs(lattice);
  flux[3] /= get_num_procs(lattice);

} /* void compute_flux( lattice_ptr lattice, double *flux, int subs) */

#if STORE_UEQ
//void compute_max_ueq( lattice_ptr lattice, double *max_u)
void compute_max_ueq( lattice_ptr lattice, double *max_u)
{
  int n;
  *max_u = 0.;
  *(max_u+1) = 0.;
  *(max_u+2) = 0.;
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      if( fabs( lattice->ueq[n].u[0]) > *(max_u))
      {
        *max_u = fabs( lattice->ueq[n].u[0]);
      }
      if( fabs( lattice->ueq[n].u[1]) > *(max_u+1))
      {
        *(max_u+1) = fabs( lattice->ueq[n].u[1]);
      }
      if( fabs( lattice->ueq[n].u[2]) > *(max_u+2))
      {
        *(max_u+2) = fabs( lattice->ueq[n].u[2]);
      }
    }
  }

  process_allreduce_double_max( lattice, max_u+0);
  process_allreduce_double_max( lattice, max_u+1);
  process_allreduce_double_max( lattice, max_u+2);

} /* void compute_max_ueq( lattice_ptr lattice, double *max_u) */

//void compute_min_ueq( lattice_ptr lattice, double *min_u)
void compute_min_ueq( lattice_ptr lattice, double *min_u)
{
  int n;
  compute_max_ueq( lattice, min_u);
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      if( fabs( lattice->ueq[n].u[0]) < *(min_u))
      {
        *min_u = fabs( lattice->ueq[n].u[0]);
      }
      if( fabs( lattice->ueq[n].u[1]) < *(min_u+1))
      {
        *(min_u+1) = fabs( lattice->ueq[n].u[1]);
      }
      if( fabs( lattice->ueq[n].u[2]) < *(min_u+2))
      {
        *(min_u+2) = fabs( lattice->ueq[n].u[2]);
      }
    }
  }

  process_allreduce_double_min( lattice, min_u+0);
  process_allreduce_double_min( lattice, min_u+1);
  process_allreduce_double_min( lattice, min_u+2);

} /* void compute_min_ueq( lattice_ptr lattice, double *min_u) */

//void compute_ave_ueq( lattice_ptr lattice, double *ave_u)
void compute_ave_ueq( lattice_ptr lattice, double *ave_u)
{
  int n, nn;
  *ave_u = 0.;
  *(ave_u+1) = 0.;
  *(ave_u+2) = 0.;
  nn = 0;
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      *ave_u += fabs( lattice->ueq[n].u[0]);
      *(ave_u+1) += fabs( lattice->ueq[n].u[1]);
      *(ave_u+2) += fabs( lattice->ueq[n].u[2]);
      nn++;
    }
  }

  process_allreduce_double_sum( lattice, ave_u+0);
  process_allreduce_double_sum( lattice, ave_u+1);
  process_allreduce_double_sum( lattice, ave_u+2);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_u     = (*ave_u)/nn;
    *(ave_u+1) = (*(ave_u+1))/nn;
    *(ave_u+2) = (*(ave_u+2))/nn;
  }

} /* void compute_ave_ueq( lattice_ptr lattice, double *ave_u) */

#endif /* STORE_UEQ */

#if 0
//void compute_vorticity(
//       lattice_ptr lattice, int i, int j, int n, double *vor, int subs)
void compute_vorticity(
       lattice_ptr lattice, int i, int j, int n, double *vor, int subs)
{
  double duyx, duxy;
  int    nn[5]; // Indices of neighbors;

  int    ni=lattice->param.LX,
         nj=lattice->param.LY;

  int    ip, in,
         jp, jn;

  ip = ( i<ni-1)?(i+1):(0   );
  in = ( i>0   )?(i-1):(ni-1);
  jp = ( j<nj-1)?(j+1):(0   );
  jn = ( j>0   )?(j-1):(nj-1);

  nn[1] = j *ni + ip;
  nn[2] = jp*ni + i ;
  nn[3] = j *ni + in;
  nn[4] = jn*ni + i ;

  // NOTE: Assuming dx=1 .  TODO: Generalize dx?

  // Derivative of uy wrt x.
  if(    lattice->solids[subs][nn[1]].is_solid == BC_FLUID_NODE
      && lattice->solids[subs][nn[3]].is_solid == BC_FLUID_NODE)
  {
    // Forward difference.
    duyx = lattice->macro_vars[subs][nn[1]].u[1]
        - lattice->macro_vars[subs][n].u[1];
  }
//else if(    lattice->solids[subs][nn[3]].is_solid == BC_FLUID_NODE)
//{
//  // Backward difference.
//  duyx = lattice->macro_vars[subs][n].u[1]
//      - lattice->macro_vars[subs][nn[3]].u[1];
//}
  else
  {
    duyx = 0.;
  }

//printf("compute_vorticity() -- "
//    "n=%d, (i,j)=(%d,%d), duyx=%f, "
//    "nn1 = %d, nn3 = %d,"
//    "solids1 = %d, solids3 = %d"
//    "\n",
//    n, i, j, duyx,
//    nn[1], nn[3],
//    lattice->solids[subs][nn[1]].is_solid,
//    lattice->solids[subs][nn[3]].is_solid);

  // Derivative of ux wrt y.
  if(    lattice->solids[subs][nn[2]].is_solid == BC_FLUID_NODE
      && lattice->solids[subs][nn[4]].is_solid == BC_FLUID_NODE)
  {
    // Forward difference.
    duxy = lattice->macro_vars[subs][nn[2]].u[0]
        - lattice->macro_vars[subs][n].u[0];
  }
//else if(    lattice->solids[subs][nn[4]].is_solid == BC_FLUID_NODE)
//{
//  // Backward difference.
//  duxy = lattice->macro_vars[subs][n].u[0]
//      - lattice->macro_vars[subs][nn[4]].u[0];
//}
  else
  {
    duxy = 0.;
  }

  if( duxy*duyx != 0.)
  {
    *vor = duyx - duxy;
  }
  else
  {
    *vor = 0.;
  }

} /* void compute_vorticity( lattice_ptr lattice, int i, int j, int n, ... */

//void compute_max_vor(
//       lattice_ptr lattice,double *max_vor_p,double *max_vor_n, int subs)
void compute_max_vor(
       lattice_ptr lattice, double *max_vor_p, double *max_vor_n, int subs)
{
  int n;
  double vor;
  int nnz;

  *max_vor_p = 0.;
  *max_vor_n = 0.;
  nnz = 0;

  for( n=0; n<=lattice->NumNodes; n++)
  {
    if(    lattice->solids[subs][n].is_solid == BC_FLUID_NODE)
    {
      compute_vorticity( lattice,
                         n%lattice->param.LX,
                         n/lattice->param.LX,
                         n,
                         &vor,
                         subs  );
      if( vor != 0.) { nnz++;}

      if( vor > *max_vor_p)
      {
        *max_vor_p = vor;

      } /* if( vor > *max_vor_p) */

      else if( vor < *max_vor_n)
      {
        *max_vor_n = vor;

      } /* if( vor > *max_vor_p) */

    } /* if( lattice->solids[subs][n].is_solid == 0) */

  } /* for( n=0; n<=lattice->NumNodes; n++) */

#if 0 && VERBOSITY_LEVEL > 0
  printf("compute_max_vor() -- nnz = %d.  nnz/NumNodes = %f\n",
    nnz, (double)nnz/(double)lattice->NumNodes);
#endif /* 0 && VERBOSITY_LEVEL > 0 */

} /* void compute_max_vor( lattice_ptr lattice, double *max_vor_p, ... */

//void compute_ave_vor(
//       lattice_ptr lattice,double *ave_vor_p,double *ave_vor_n, int subs)
void compute_ave_vor(
       lattice_ptr lattice, double *ave_vor_p, double *ave_vor_n, int subs)
{
  int n;
  double vor;
  int nnz;
  int num_p;
  int num_n;

  *ave_vor_p = 0.;
  *ave_vor_n = 0.;
  nnz = 0;
  num_p = 0;
  num_n = 0;

  for( n=0; n<=lattice->NumNodes; n++)
  {
    if(    lattice->solids[subs][n].is_solid == BC_FLUID_NODE)
    {
      compute_vorticity( lattice,
                         n%lattice->param.LX,
                         n/lattice->param.LX,
                         n,
                         &vor,
                         subs  );
      if( vor != 0.) { nnz++;}

      if( vor > *ave_vor_p)
      {
        *ave_vor_p += vor;
        num_p++;

      } /* if( vor > *ave_vor_p) */

      else if( vor < *ave_vor_n)
      {
        *ave_vor_n += vor;
        num_n++;

      } /* if( vor > *ave_vor_p) */

    } /* if( lattice->solids[subs][n].is_solid == 0) */

  } /* for( n=0; n<=lattice->NumNodes; n++) */

  if( num_p > 0) { *ave_vor_p /= num_p;}
  if( num_n > 0) { *ave_vor_n /= num_n;}

#if 0 && VERBOSITY_LEVEL > 0
  printf("compute_ave_vor() -- nnz = %d.  nnz/NumNodes = %f\n",
    nnz, (double)nnz/(double)lattice->NumNodes);
#endif /* VERBOSITY_LEVEL > 0 */

} /* void compute_ave_vor( lattice_ptr lattice, double *ave_vor_p, ... */
#endif


