//##############################################################################
//
// all_nodes_latman.c
//
//  - Lattice Manager.
//
//  - Routines for managing a lattice:
//
//    - construct
//    - init
//    - destruct
//
//  - This file, with prefix "all_nodes_", is for the version of the code
//    that stores all nodes of the domain even if they are interior solid
//    nodes that are not involved in any computations.  This is more
//    efficient, in spite of storing unused nodes, if the ratio of
//    interior solid nodes to total nodes is sufficiently low.  How
//    low is sufficient? is a difficult question.  If storage is the
//    only consideration, then the two approaches balance at somewhere
//    around .75 .  But more observations need to be made to characterize
//    the trade-off in terms of computational efficiency.
//

// void construct_lattice( struct lattice_struct *lattice)
//##############################################################################
//
// C O N S T R U C T   L A T T I C E 
//
//  - Construct lattice.
//
void construct_lattice( lattice_ptr *lattice, int argc, char **argv)
{
  // Variable declarations
  int    **matrix;
  int    i, 
         j;
  int    n;
  int    subs;
  int    width, 
         height;
  char   filename[1024];
  char   dirname[1024];

  // Allocate the lattice structure.
  *lattice = ( struct lattice_struct*)malloc( sizeof(struct lattice_struct));

  assert(*lattice!=NULL);

  process_init( *lattice, argc, argv);

  // Read problem parameters
  if( argc == 2)
  {
    sprintf( dirname, "in_%s", argv[1]);
  }
  else
  {
    sprintf( dirname, "in");
  }
  sprintf(filename, "./%s/%s", dirname, "params.in");
  printf("%s %d %04d >> Reading params from file \"%s\"\n", 
    __FILE__, __LINE__, get_proc_id( *lattice), filename);
  read_params( *lattice, filename);

  process_compute_local_params( *lattice);

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  // Allocate space for solids.
  (*lattice)->solids[subs] = 
    (struct solids_struct*)malloc( 
       (*lattice)->NumNodes*sizeof(struct solids_struct));
 }

  // Read solids from .raw file.
  sprintf(filename, "./in/%dx%dx%d.raw", 
      get_g_LX(*lattice),
      get_g_LY(*lattice),
      get_g_LZ(*lattice) );
  printf("%s %d >> Reading solids from file \"%s\"\n", 
    __FILE__, __LINE__, filename);
  read_solids( *lattice, filename);

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  // Allocate NumNodes particle distribution functions.
  (*lattice)->pdf[subs] = 
    ( struct pdf_struct*)malloc( 
        (*lattice)->NumNodes*sizeof( struct pdf_struct));
  if( (*lattice)->pdf[subs] == NULL)
  {
    printf(
      "%s %d %04d >> "
      "construct_lattice() -- ERROR:  "
      "Attempt to allocate %d struct pdf_struct types failed.  "
      "Exiting!\n",
      __FILE__,__LINE__,get_proc_id(*lattice), (*lattice)->NumNodes );
    process_exit(1);
  }

  // Allocate NumNodes macroscopic variables.
  (*lattice)->macro_vars[subs] = 
    ( struct macro_vars_struct*)malloc( 
        (*lattice)->NumNodes*sizeof( struct macro_vars_struct));
  if( (*lattice)->macro_vars[subs]==NULL)
  {
    printf(
      "construct_lattice() -- ERROR:  "
      "Attempt to allocate %d struct macro_vars_struct types failed.  "
      "Exiting!\n",
      (*lattice)->NumNodes
      );
    process_exit(1);
  }

#if NON_LOCAL_FORCES
  // Allocate NumNodes elements for force.
  (*lattice)->force[subs] = 
    ( struct force_struct*)malloc( 
        (*lattice)->NumNodes*sizeof( struct force_struct));
  if( (*lattice)->force[subs]==NULL)
  {
    printf(
      "construct_lattice() -- ERROR:  "
      "Attempt to allocate %d struct force_struct types failed.  "
      "Exiting!\n",
      (*lattice)->NumNodes
      );
    process_exit(1);
  }
#endif /* NON_LOCAL_FORCES */
 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if STORE_UEQ
  // Allocate NumNodes elements for ueq.
  (*lattice)->ueq = 
    ( struct ueq_struct*)malloc( 
        (*lattice)->NumNodes*sizeof( struct ueq_struct));
  if( (*lattice)->ueq==NULL)
  {
    printf(
      "construct_lattice() -- ERROR:  "
      "Attempt to allocate %d struct ueq_struct types failed.  "
      "Exiting!\n",
      (*lattice)->NumNodes
      );
    process_exit(1);
  }
#endif /* STORE_UEQ */

  dump_params( *lattice);

} /* void construct_lattice( struct lattice_struct **lattice) */

// void init_problem( struct lattice_struct *lattice)
//##############################################################################
//
// I N I T   P R O B L E M 
//
//  - Initialize the problem on the lattice.
//
//    - Set the initial density and velocity.
//
//    - Compute the initial feq.
//
void init_problem( struct lattice_struct *lattice)
{
  int    a, n, i, j, k, ni, nj, nk, id1, id2, newz1, newz2, k1,k2;
  double x, u_max, K, drho, m;
  double *macro_var_ptr;
  double *f, *feq, *ftemp, ftemp2[19];
  double *rho, *u_x, *u_y, *u_z;

  double fcountone; 

#if STORE_UEQ
  double *ueq_x, *ueq_y, *ueq_z;
#endif /* STORE_UEQ */
#if NON_LOCAL_FORCES
  double *force;
#endif /* NON_LOCAL_FORCES */
  int    is_solid;
  int    subs;
  double kappa;
  double ti;
  double y;

#if VERBOSITY_LEVEL > 0
  printf("init_problem() -- Initilizing problem...\n");
#endif /* VERBOSITY_LEVEL > 0 */

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  lattice->time = 0;

#if PARALLEL
 if(lattice->param.initial_condition ==1)
 {
 
 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {

  for( n=0; n<get_proc_id(lattice)* lattice->NumNodes; n++)
  {
 
      if( (double)rand()/(double)RAND_MAX < lattice->param.cut)
              {
                //Non to test RAND is not a real random;
              }
  }

 }

 }//it works !! Great!!, no problem now , it gives same results as PC version!! 
#endif


 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {

  for( n=0; n<lattice->NumNodes; n++)
  {
    i = N2X(n,ni,nj,nk);//n%lattice->param.LX;
    j = N2Y(n,ni,nj,nk);//n/lattice->param.LX;
    k = N2Z(n,ni,nj,nk);

    rho = &( lattice->macro_vars[subs][n].rho);
    u_x = &( lattice->macro_vars[subs][n].u[0]);
    u_y = &( lattice->macro_vars[subs][n].u[1]);
    u_z = &( lattice->macro_vars[subs][n].u[2]);
#if STORE_UEQ
    ueq_x = &( lattice->ueq[n].u[0]);
    ueq_y = &( lattice->ueq[n].u[1]);
    ueq_z = &( lattice->ueq[n].u[2]);
#endif /* STORE_UEQ */


    is_solid = lattice->solids[subs][n].is_solid;

    // Set initial density.
    if( ( 1 || !( is_solid)) )
//    if(   !( is_solid) )
//rev_Huang
    {
      switch( lattice->param.initial_condition)
      {
		
        case 0: // Uniform
        {
           switch( NUM_FLUID_COMPONENTS)
	   { 
		case 1:  
		{
		*rho = lattice->param.rho_A[subs];
                break;
		}
		
		case 2:  
		{
		*rho = lattice->param.rho_B[subs];
		break;
		}
	   }
		break;
        }

        case 1: // Random
        {
          switch( NUM_FLUID_COMPONENTS)
          {
            case 1:
            {
              *rho = 
                lattice->param.rho_A[0] + 10.*(double)rand()/(double)RAND_MAX;
              break;
            }

            case 2:
            {
              if( (double)rand()/(double)RAND_MAX < lattice->param.cut)
              {
                *rho = lattice->param.rho_A[subs];
              }
              else
              {
                *rho = lattice->param.rho_B[subs];
              }
              break;
            }

            default:
            {
              printf("%s %d >> ERROR: Unhandled case %d.\n",
                __FILE__,__LINE__,
                NUM_FLUID_COMPONENTS);
              process_exit(1);
              break;
            }
          }
          break;
        }

        case 2: // Cube
        {
#if PARALLEL
#if 0
		id1 =  (int)floor(lattice->param.z1/(double)(nk);
		newz1 =  (int)(lattice->param.z1) % (nk); 
		id2 =  (int)floor(lattice->param.z2/(double)(nk));
		newz2 =  (int)(lattice->param.z2) % (nk);
		if(id2>get_num_procs(lattice)) 
		{
			id2 = get_num_procs(lattice);
			newz2 = nk ;
		}
#endif
#endif
          switch( NUM_FLUID_COMPONENTS)
          {
            case 1:
            {
#if PARALLEL
#if 0
	     if(get_proc_id(lattice) ==id1 ) {k1=newz1; k2 = nk;}
	     if(get_proc_id(lattice) ==id2 ) {k1=0; k2 = newz2;}
	     if(get_proc_id(lattice) >id1 && get_proc_id(lattice) <id2 ) {k1=0; k2 = nk;}
#endif
#endif	     
             if( i > lattice->param.x1 && i < lattice->param.x2 &&
                 j > lattice->param.y1 && j < lattice->param.y2 && 
#if PARALLEL
                 k+ get_proc_id(lattice)*nk > lattice->param.z1 && k+ get_proc_id(lattice)*nk < lattice->param.z2 )
//if 0		  k >=k1 && k<=k2 )
#else		 
		 k > lattice->param.z1 && k < lattice->param.z2  )
#endif
              {
                *rho = lattice->param.rho_A[0];
              }
              else
              {
                *rho = lattice->param.rho_B[0];
              }
              break;
            }

            case 2:
            {
#if PARALLEL
#if 0		    
	     if(get_proc_id(lattice) ==id1 ) {k1=newz1; k2 = nk;}
	     if(get_proc_id(lattice) ==id2 ) {k1=0; k2 = newz2;}
	     if(get_proc_id(lattice) >id1 && get_proc_id(lattice) <id2 ) {k1=0; k2 = nk;}
#endif
#endif
	     if( i >= lattice->param.x1 && i <= lattice->param.x2 &&
                 j >= lattice->param.y1 && j <= lattice->param.y2 && 
#if PARALLEL
		k+ get_proc_id(lattice)*nk > lattice->param.z1 && k+ get_proc_id(lattice)*nk < lattice->param.z2)
//		  k >=k1 && k<=k2 )
#else		 
		 k > lattice->param.z1 && k < lattice->param.z2  )
#endif
              {
                *rho = lattice->param.rho_A[subs];
              }
              else
              {
                *rho = lattice->param.rho_B[subs];
              }
              break;
            }

            default:
            {
              printf("%s %d >> ERROR: Unhandled case %d.\n",
                __FILE__,__LINE__,
                NUM_FLUID_COMPONENTS);
              process_exit(1);
              break;
            }
          }

          break;
        }

        case 3: // Sphere
        {
          switch( NUM_FLUID_COMPONENTS)
          {
            case 1:
            {
              if( (i-lattice->param.x0)*(i-lattice->param.x0)+(j-lattice->param.y0)*
		  (j-lattice->param.y0)+(k-lattice->param.z0)*(k-lattice->param.z0)<= lattice->param.r0*lattice->param.r0 )
              {
                *rho = lattice->param.rho_A[0];
              }
              else
              {
                *rho = lattice->param.rho_B[0];
              }
              break;
            }
            case 2:
            {
              if( (i-lattice->param.x0)*(i-lattice->param.x0)+(j-lattice->param.y0)*
		  (j-lattice->param.y0)+(k-lattice->param.z0)*(k-lattice->param.z0)<= lattice->param.r0*lattice->param.r0 )
              {
                *rho = lattice->param.rho_A[subs];
              }
              else
              {
                *rho = lattice->param.rho_B[subs];
              }
              break;
            }
            default:
            {
              printf( 
                "%s (%d) >> init_problem() -- Unhandled case  "
                "NUM_FLUID_COMPONENTS = %d.  "
                "Exiting!\n", __FILE__, __LINE__,
                NUM_FLUID_COMPONENTS);
              process_exit(1);
            }
          }
          break;
        }
      case 4: // Linear Density Gradient
        {
          switch( NUM_FLUID_COMPONENTS)
          {
            case 1:
            {                       
              
                *rho = ((lattice->param.rho_out-lattice->param.rho_in)/nk)*k + lattice->param.rho_in;//lattice->param.rho_A[0];
              
              
              break;
            }

            case 2:
            {
              if( i > 20 && i < 32 &&
                  j > 10 && j < 22
               && k >  8 && k < 18     
                  )
              {
                *rho = lattice->param.rho_A[subs];
              }
              else
              {
                *rho = lattice->param.rho_B[subs];
              }
              break;
            }

            default:
            {
              printf("%s %d >> ERROR: Unhandled case %d.\n",
                __FILE__,__LINE__,
                NUM_FLUID_COMPONENTS);
              process_exit(1);
              break;
            }
          }

          break;
        }

        default:
        {
          printf( 
            "%s (%d) >> init_problem() -- Unhandled case  "
            "lattice->param.initial_condition = %d.  "
            "Exiting!\n", __FILE__, __LINE__,
            lattice->param.initial_condition );
          process_exit(1);
          break;
        }
      } /* switch( lattice->param.initial_condition) */
    } /* if( ( 1 || !( solids->is_solid & BC_SOLID_NODE)) ) */
    else
    {
    //if( solids->is_solid)
    //{
    //  *macro_var_ptr++ = lattice->param.rho_A[subs];
    //  *macro_var_ptr++ = lattice->param.rho_in;
    //}
    //else
    //{
        *rho = 0.;
    //}
    } /* if( ( 1 || !( solids->is_solid & BC_SOLID_NODE)) ) else */

#if SPONGE
	     if(// i >= lattice->param.x1 && i <= lattice->param.x2 &&
                // j >= lattice->param.y1 && j <= lattice->param.y2 && 
#if PARALLEL
		k+ get_proc_id(lattice)*nk < lattice->param.z1 || k+ get_proc_id(lattice)*nk > lattice->param.z2)
//		  k >=k1 && k<=k2 )
#else		 
		 k < lattice->param.z1 || k > lattice->param.z2  )
#endif
              {
                *rho = lattice->param.rho_A[subs];
              }
#endif

    // Set initial velocity.
    *u_x = 0.;
    *u_y = 0.;
    *u_z = 0.;//lattice->param.ux_in; //0.;

#if STORE_UEQ
    *ueq_x = 0.;
    *ueq_y = 0.;
    *ueq_z = 0.;
#endif /* STORE_UEQ */

  } /* for( n=0; n<lattice->NumNodes; n++) */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

printf("%s %d >> Before compute_feq...\n",__FILE__,__LINE__);

  // Compute initial feq.
#if SAVE_MEMO

#else
compute_feq( lattice);
#endif

printf("%s %d >> After compute_feq!\n",__FILE__,__LINE__);

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  for( n=0; n<lattice->NumNodes; n++)
  {
#if SAVE_MEMO
    compute_single_feq(  lattice, n, subs, &(lattice->pdf[subs][n].ftemp) );
    
#else

      	  f = lattice->pdf[subs][n].f;
    feq = lattice->pdf[subs][n].feq;
    ftemp = lattice->pdf[subs][n].ftemp;
    is_solid = lattice->solids[subs][n].is_solid;
//***********************************************************************************
      if(!lattice->time)
	  {
		  fcountone = is_solid;
	  }
	  else
	  {
		  fcountone = 0.;
	  }
//***************************************************************p-edit*************
    if( /*ALLOW_INIT_ON_SOLID_NODES*/fcountone || !( is_solid))
    {
      // Copy feq to f.
      for( a=0; a<Q; a++)
      {
		 // printf("feq = %d, is_solid = %d, ftemp = %d, f[ %d ] = %d \n", lattice->pdf[subs][n].feq,
		 // lattice->solids[subs][n].is_solid,lattice->pdf[subs][n].ftemp, a, f[a]); //****Problem here?******************
        f[a] = feq[a];
      }

      // Initialize ftemp.
      for( a=0; a<Q; a++)
      {
        ftemp[a]= 0.;
      }

    }
    else
    {
      // f
      for( a=0; a<Q; a++)
      {
        f[a] = 0.;
      }
          
      // ftemp
      for( a=0; a<Q; a++)
      {
        ftemp[a] = 0.;
      }
   
    }
#endif
  } /* for( n=0; n<lattice->NumNodes; n++) */

#if NON_LOCAL_FORCES
  for( n=0; n<lattice->NumNodes; n++)
  {
    lattice->force[subs][n].force[0] = 0.;
    lattice->force[subs][n].force[1] = 0.;
    lattice->force[subs][n].force[2] = 0.;
    lattice->force[subs][n].sforce[0] = 0.;
    lattice->force[subs][n].sforce[1] = 0.;
    lattice->force[subs][n].sforce[2] = 0.;
  }
#endif /* NON_LOCAL_FORCES */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */


  //compute_macro_vars( lattice);

  //dump_macro_vars( lattice, /*time=*/ 0);


#if VERBOSITY_LEVEL > 0
  printf("init_problem() -- Problem initialized.\n");
#endif /* VERBOSITY_LEVEL > 0 */

//compute_macro_vars( lattice);//*********************************************p-edit***

} /* void init_problem( struct lattice_struct *lattice) */

// void destruct_lattice( struct lattice_struct *lattice)
//##############################################################################
//
// D E S T R U C T   L A T T I C E 
//
//  - Destruct lattice.
//
void destruct_lattice( struct lattice_struct *lattice)
{
  int subs;

  assert( lattice!=NULL);

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  assert( lattice->pdf[subs]!=NULL);
  free(   lattice->pdf[subs]);

  assert( lattice->macro_vars[subs]!=NULL);
  free(   lattice->macro_vars[subs]);

  assert( lattice->solids[subs]!=NULL);
  free(   lattice->solids[subs]);
#if NON_LOCAL_FORCES
  assert( lattice->force[subs]!=NULL);
  free(   lattice->force[subs]);
#endif /* NON_LOCAL_FORCES */
 }

#if STORE_UEQ
  free(   lattice->ueq);
#endif /* STORE_UEQ */

  process_finalize();

  free(   lattice);

} /* void destruct_lattice( struct lattice_struct *lattice) */

//##############################################################################
void dump_north_pointing_pdfs(
       lattice_ptr lattice,
       const int subs,
       const int z_slice,
       char *comment_str,
       const int which_pdf)
{
  int n, p;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  double pdfs[5*get_LX(lattice)*get_LY(lattice)];

  if( domain_is_too_big_to_display( lattice)) { return;}

 for( p=0; p<get_num_procs( lattice); p++)
 {
 process_barrier();
 if( p == get_proc_id(lattice))
 {
  if( z_slice >= 0)
  {
    k = z_slice;

    gather_north_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

    printf("\n\n// Proc %d, k=%d >> North pointing, \"%s\".",
      get_proc_id( lattice), k, comment_str);
    printf("\n ");
    for( i=0; i<ni; i++)
    {
      printf("+");
      printf("---");
      printf("---");
      printf("---");
      printf("-");
    }
    printf("+");

    for( j=0; j<nj; j++)
    {
      // 0 1 2 3 4
      // O W E N S
    
      // South
      n = 5*j*ni + 4;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf("   ");
        printf(" %2.0f", pdfs[n]);
        printf("   ");
        printf(" ");
        n+=5;
      }
      printf("|");
    
      // West/O/East
      n = 5*j*ni + 0;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf(" %2.0f", pdfs[n+1]);
        printf(" %2.0f", pdfs[n]);
        printf(" %2.0f", pdfs[n+2]);
        printf(" ");
        n+=5;
      }
      printf("|");
    
      // North
      n = 5*j*ni + 3;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf("   ");
        printf(" %2.0f", pdfs[n]);
        printf("   ");
        printf(" ");
        n+=5;
      }
      printf("|");
    
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("+");
        printf("---");
        printf("---");
        printf("---");
        printf("-");
      }
      printf("+");
    
    } /* if( j=0; j<nj; j++) */

  } /* if( z_slice >= 0) */

  else
  {
    for( k=0; k<get_LZ(lattice); k++)
    {

      gather_north_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

      printf("\n\n// Proc %d, k=%d >> North pointing, \"%s\".",
        get_proc_id( lattice), k, comment_str);
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("+");
        printf("---");
        printf("---");
        printf("---");
        printf("-");
      }
      printf("+");

      for( j=0; j<nj; j++)
      {
        // 0 1 2 3 4
        // O W E N S
      
        // South
        n = 5*j*ni + 4;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf("   ");
          printf(" %2.0f", pdfs[n]);
          printf("   ");
          printf(" ");
          n+=5;
        }
        printf("|");
      
        // West/O/East
        n = 5*j*ni + 0;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf(" %2.0f", pdfs[n+1]);
          printf(" %2.0f", pdfs[n]);
          printf(" %2.0f", pdfs[n+2]);
          printf(" ");
          n+=5;
        }
        printf("|");
      
        // North
        n = 5*j*ni + 3;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf("   ");
          printf(" %2.0f", pdfs[n]);
          printf("   ");
          printf(" ");
          n+=5;
        }
        printf("|");
      
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("+");
          printf("---");
          printf("---");
          printf("---");
          printf("-");
        }
        printf("+");
      
      } /* if( j=0; j<nj; j++) */

    }

  } /* if( z_slice >= 0) else */

 }

 } /* for( p=0; p<get_num_procs( lattice); p++) */
 process_barrier();

} /* void dump_north_pointing_pdfs( lattice_ptr lattice, const int subs, ... */

//##############################################################################
void dump_south_pointing_pdfs(
       lattice_ptr lattice,
       const int subs,
       const int z_slice,
       char *comment_str,
       const int which_pdf)
{
  int n, p;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  double pdfs[5*get_LX(lattice)*get_LY(lattice)];

  if( domain_is_too_big_to_display( lattice)) { return;}

 for( p=0; p<get_num_procs( lattice); p++)
 {
 process_barrier();
 if( p == get_proc_id(lattice))
 {
  if( z_slice >= 0)
  {
    k = z_slice;

    gather_south_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

    printf("\n\n// Proc %d, k=%d >> South pointing, \"%s\".",
      get_proc_id( lattice), k, comment_str);
    printf("\n ");
    for( i=0; i<ni; i++)
    {
      printf("+");
      printf("---");
      printf("---");
      printf("---");
      printf("-");
    }
    printf("+");

    for( j=0; j<nj; j++)
    {
      // 0 1 2 3 4
      // O W E N S
    
      // South
      n = 5*j*ni + 4;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf("   ");
        printf(" %2.0f", pdfs[n]);
        printf("   ");
        printf(" ");
        n+=5;
      }
      printf("|");
    
      // West/O/East
      n = 5*j*ni + 0;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf(" %2.0f", pdfs[n+1]);
        printf(" %2.0f", pdfs[n]);
        printf(" %2.0f", pdfs[n+2]);
        printf(" ");
        n+=5;
      }
      printf("|");
    
      // North
      n = 5*j*ni + 3;
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("|");
        printf("   ");
        printf(" %2.0f", pdfs[n]);
        printf("   ");
        printf(" ");
        n+=5;
      }
      printf("|");
    
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("+");
        printf("---");
        printf("---");
        printf("---");
        printf("-");
      }
      printf("+");
    
    } /* if( j=0; j<nj; j++) */

  } /* if( z_slice >= 0) */

  else
  {
    for( k=0; k<get_LZ(lattice); k++)
    {

      gather_south_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

      printf("\n\n// Proc %d, k=%d >> South pointing, \"%s\".",
        get_proc_id( lattice), k, comment_str);
      printf("\n ");
      for( i=0; i<ni; i++)
      {
        printf("+");
        printf("---");
        printf("---");
        printf("---");
        printf("-");
      }
      printf("+");

      for( j=0; j<nj; j++)
      {
        // 0 1 2 3 4
        // O W E N S
      
        // South
        n = 5*j*ni + 4;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf("   ");
          printf(" %2.0f", pdfs[n]);
          printf("   ");
          printf(" ");
          n+=5;
        }
        printf("|");
      
        // West/O/East
        n = 5*j*ni + 0;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf(" %2.0f", pdfs[n+1]);
          printf(" %2.0f", pdfs[n]);
          printf(" %2.0f", pdfs[n+2]);
          printf(" ");
          n+=5;
        }
        printf("|");
      
        // North
        n = 5*j*ni + 3;
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("|");
          printf("   ");
          printf(" %2.0f", pdfs[n]);
          printf("   ");
          printf(" ");
          n+=5;
        }
        printf("|");
      
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("+");
          printf("---");
          printf("---");
          printf("---");
          printf("-");
        }
        printf("+");
      
      } /* if( j=0; j<nj; j++) */

    }

  } /* if( z_slice >= 0) else */

 }

 } /* for( p=0; p<get_num_procs( lattice); p++) */
 process_barrier();

} /* void dump_south_pointing_pdfs( lattice_ptr lattice, const int subs, ... */

int domain_is_too_big_to_display( lattice_ptr lattice)
{
  if(   get_LX(lattice) > 12
     || get_LX(lattice) > 12
     || get_LX(lattice) > 12 ) { return 1;} else { return 0;}
}

int domain_is_not_too_big_to_display( lattice_ptr lattice)
{
  return !domain_is_too_big_to_display(lattice);
}

void display_warning_about_contrived_data( lattice_ptr lattice)
{
  printf("\n");
  printf("\n");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("\n");
  printf("\n");
  printf(" W A R N I N G:");
  printf("\n");
  printf("\n");
  printf("      Contrived data! This run good for DEBUG purposes only.");
  printf("\n");
  printf("\n");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("\n");
  printf("\n");
}

// int get_sizeof_lattice_structure( lattice_ptr lattice)
//##############################################################################
//
// G E T _ S I Z E O F _ L A T T I C E _ S T R U C T U R E
//
//  - Return size of struct lattice_struct in bytes.
//
int get_sizeof_lattice_structure( lattice_ptr lattice)
{
  return sizeof( struct lattice_struct);

} /* int get_sizeof_lattice_structure( lattice_ptr lattice) */

// int get_sizeof_lattice( lattice_ptr lattice)
//##############################################################################
//
// G E T _ S I Z E O F _ L A T T I C E
//
//  - Return size of lattice in bytes.
//
int get_sizeof_lattice( lattice_ptr lattice)
{
  return
      sizeof(int)
    + sizeof(int)
    + sizeof(int)
    + sizeof(int)*NUM_FLUID_COMPONENTS
    + sizeof(int)*NUM_FLUID_COMPONENTS

    + sizeof(struct param_struct)

    + lattice->NumNodes
    * ( 
        NUM_FLUID_COMPONENTS*sizeof(struct pdf_struct)
      + NUM_FLUID_COMPONENTS*sizeof(struct macro_vars_struct)
      + NUM_FLUID_COMPONENTS*sizeof(struct solids_struct) 
#if NON_LOCAL_FORCES
      + NUM_FLUID_COMPONENTS*sizeof(struct force_struct) 
#endif /* NON_LOCAL_FORCES */
#if STORE_UEQ
      + sizeof(struct ueq_struct) 
#endif /* STORE_UEQ */
#if POROUS_MEDIA
      + sizeof(struct ns_struct) 
#endif /* POROUS_MEDIA */
      );

} /* int get_sizeof_lattice( lattice_ptr lattice) */

