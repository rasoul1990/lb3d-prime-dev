//##############################################################################
//
// collide.c
//

//##############################################################################
//
// void collide( lattice_ptr lattice)
#if SAVE_MEMO
void collide( lattice_ptr lattice)
{
  double *feq;
  double f;
  double *ftemp;

  double omega;

  int    is_solid;

  int    n, a;

  int    subs;

   feq= (double*)malloc(19*sizeof(double));
#if SAY_HI
  printf("collide() -- Hi!\n");
#endif /* SAY_HI */

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  for( n=0; n<lattice->NumNodes; n++)
  {
    compute_single_feq(  lattice, n, subs, feq);
    ftemp    = lattice->pdf[subs][n].ftemp;
    is_solid = lattice->solids[subs][n].is_solid;

    if( !is_solid)
    {
        // C O L L I D E 

        // f = ftemp - (1/tau[subs])( ftemp - feq)
        for( a=0; a<Q; a++)
        {
#if 1
          ftemp[a] = ftemp[a]
               - ftemp[a] / lattice->param.tau[subs]
               +   feq[a] / lattice->param.tau[subs] ;
#else
          ftemp[a] = ftemp[a] - ( ( ftemp[a] )
                            - ( feq[a]   ) ) / lattice->param.tau[subs];
#endif
        } /* for( a=0; a<=8; a++) */

    } /* if( !( is_solid & BC_SOLID_NODE)) */

    else // is_solid & BC_SOLID_NODE
    {
      // B O U N C E B A C K 
        if( subs==0)
        {
          // Usual non-slip bounce-back condition.
          f = ftemp[ E];
      	  ftemp[ E] = ftemp[ W];
          ftemp[ W] = f;
          f = ftemp[ N];
      	  ftemp[ N] = ftemp[ S];
          ftemp[ S] = f;
          f = ftemp[ T];
      	  ftemp[ T] = ftemp[ B];
          ftemp[ B] = f;
          f = ftemp[NW];
      	  ftemp[NW] = ftemp[SE];
          ftemp[SE] = f;
          f = ftemp[NE];
      	  ftemp[NE] = ftemp[SW];
          ftemp[SW] = f;
          f = ftemp[TW];
      	  ftemp[TW] = ftemp[BE];
          ftemp[BE] = f;
          f = ftemp[TE];
      	  ftemp[TE] = ftemp[BW];
          ftemp[BW] = f;
          f = ftemp[TN];
      	  ftemp[TN] = ftemp[BS];
          ftemp[BS] = f;
          f = ftemp[TS];
      	  ftemp[TS] = ftemp[BN];
          ftemp[BN] = f;

        } /* if( subs==0) */

#if NUM_FLUID_COMPONENTS==2
        else // subs==1
        {
		
           // Usual non-slip bounce-back condition.

          f = ftemp[ E];
      	  ftemp[ E] = ftemp[ W];
          ftemp[ W] = f;
          f = ftemp[ N];
      	  ftemp[ N] = ftemp[ S];
          ftemp[ S] = f;
          f = ftemp[ T];
      	  ftemp[ T] = ftemp[ B];
          ftemp[ B] = f;
          f = ftemp[NW];
      	  ftemp[NW] = ftemp[SE];
          ftemp[SE] = f;
          f = ftemp[NE];
      	  ftemp[NE] = ftemp[SW];
          ftemp[SW] = f;
          f = ftemp[TW];
      	  ftemp[TW] = ftemp[BE];
          ftemp[BE] = f;
          f = ftemp[TE];
      	  ftemp[TE] = ftemp[BW];
          ftemp[BW] = f;
          f = ftemp[TN];
      	  ftemp[TN] = ftemp[BS];
          ftemp[BS] = f;
          f = ftemp[TS];
      	  ftemp[TS] = ftemp[BN];
          ftemp[BN] = f;

        } /* if( subs==0) else*/
#endif /* NUM_FLUID_COMPONENTS==2 */

    } /* if( !( is_solid & BC_SOLID_NODE)) else */


  } /* for( n=0; n<lattice_NumNodes; n++) */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if SAY_HI
  printf("collide() -- Bye!\n");
#endif /* SAY_HI */

} /* void collide( lattice_ptr lattice) */


#else
void collide( lattice_ptr lattice)
{
  double *feq;
  double *f;
  double *ftemp;

  double omega;

  int    is_solid;

  int    n, a;

  int    subs;

#if SAY_HI
  printf("collide() -- Hi!\n");
#endif /* SAY_HI */

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
  for( n=0; n<lattice->NumNodes; n++)
  {
    feq      = lattice->pdf[subs][n].feq;
    f        = lattice->pdf[subs][n].f;
    ftemp    = lattice->pdf[subs][n].ftemp;
    is_solid = lattice->solids[subs][n].is_solid;

    if( !is_solid)
    {
        // C O L L I D E 

        // f = ftemp - (1/tau[subs])( ftemp - feq)
        for( a=0; a<Q; a++)
        {
#if 1
          f[a] = ftemp[a]
               - ftemp[a] / lattice->param.tau[subs]
               +   feq[a] / lattice->param.tau[subs] ;
#else
          f[a] = ftemp[a] - ( ( ftemp[a] )
                            - ( feq[a]   ) ) / lattice->param.tau[subs];
#endif
        } /* for( a=0; a<=8; a++) */

#if PUKE_NEGATIVE_DENSITIES
        for( a=0; a<Q; a++)
        {
          if( *f < 0.)
          {
            printf("\n");
            printf(
              "collide() -- Node %d (%d,%d), subs %d, "
              "has negative density %20.17f "
              "in direction %d "
              "at timestep %d. Exiting!\n", 
              n, n%lattice->param.LX, 
                 n/lattice->param.LX,
                 subs,
                 f[a], a,
                 lattice->time             );
            printf("\n");
            process_exit(1);
          }
        } /* for( a=0; a<Q; a++) */
#endif /* PUKE_NEGATIVE_DENSITIES */

    } /* if( !( is_solid & BC_SOLID_NODE)) */

    else // is_solid & BC_SOLID_NODE
    {
      // B O U N C E B A C K 

      if(   lattice->param.bc_slip_north
         && n >= lattice->NumNodes - lattice->param.LX)
      {
        // Slip condition on north boundary.

        f[1] = ftemp[1];
        f[2] = ftemp[4];
        f[3] = ftemp[3];
        f[4] = ftemp[2];
        f[5] = ftemp[8];
        f[6] = ftemp[7];
        f[7] = ftemp[6];
        f[8] = ftemp[5];

      } /* if(   lattice->param.bc_slip_north && ... ) */

      else
      {
        if( subs==0)
        {
          // Usual non-slip bounce-back condition.
          /*
          //   N -> S
          //   S -> N
          //   E -> W
          //   W -> E
          //   T -> B
          //   B -> T
          */
          f[ E] = ftemp[ W];
          f[ W] = ftemp[ E];
          f[ N] = ftemp[ S];
          f[ S] = ftemp[ N];
          f[ T] = ftemp[ B];
          f[ B] = ftemp[ T];
          f[NW] = ftemp[SE];   
          f[NE] = ftemp[SW];   
          f[SW] = ftemp[NE];   
          f[SE] = ftemp[NW];   
          f[TW] = ftemp[BE];   
          f[TE] = ftemp[BW];   
          f[BW] = ftemp[TE];   
          f[BE] = ftemp[TW];   
          f[TN] = ftemp[BS];   
          f[TS] = ftemp[BN];   
          f[BN] = ftemp[TS];   
          f[BS] = ftemp[TN];   

        } /* if( subs==0) */

#if NUM_FLUID_COMPONENTS==2
        else // subs==1
        {
		
#if INAMURO_SIGMA_COMPONENT
          if( lattice->param.bc_sigma_slip)
          {
            //
            // Slip BC for solute on side walls.
            // Will this make a difference on Taylor dispersion?
            //
            if( lattice->FlowDir == /*Vertical*/2)
            {
              if(   /*west*/(n  )%lattice->param.LX == 0
                 || /*east*/(n+1)%lattice->param.LX == 0)
              {
                // Slip condition on east/west boundary.

                f[1] = ftemp[3];
                f[2] = ftemp[2];
                f[3] = ftemp[1];
                f[4] = ftemp[4];
                f[5] = ftemp[6];
                f[6] = ftemp[5];
                f[7] = ftemp[8];
                f[8] = ftemp[7];

              }
            }
            else if( lattice->FlowDir == /*Horizontal*/1)
            {
              if(   /*north*/ n >= lattice->NumNodes - lattice->param.LX
                 || /*south*/ n <  lattice->param.LX )
              {
                // Slip condition on north/south boundary.

                f[1] = ftemp[1];
                f[2] = ftemp[4];
                f[3] = ftemp[3];
                f[4] = ftemp[2];
                f[5] = ftemp[8];
                f[6] = ftemp[7];
                f[7] = ftemp[6];
                f[8] = ftemp[5];

              }
              else
              {
                // ERROR: Solid exists somewhere other than as side walls.
                printf("%s (%d) >> "
                  "ERROR: "
                  "bc_sigma_slip is on. "
                  "FlowDir is determined to be horizontal. "
                  "Encountered solid node somewhere other than side walls. "
                  "That situation is not supported. "
                  "Exiting!", __FILE__, __LINE__);
                process_exit(1);
              }
            }
            else
            {
              printf("%s (%d) >> "
                "FlowDir is indeterminate. "
                "Cannot apply slip BC (bc_sigma_slip). "
                "Exiting!", __FILE__, __LINE__);
              process_exit(1);
            }

          } /* if( lattice->param.bc_sigma_slip) */

          else
          {
#endif /* INAMURO_SIGMA_COMPONENT */
            // Usual non-slip bounce-back condition.

            f[ E] = ftemp[ W];
            f[ W] = ftemp[ E];
            f[ N] = ftemp[ S];
            f[ S] = ftemp[ N];
            f[ T] = ftemp[ B];
            f[ B] = ftemp[ T];
            f[NW] = ftemp[SE];   
            f[NE] = ftemp[SW];   
            f[SW] = ftemp[NE];   
            f[SE] = ftemp[NW];   
            f[TW] = ftemp[BE];   
            f[TE] = ftemp[BW];   
            f[BW] = ftemp[TE];   
            f[BE] = ftemp[TW];   
            f[TN] = ftemp[BS];   
            f[TS] = ftemp[BN];   
            f[BN] = ftemp[TS];   
            f[BS] = ftemp[TN];   

#if INAMURO_SIGMA_COMPONENT
          } /* if( lattice->param.bc_sigma_slip) else */
#endif /* INAMURO_SIGMA_COMPONENT */

        } /* if( subs==0) else*/
#endif /* NUM_FLUID_COMPONENTS==2 */

      } /* if(   lattice->param.bc_slip_north && ... ) else */

    } /* if( !( is_solid & BC_SOLID_NODE)) else */


  } /* for( n=0; n<lattice_NumNodes; n++) */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if SAY_HI
  printf("collide() -- Bye!\n");
#endif /* SAY_HI */

} /* void collide( lattice_ptr lattice) */


#endif
///ccccccccccccccccccccccccccccccccccccccccccccccc
//ccccccccccccccccccccccccccccccccccccccccccccccccc




