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
            f[ C] = ftemp[ C];
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

#if POROUS_MEDIA
    if( INAMURO_SIGMA_COMPONENT!=0 || subs==0)
    {
      int LX = get_LX( lattice);
      int LY = get_LY( lattice);
      int LZ = get_LZ( lattice);

      double ns;
      double* nsterm = lattice->nsterm;

      // Compute the solid density term for fluid component.
      for( n=0; n<lattice->NumNodes; n++)
      {
        if( is_not_solid(lattice, n))
        {
          // ns_flag = 0 ==> uniform ns value read from scalar
          // ns_flag in {1,2} ==> variable ns_value read from array
          if( lattice->param.ns_flag == 0) { ns = lattice->param.ns; }
          else { ns = lattice->ns[n].ns; }
/*  E */nsterm[Q*n+ E] = ns*( lattice->pdf[subs][n].ftemp[ W]
                            - lattice->pdf[subs][n].f    [ E]);
/*  W */nsterm[Q*n+ W] = ns*( lattice->pdf[subs][n].ftemp[ E]
                            - lattice->pdf[subs][n].f    [ W]);
/*  N */nsterm[Q*n+ N] = ns*( lattice->pdf[subs][n].ftemp[ S]
                            - lattice->pdf[subs][n].f    [ N]);
/*  S */nsterm[Q*n+ S] = ns*( lattice->pdf[subs][n].ftemp[ N]
                            - lattice->pdf[subs][n].f    [ S]);
/*  T */nsterm[Q*n+ T] = ns*( lattice->pdf[subs][n].ftemp[ B]
                            - lattice->pdf[subs][n].f    [ T]);
/*  B */nsterm[Q*n+ B] = ns*( lattice->pdf[subs][n].ftemp[ T]
                            - lattice->pdf[subs][n].f    [ B]);
/* NW */nsterm[Q*n+NW] = ns*( lattice->pdf[subs][n].ftemp[NE]
                            - lattice->pdf[subs][n].f    [NW]);
/* NE */nsterm[Q*n+NE] = ns*( lattice->pdf[subs][n].ftemp[NW]
                            - lattice->pdf[subs][n].f    [NE]);
/* SW */nsterm[Q*n+SW] = ns*( lattice->pdf[subs][n].ftemp[SE]
                            - lattice->pdf[subs][n].f    [SW]);
/* SE */nsterm[Q*n+SE] = ns*( lattice->pdf[subs][n].ftemp[SW]
                            - lattice->pdf[subs][n].f    [SE]);
/* TW */nsterm[Q*n+TW] = ns*( lattice->pdf[subs][n].ftemp[TE]
                            - lattice->pdf[subs][n].f    [TW]);
/* TE */nsterm[Q*n+TE] = ns*( lattice->pdf[subs][n].ftemp[TW]
                            - lattice->pdf[subs][n].f    [TE]);
/* BW */nsterm[Q*n+BW] = ns*( lattice->pdf[subs][n].ftemp[BE]
                            - lattice->pdf[subs][n].f    [BW]);
/* BE */nsterm[Q*n+BE] = ns*( lattice->pdf[subs][n].ftemp[BW]
                            - lattice->pdf[subs][n].f    [BE]);
/* TN */nsterm[Q*n+TN] = ns*( lattice->pdf[subs][n].ftemp[TS]
                            - lattice->pdf[subs][n].f    [TN]);
/* TS */nsterm[Q*n+TS] = ns*( lattice->pdf[subs][n].ftemp[TN]
                            - lattice->pdf[subs][n].f    [TS]);
/* BN */nsterm[Q*n+BN] = ns*( lattice->pdf[subs][n].ftemp[BS]
                            - lattice->pdf[subs][n].f    [BN]);
/* BS */nsterm[Q*n+BS] = ns*( lattice->pdf[subs][n].ftemp[BN]
                            - lattice->pdf[subs][n].f    [BS]);
        }
      } /* for( n=0; n<lattice_NumNodes; n++) */

      for( n=0; n<lattice->NumNodes; n++)
      {
        f = lattice->pdf[subs][n].f;

        if( is_not_solid(lattice, n))
        {
          for( a=1; a<Q; a++)
          {
#if 0
            // Store f in ftemp, because the compute_macro vars before
            // output_frames needs to have the pre-ns version of f.
            if( is_last_step_of_frame(lattice))
            {
              // This temp copy of f in ftemp prior to application
              // of ns only needs to be done before outputting a frame,
              // not after each timestep.
              lattice->pdf[subs][n].ftemp[a] = f[a];
            }
#endif
        // OLD: fout = fc + ns*( fc_(x+cdt) - fc(x))
        //        f +=      ns*( f_(x+cdt)  - f(x))
        //
        //  - - - >o< - - -   <------o------>   <--<---o--->-->
        //        1 2         3             4   3  5       6  4
        //                    ------>o<------
        //                          5 6
        //
        // CUR: fout = fc + ns*( fc_(x)     - fc(x)) = (1-ns)*fc + ns*fc_(x)
        //        f +=      ns*( f_(x)      - f(x))
        //
        //  - - - >o< - - -   <------o------>   <--<---o--->-->
        //        1 2         3             4   3  4       3  4
        //
        // NEW: fout = fc + ns*( fin_(x)    - fc(x)) = (1-ns)*fc + ns*fin_(x)
        //        f +=      ns*( ftemp_(x)  - f(x))
        //
        //  - - - >o< - - -   <------o------>   <--<- -o- ->-->
        //        1 2         3             4   3  1       2  4
          //
            // c.f., Walsh, Stuart D. C. and Burwinkle, Holly and Saar, Martin
            // O., A new partial-bounceback lattice-Boltzmann method for fluid
            // flow through heterogeneous media, COMPUTERS & GEOSCIENCES, 2009,
            // 35, 6, 1186-1193, JUN, ISI:000266544700013
            f[a] += nsterm[Q*n+a];

          } /* for( a=1; a<Q; a++) */
        } /* if( !( bc_type & BC_SOLID_NODE)) */
      } /* for( n=0; n<lattice->NumNodes; n++, f+=18) */

    } /* if( INAMURO_SIGMA_COMPONENT!=0 || subs==0) */
#endif

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if SAY_HI
  printf("collide() -- Bye!\n");
#endif /* SAY_HI */

} /* void collide( lattice_ptr lattice) */


#endif
///ccccccccccccccccccccccccccccccccccccccccccccccc
//ccccccccccccccccccccccccccccccccccccccccccccccccc




