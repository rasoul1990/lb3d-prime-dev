//##############################################################################
//
// stream.c
//
#if SAVE_MEMO


#else
void stream( lattice_ptr lattice)
{
  double *f;

  int     subs;

  int     a,  n,
          ni, nj, nk,
          i,  ip, in,
          j,  jp, jn,
          k,  kp, kn;

#if SAY_HI
  printf("%s %d >> stream() -- Hi!\n", __FILE__,__LINE__);
#endif /* SAY_HI */

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
   ni = get_LX( lattice);
   nj = get_LY( lattice);
   nk = get_LZ( lattice);

#if 0
   // Contrived f data for debugging...
  if( lattice->time == 1)
  {
   display_warning_about_contrived_data( lattice);
   for( n=0; n<get_NumNodes(lattice); n++)
   {
     //for( a=0; a<19; a++)
     //{
     //  lattice->pdf[subs][ n].f[a] = get_g_StartNode(lattice)+n+1;
     //}

     lattice->pdf[subs][n].f[ C] = 0.;
     lattice->pdf[subs][n].f[ E] = 0.;
     lattice->pdf[subs][n].f[ W] = 0.;
     lattice->pdf[subs][n].f[ N] = 0.;
     lattice->pdf[subs][n].f[ S] = 0.;
     lattice->pdf[subs][n].f[NW] = 0.;
     lattice->pdf[subs][n].f[NE] = 0.;
     lattice->pdf[subs][n].f[SW] = 0.;
     lattice->pdf[subs][n].f[SE] = 0.;

     lattice->pdf[subs][n].f[ T] = (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[TW] = (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[TE] = (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[TN] = (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[TS] = (get_g_StartNode(lattice) + n + 1);

     lattice->pdf[subs][n].f[ B] = get_g_NumNodes(lattice)
                                   - (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[BW] = get_g_NumNodes(lattice)
                                   - (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[BE] = get_g_NumNodes(lattice)
                                   - (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[BN] = get_g_NumNodes(lattice)
                                   - (get_g_StartNode(lattice) + n + 1);
     lattice->pdf[subs][n].f[BS] = get_g_NumNodes(lattice)
                                   - (get_g_StartNode(lattice) + n + 1);

     for( a=0; a<19; a++)
     {
       lattice->pdf[subs][ n].f[a] += a;
     }

   }
  }
#endif

   //dump_north_pointing_pdfs( lattice, subs, -1, "f before streaming", 1);
   //dump_south_pointing_pdfs( lattice, subs, -1, "f before streaming", 1);

   process_send_recv_begin( lattice, subs);

   // S T R E A M
   //###########################################################################
   for( k=0; k<nk; k++)
   {
//!lattice->param.GZL
//     if(!lattice->param.GZL )
     {
     kp = ( k<nk-1)?( k+1):( 0   );
     kn = ( k>0   )?( k-1):( nk-1);
     }
//     else
//     {
//     kp = ( k<nk-1)?( k+1):( nk-1);

//     kn = ( k>0   )?( k-1):( 0   );
//     }

     for( j=0; j<nj; j++)
     {
//     if(!lattice->param.GZL)
            {
        jp = ( j<nj-1)?( j+1):( 0   );
                    jn = ( j>0   )?( j-1):( nj-1);
      }
//     else
//            {
//        jp = ( j<nj-1)?( j+1):( nj-1   );
//                    jn = ( j>0   )?( j-1):( 0);
//      }

       for( i=0; i<ni; i++)
       {
//        if(!lattice->param.GZL)
    {
       ip = ( i<ni-1)?( i+1):( 0   );
             in = ( i>0   )?( i-1):( ni-1);
    }
//  else
//    {
//       ip = ( i<ni-1)?( i+1):( ni-1   );
//             in = ( i>0   )?( i-1):( 0);
//    }


         f = lattice->pdf[subs][ XYZ2N( i , j , k , ni, nj)].f;

         lattice->pdf[subs][ XYZ2N( i , j , k , ni, nj)].ftemp[ C] = f[ C];

         lattice->pdf[subs][ XYZ2N( ip, j , k , ni, nj)].ftemp[ E] = f[ E];
         lattice->pdf[subs][ XYZ2N( in, j , k , ni, nj)].ftemp[ W] = f[ W];
         lattice->pdf[subs][ XYZ2N( i , jp, k , ni, nj)].ftemp[ N] = f[ N];
         lattice->pdf[subs][ XYZ2N( i , jn, k , ni, nj)].ftemp[ S] = f[ S];

         lattice->pdf[subs][ XYZ2N( in, jp, k , ni, nj)].ftemp[NW] = f[NW];
         lattice->pdf[subs][ XYZ2N( ip, jp, k , ni, nj)].ftemp[NE] = f[NE];
         lattice->pdf[subs][ XYZ2N( in, jn, k , ni, nj)].ftemp[SW] = f[SW];
         lattice->pdf[subs][ XYZ2N( ip, jn, k , ni, nj)].ftemp[SE] = f[SE];
//  if(index !=1)
  {
         lattice->pdf[subs][ XYZ2N( i , j , kp, ni, nj)].ftemp[ T] = f[ T];
   lattice->pdf[subs][ XYZ2N( in, j , kp, ni, nj)].ftemp[TW] = f[TW];
         lattice->pdf[subs][ XYZ2N( ip, j , kp, ni, nj)].ftemp[TE] = f[TE];
         lattice->pdf[subs][ XYZ2N( i , jp, kp, ni, nj)].ftemp[TN] = f[TN];
         lattice->pdf[subs][ XYZ2N( i , jn, kp, ni, nj)].ftemp[TS] = f[TS];
  }
//  if(index !=2)
  {
         lattice->pdf[subs][ XYZ2N( i , j , kn, ni, nj)].ftemp[ B] = f[ B];
   lattice->pdf[subs][ XYZ2N( in, j , kn, ni, nj)].ftemp[BW] = f[BW];
         lattice->pdf[subs][ XYZ2N( ip, j , kn, ni, nj)].ftemp[BE] = f[BE];
         lattice->pdf[subs][ XYZ2N( i , jp, kn, ni, nj)].ftemp[BN] = f[BN];
         lattice->pdf[subs][ XYZ2N( i , jn, kn, ni, nj)].ftemp[BS] = f[BS];
  }
       } /* if( i=0; i<ni; i++, n++) */
     } /* if( j=0; j<nj; j++) */
   } /* if( k=0; k<nk; k++) */

   process_send_recv_end( lattice, subs);

   //dump_north_pointing_pdfs( lattice, subs, -1, "ftemp after Streaming", 2);
   //dump_south_pointing_pdfs( lattice, subs, -1, "ftemp after Streaming", 2);

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if SAY_HI
  printf("%s %d >> stream() -- Bye!\n", __FILE__,__LINE__);
#endif /* SAY_HI */

} /* void stream( lattice_ptr lattice) */


#endif
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


