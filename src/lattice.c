//##############################################################################
//
// lattice.c
//

int get_LX( lattice_ptr lattice) { return lattice->param.LX;}
int get_LY( lattice_ptr lattice) { return lattice->param.LY;}
int get_LZ( lattice_ptr lattice) { return lattice->param.LZ;}
void set_LX( lattice_ptr lattice, const int arg_LX)
{
  lattice->param.LX = arg_LX;
}
void set_LY( lattice_ptr lattice, const int arg_LY)
{
  lattice->param.LY = arg_LY;
}
void set_LZ( lattice_ptr lattice, const int arg_LZ)
{
  lattice->param.LZ = arg_LZ;
}

int get_NumNodes( lattice_ptr lattice) { return lattice->NumNodes;}
void set_NumNodes( lattice_ptr lattice)
{
  lattice->NumNodes = get_LX( lattice)*get_LY( lattice)*get_LZ( lattice);
}

int is_solid( lattice_ptr lattice, const int n)
{
  if( lattice->solids[0][n].is_solid != 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int is_not_solid( lattice_ptr lattice, const int n)
{
  return !is_solid( lattice, n);
}

double get_rho( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->macro_vars[subs][n].rho;
}

double get_ux( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->macro_vars[subs][n].u[0];
}
double get_uy( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->macro_vars[subs][n].u[1];
}
double get_uz( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->macro_vars[subs][n].u[2];
}

