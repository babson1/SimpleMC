#include "header.h"

// Simple flux tally
void score_tally(Parameters *parameters, Material *material, Tally *tally, Particle *p)
{
  int ix, iy, iz;
  double vol;

  // Volume
  vol = tally->dx * tally->dy * tally->dz;

  // Find the indices of the grid box of the particle
  ix = p->x/tally->dx;
  iy = p->y/tally->dy;
  iz = p->z/tally->dz;

  // Scalar flux
  tally->flux[ix + tally->nx*iy + tally->ny*tally->nz*iz] += 1./(vol * material->xs_t * parameters->n_particles);

  return;
}

void reset_tally(Tally *tally)
{
  memset(tally->flux, 0, tally->nx*tally->ny*tally->nz*sizeof(double));

  return;
}
