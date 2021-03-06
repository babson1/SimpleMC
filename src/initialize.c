#include "header.h"

Parameters *init_parameters(void)
{
  Parameters *p = malloc(sizeof(Parameters));

  p->n_particles = 1000000;
  p->n_batches = 10;
  p->n_generations = 1;
  p->n_active = 10;
  p->bc = REFLECT;
  p->n_nuclides = 1;
  p->tally = TRUE;
  p->n_bins = 16;
  p->seed = 1;
  p->nu = 2.5;
  p->xs_f = 0.012;
  p->xs_a = 0.03;
  p->xs_s = 0.27;
  p->gx = 400;
  p->gy = 400;
  p->gz = 400;
  p->write_tally = FALSE;
  p->write_keff = FALSE;
  p->tally_file = NULL;
  p->keff_file = NULL;
  
	//Quickly Added Stuff 6-23-2020
  MPI_Datatype dots, base[2]={MPI_INT, MPI_DOUBLE};
  int blocks[2]={2,8};
  MPI_Aint offsets[2], lb, extent;
  
  MPI_Type_get_extent(MPI_INT, &lb, &extent);
  offsets[0]=lb; offsets[1]=blocks[0]*extent+lb;
  MPI_Type_create_struct(2, blocks, offsets, base, &dots); 
  MPI_Type_commit(&dots);
  p->type=dots; 
  return p;
}

Geometry *init_geometry(Parameters *parameters)
{
  Geometry *g = malloc(sizeof(Geometry));

	int src, dest; 

  g->x = parameters->gx;
  g->y = parameters->gy;
  g->z = parameters->gz;
  g->bc = parameters->bc;
	
	for(int i=0; i<3; i++){
	MPI_Cart_shift(parameters->comm, i, 1, &src, &dest);
	g->nay[i*2]=src; g->nay[(i*2)+1]=dest;
}
  return g;
}

Tally *init_tally(Parameters *parameters)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->nx = parameters->n_bins/parameters->pX;
	t->ny = parameters->n_bins/parameters->pY;
	t->nz = parameters->n_bins/parameters->pZ;
  
	t->dx = parameters->gx/parameters->n_bins;
  t->dy = parameters->gy/parameters->n_bins;
  t->dz = parameters->gz/parameters->n_bins;
  t->flux = calloc(t->nx*t->ny*t->nz, sizeof(double));

  return t;
}

Material *init_material(Parameters *parameters)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, scattering, total, atomic density)
  Nuclide macro = {parameters->xs_f, parameters->xs_a, parameters->xs_s,
     parameters->xs_f + parameters->xs_a + parameters->xs_s, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = parameters->n_nuclides;
  m->nuclides = malloc(m->n_nuclides*sizeof(Nuclide));

  // Generate some arbitrary microscopic cross section values and atomic
  // densities for each nuclide in the material such that the total macroscopic
  // cross sections evaluate to what is hardwired above
  for(i=0; i<m->n_nuclides; i++){
    if(i<m->n_nuclides-1){
      m->nuclides[i].atom_density = rn()*macro.atom_density;
      macro.atom_density -= m->nuclides[i].atom_density;
    }
    else{
      m->nuclides[i].atom_density = macro.atom_density;
    }
    m->nuclides[i].xs_a = rn();
    sum.xs_a += m->nuclides[i].xs_a * m->nuclides[i].atom_density;
    m->nuclides[i].xs_f = rn();
    sum.xs_f += m->nuclides[i].xs_f * m->nuclides[i].atom_density;
    m->nuclides[i].xs_s = rn();
    sum.xs_s += m->nuclides[i].xs_s * m->nuclides[i].atom_density;
  }
  for(i=0; i<m->n_nuclides; i++){
    m->nuclides[i].xs_a /= sum.xs_a/macro.xs_a;
    m->nuclides[i].xs_f /= sum.xs_f/macro.xs_f;
    m->nuclides[i].xs_s /= sum.xs_s/macro.xs_s;
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_s;
  }

  m->xs_f = parameters->xs_f;
  m->xs_a = parameters->xs_a;
  m->xs_s = parameters->xs_s;
  m->xs_t = parameters->xs_a + parameters->xs_s;

  return m;
}

Bank *init_source_bank(Parameters *parameters, Geometry *geometry)
{
  unsigned long i_p; // index over particles
  Bank *source_bank;

  // Initialize source bank
	if(parameters->rank == 0)
  	source_bank = init_bank(parameters->n_particles);
	else
		source_bank = init_bank(parameters->n_particles/parameters->size);

  // Sample source particles
  for(i_p=0; i_p<parameters->n_particles; i_p++){
    sample_source_particle(geometry, &(source_bank->p[i_p]));
    source_bank->n++;
  }

  return source_bank;
}

Bank *init_fission_bank(Parameters *parameters)
{
  Bank *fission_bank;
	if(parameters->rank == 0) 
		fission_bank = init_bank(2*parameters->n_particles);
	else
		fission_bank = init_bank(2*parameters->n_particles/parameters->size);

  return fission_bank;
}
// old todo  
Bank *init_bank(unsigned long n_particles)
{
  Bank *b = malloc(sizeof(Bank));
  b->p = malloc(n_particles*sizeof(Particle));
  b->sz = n_particles;
  b->n = 0;
  b->resize = resize_particles;

  return b;
}

void sample_source_particle(Geometry *geometry, Particle *p)
{
  p->alive = TRUE;
	p->hit = FALSE; 
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*geometry->x;
  p->y = rn()*geometry->y;
  p->z = rn()*geometry->z;
	p->coord[0] = p->x/geometry->xl;
	p->coord[1] = p->y/geometry->yl;
	p->coord[2] = p->z/geometry->zl;
	p->lx = p->x - (geometry->xl*p->coord[0]);
	p->ly = p->y - (geometry->yl*p->coord[1]);
	p->lz = p->z - (geometry->zl*p->coord[2]);
 
 return;
}

void sample_fission_particle(Particle *p, Particle *p_old)
{
  p->alive = TRUE;
  p->hit = FALSE;
	p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = p_old->x;
  p->y = p_old->y;
  p->z = p_old->z;

  return;
}

void resize_particles(Bank *b)
{
  b->p = realloc(b->p, sizeof(Particle)*2*b->sz);
  b->sz = 2*b->sz;

  return;
}

void free_bank(Bank *b)
{
  free(b->p);
  b->p = NULL;
  free(b);
  b = NULL;

  return;
}

void free_material(Material *material)
{
  free(material->nuclides);
  material->nuclides = NULL;
  free(material);
  material = NULL;

  return;
}

void free_tally(Tally *tally)
{
  free(tally->flux);
  tally->flux = NULL;
  free(tally);
  tally = NULL;

  return;
}
