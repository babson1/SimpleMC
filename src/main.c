#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *parameters; // user defined parameters
  Geometry *geometry; // homogenous cube geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *tally; // scalar flux tally
  double *keff; // effective multiplication factor
  double t1, t2; // timers
  MPI_Comm cart;
  int rank, size; 
  int dims[3], wrap[3]={TRUE, TRUE, TRUE};

  MPI_Init(&argc, &argv);

  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
  parse_parameters(parameters);
  read_CLI(argc, argv, parameters);

	parameters->comm=cart;
  dims[0]=parameters->pX ; dims[1]=parameters->pY; dims[2]=parameters->pZ;
  MPI_Cart_create(parameters->comm, 3, dims, wrap, TRUE, &cart);
  MPI_Comm_rank(cart, &rank);
  MPI_Comm_size(cart, &size);
 
  parameters->comm = cart;
  parameters->comm = rank;
  parameters->comm = size;

if(size != dims[0]*dims[1]*dims[2]) {
  if(parameters->rank == 0)
    printf("Product of procs per dimension != total proc count\n");
  MPI_Barrier(cart);
  MPI_Abort(cart, 1);
}  

 MPI_Cart_coords(parameters->comm, rank, 3, parameters->coords);
  if(rank==0){ 
   print_parameters(parameters);
 }

 // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Create files for writing results to
  if(rank==0) 
     init_output(parameters);



  // Set up geometry
  geometry = init_geometry(parameters);

  // Set up material
  material = init_material(parameters);

  // Set up tallies
  tally = init_tally(parameters);

  // Create source bank and initial source distribution
  source_bank = init_source_bank(parameters, geometry);

  // Create fission bank
  fission_bank = init_fission_bank(parameters);

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));
 if(rank==0){
  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");
}

  // Start time

 MPI_Barrier(cart);
 t1 = MPI_Wtime();

  run_eigenvalue(parameters, geometry, material, source_bank, fission_bank, tally, keff);

  // Stop time
 MPI_Barrier(cart);
 t2 = MPI_Wtime();
 

 if(rank==0)
  printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free(keff);
  free_tally(tally);
  free_bank(fission_bank);
  free_bank(source_bank);
  free_material(material);
  free(geometry);
  free(parameters);



  MPI_Finalize();
  return 0;}
