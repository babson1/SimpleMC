#include "header.h"

void run_eigenvalue(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, double *keff)
{
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p, i_s; // index over particles
  double keff_gen = 1; // keff of generation
  double keff_batch; // keff of batch
  double keff_mean; // keff mean over active batches
  double keff_std; // keff standard deviation over active batches


  // Loop over batches
  for(i_b=0; i_b<parameters->n_batches; i_b++){

    keff_batch = 0;

    // Turn on tallying and increment index in active batches
    if(i_b >= parameters->n_batches - parameters->n_active){
      i_a++;
      if(parameters->tally == TRUE){
        tally->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<parameters->n_generations; i_g++){
		

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);
			
			//TODO Disperse Function Call, MUST BE CREATED 
			disperse(parameters, geometry, source_bank);
			
			MPI_Scan(&source_bank->n, &i_s, 1, MPI_UNSIGNED_LONG, MPI_SUM, parameters->comm);
			i_s -=  source_bank->n; 

      // Loop over particles
			parameters->dead = 0;
			while(parameters->dead < parameters->n_particles){
      for(i_p=0; i_p<source_bank->n; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers from the initial
	// seed. This allows for reproducibility of the particle history.
        rn_skip((i_b*parameters->n_generations + i_g)*parameters->n_particles + i_s);
				i_s++;

        // Transport the next particle
        transport(parameters, geometry, material, source_bank, fission_bank, tally, &(source_bank->p[i_p]));
      }
	//TODO Create disperse function
			disperse(parameters, geometry, source_bank);
      // Find Collision Locations
			for(i_p=0; i_p<source_bank->n; i_p++){
				if(source_bank->p[i_p].hit == TRUE) {
					if(tally->tallies_on == TRUE)
						score_tally(parameters, material, tally, &source_bank->p[i_p]); 
					// livers live
						source_bank->p[i_p].hit = FALSE; 
				}
			}
			// let processors sum count of dead particles
			MPI_Allreduce(&source_bank->dead, &parameters->dead, 1, MPI_UNSIGNED_LONG, MPI_SUM, parameters->comm);
		}
			// Switch RNG stream off tracking
      set_stream(STREAM_OTHER);
      rn_skip(i_b*parameters->n_generations + i_g);
		
      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      i_p = synchronize_bank(source_bank, fission_bank,parameters);
     	 
			// Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) i_p / source_bank->n;
      keff_batch += keff_gen;
		}

    // Calculate k effective
    keff_batch /= parameters->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }
    calculate_keff(keff, &keff_mean, &keff_std, i_a+1);

    // Tallies for this realization
    if(tally->tallies_on == TRUE){
      if(parameters->write_tally == TRUE){
        write_tally(tally, parameters->tally_file);
      }
      reset_tally(tally);
    }

    // Status text
		if(parameters->rank == 0)
    print_status(i_a, i_b, keff_batch, keff_mean, keff_std);
  }

  // Write out keff
	if(parameters->rank == 0){
		if(parameters->write_keff == TRUE){
    write_keff(keff, parameters->n_active, parameters->keff_file);
  	}	
	}
  return;
}

unsigned long	synchronize_bank(Bank *source_bank, Bank *fission_bank, Parameters *parameters)
{
  unsigned long i, j;
  unsigned long n_s = source_bank->n;
  unsigned long n_f = fission_bank->n;
	unsigned long *n_ptr = malloc(parameters->size*sizeof(unsigned long));
	MPI_Status stat; 

	MPI_Gather(&fission_bank->n,1, MPI_UNSIGNED_LONG, n_ptr,1, MPI_UNSIGNED_LONG, 0, parameters->comm);  
		if(parameters->rank == 0){
				for(i=1; i<parameters->size; i++)
					n_f += n_ptr[i]; 
 				if(fission_bank->sz<n_f)
					fission_bank->resize(fission_bank);
				for(i=1; i<parameters->size; i++){
					MPI_Recv(&fission_bank->p[fission_bank->n],n_ptr[i], parameters->type, i, 0, parameters->comm, &stat);
					fission_bank->n += n_ptr[i]; 
				}
				n_f = fission_bank->n;
		}
		else
			MPI_Send(fission_bank->p, fission_bank->n, parameters->type, 0, 0, parameters->comm);
	
		if(parameters->rank == 0){
	 // If the fission bank is larger than the source bank, randomly select
  // n_particles sites from the fission bank to create the new source bank
				if(n_f >=n_s){
				
				// Copy first n_particles sites from fission bank to source bank
    		memcpy(source_bank->p, fission_bank->p, n_s*sizeof(Particle));

    // Replace elements with decreasing probability, such that after final
    // iteration each particle in fission bank will have equal probability of
    // being selected for source bank
    for(i=n_s; i<n_f; i++){
      j = rni(0, i+1);
      if(j<n_s){
        memcpy(&(source_bank->p[j]), &(fission_bank->p[i]), sizeof(Particle));
      }
    }
  }

  // If the fission bank is smaller than the source bank, use all fission bank
  // sites for the source bank and randomly sample remaining particles from
  // fission bank
  else{

    // First randomly sample particles from fission bank
    for(i=0; i<(n_s-n_f); i++){
      j = rni(0, n_f);
      memcpy(&(source_bank->p[i]), &(fission_bank->p[j]), sizeof(Particle));
    }

    // Fill remaining source bank sites with fission bank
    memcpy(&(source_bank->p[n_s-n_f]), fission_bank->p, n_f*sizeof(Particle));
  	}	
 }
  MPI_Barrier(parameters->comm);
	//distribute to source banks
	if(parameters->rank == 0)
		MPI_Scatter(source_bank->p, n_s/parameters->size, parameters->type,MPI_IN_PLACE, n_s/parameters->size, parameters->type, 0, parameters->comm); 
	else
		MPI_Scatter(NULL, n_s/parameters->size,parameters->type, MPI_IN_PLACE, n_s/parameters->size, parameters->type, 0, parameters->comm);

	fission_bank->n = 0;
	source_bank->n = n_s/parameters->size;
	source_bank->dead=0;

  return n_f;
}

void calculate_keff(double *keff, double *mean, double *std, int n)
{
  int i;

  *mean = 0;
  *std = 0;

  // Calculate mean
  for(i=0; i<n; i++){
    *mean += keff[i];
  }
  *mean /= n;

  // Calculate standard deviation
  for(i=0; i<n; i++){
    *std += pow(keff[i] - *mean, 2);
  }
  *std = sqrt(*std/(n-1));

  return;
}
