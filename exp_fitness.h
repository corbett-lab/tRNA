#ifndef __EXP_FITNESS_H
#define __EXP_FITNESS_H

/// compute fitness
void compute_exp_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

	//////////////////////////////////
    //////////////////////////////////
    ///////// REGULAR METHOD ///////// 
    //////////////////////////////////
    //////////////////////////////////
    //
    // tRNA MODEL BASED ON OBSERVED DATA:
    // fitness = fitness_seq * fitness_exp
    // fitness_seq = 1 - e^((options.lambda_seq) * x)
    // fitness_exp = 1 - e^(-10.5 * x)

    // ^ no don't do this
    // gaussian with set mean and sd
    // float x = sum_all_tRNAs(tRNA.sequence * tRNA.expression)
    // fitness = ( 1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) ) ) * exp(-1 * ( (pow((x - fitness.mean), 2)) / (2 * pow(options.fitness_sd, 2)) )
    // ^ normalize this so that highest is 1.0 and all others are a fraction of that


    double temp_fitness [options.n] ;
    double max_fitness = 0.0 ;

    for ( int i = 0 ; i < population.size() ; ++i ) {


    	double x_ind = 0.0 ;

    	// maternal x calculation block:
    	for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
    		if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
    			if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                    if ( options.mutation_pathways == false ){
            			int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
            			x_ind += ( ((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
            		}
                    else{
                        int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                        x_ind += ( ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
                    }
                }
        		// else would mean 10 muts so sequence value is 0, so add 0, so do nothing
        	}
        	else {
        		x_ind += ( ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression) ) ;
        	}
        }

        // paternal x calculation block:
    	for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
    		if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
    			if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                    if ( options.mutation_pathways == false ){
            			int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
            			x_ind += ( ((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
            		}
                    else{
                        int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                        x_ind += ( ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
                    }
                }
        		// else would mean 10 muts so sequence value is 0, so add 0, so do nothing
        	}
        	else {
        		x_ind += ( ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
        	}
        }

        temp_fitness[i] = x_ind ;
        if ( x_ind > max_fitness ){
            max_fitness = x_ind ;
            //cout << max_fitness << endl ;
        }
    }

    for ( int i = 0 ; i < population.size() ; ++i ) {
        fitness[i] = 1.0 - exp(options.fitness_lambda * (temp_fitness[i] / max_fitness)) ;
        if ( fitness[i] > 1.0 ) {
            cout << "FITNESS > 1.0:\t" << fitness[i] << "\tFITNESS NOT ALLOWED TO BE > 1.0. EXITING PROGRAM." << endl ;
            exit(0);
        }
    }


}


#endif