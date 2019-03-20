#ifndef __FITNESS_H
#define __FITNESS_H

/// compute fitness
void compute_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, double opt_fit, cmd_line &options ) {

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



    for ( int i = 0 ; i < population.size() ; i ++ ) {


    	double x_ind = 0.0 ;

    	// maternal x calculation block:
    	for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
    		if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
    			if ( (*population[i].maternal_trnas[g]).muts < 10 ) {
        			int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
        			x_ind += ( ((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
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
    			if ( (*population[i].paternal_trnas[g]).muts < 10 ) {
        			int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
        			x_ind += ( ((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
        		}
        		// else would mean 10 muts so sequence value is 0, so add 0, so do nothing
        	}
        	else {
        		x_ind += ( ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
        	}
        }

        fitness[i] = (1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) )) * ( exp(-1 * ( (pow((x_ind - options.fitness_mean), 2)) / (2 * (pow(options.fitness_sd, 2))) ) ) ) / opt_fit ;

    }
}


#endif