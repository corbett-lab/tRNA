#ifndef __SIMPLER_REPRODUCE_H
#define __SIMPLER_REPRODUCE_H

/// compute fitness
void reproduce( const double *fitness, vector<individual> &population, vector<individual> &new_population, cmd_line &options ) {
    
    // populate parent multinomial samplings
    const gsl_ran_discrete_t *g = gsl_ran_discrete_preproc( population.size(), fitness ) ;
    
    /// iterate through all individuals and draw parents + recomb
    for ( int i = 0 ; i < population.size() ; i ++ ) {
    
        /// new individual
        individual new_ind ;
        
        /// grab mom and dad
        int mom = gsl_ran_discrete( rng, g ) ;
        int dad = gsl_ran_discrete( rng, g ) ;
               
        // get their chromosomes
        if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
            new_ind.maternal_trnas = (population[mom]).maternal_trnas ;
        }
        else {
        	new_ind.maternal_trnas = (population[mom]).paternal_trnas ;
        }
        if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
            new_ind.paternal_trnas = (population[dad]).maternal_trnas ;
        }
        else {
        	new_ind.paternal_trnas = (population[dad]).paternal_trnas ;
        }

        // swap individual in
        swap( new_population[i], new_ind ) ;
    }
}

#endif