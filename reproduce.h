#ifndef __REPRODUCE_H
#define __REPRODUCE_H


/// create chromosomes for a parent using haldane's map function 

/// NOTES:
// what about unequal crossing over?

void transmit_chromosome ( individual &parent, vector<gene*> &new_chromosome ) {

    /// go through by position 
    int it_mom = 0 ; 
    int it_dad = 0 ; 

    // if mom == 1, we are currently on mom's chromosome
    float position = 0 ; 
    bool mom = gsl_ran_bernoulli( rng, 0.5 ) ; 

    /// while the positiion is lt both 
    while ( it_mom < parent.maternal_trnas.size() && it_dad < parent.paternal_trnas.size() ) {

        /// both have it at the same site 
        if ( (*parent.maternal_trnas[it_mom]).locus == (*parent.paternal_trnas[it_dad]).locus ) {

            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
            if ( mom ) { 
                new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ; 
            }
            else {
                new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
            }
            position = (*parent.maternal_trnas[it_mom]).locus ; 
            it_mom ++ ;
            it_dad ++ ; 
        } 
        /// if it's dad that's next
        else if ( (*parent.maternal_trnas[it_mom]).locus > (*parent.paternal_trnas[it_dad]).locus ) { 
            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.paternal_trnas[it_dad]).locus - position ) ) / 2 ) ) mom = !mom ; 
            if ( !mom ) { 
                new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
            }
            position = (*parent.paternal_trnas[it_dad]).locus ;
            it_dad ++ ; 
        }
        /// otherwise, mom is next 
        else  { 
            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
            if ( mom ) { 
                new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ; 
            }
            position = (*parent.maternal_trnas[it_mom]).locus ;
            it_mom ++ ; 
        }
    }

    /// if it's dad that's last
    while ( it_dad < parent.paternal_trnas.size() ) { 
        if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.paternal_trnas[it_dad]).locus - position ) ) / 2 ) ) mom = !mom ; 
        if ( !mom ) { 
            new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
        }
        position = (*parent.paternal_trnas[it_dad]).locus ;
        it_dad ++ ; 
    }
    /// otherwise, mom is next 
    while ( it_mom < parent.maternal_trnas.size() ) { 
        if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
        if ( mom ) { 
            new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ; 
        }
        position = (*parent.maternal_trnas[it_mom]).locus ;
        it_mom ++ ; 
    }

}

/// compute fitness
void reproduce( const double *fitness, vector<individual> &population, vector<individual> &new_population, cmd_line &options ) {

    // populate parent multinomial samplings
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc( population.size(), fitness ) ;

    /// iterate through all individuals and draw parents + recomb
    for ( int i = 0 ; i < population.size() ; i ++ ) {
    
        /// new individual
        individual new_ind ;
        
        /// grab mom and dad
        int mom = gsl_ran_discrete( rng, g ) ;
        int dad = gsl_ran_discrete( rng, g ) ;

        // NOW sort, and by locus, so this might fix the recombination issue:
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);
               
        // get their chromosomes 
        transmit_chromosome( population[mom], new_ind.maternal_trnas ) ; 
        transmit_chromosome( population[dad], new_ind.paternal_trnas ) ; 

        // swap individual in
        swap( new_population[i], new_ind ) ;
    }

    // de-allocate space for multinomial samplings
    if (g) {
        gsl_ran_discrete_free( g ) ;
    }
}

#endif
