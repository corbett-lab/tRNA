#ifndef __REPRODUCE_H
#define __REPRODUCE_H


/// create chromosomes for a parent using haldane's map function 

/// NOTES:
// what about unequal crossing over?

void transmit_chromosome ( Individual &parent, vector<Gene*> &new_chromosome ) {

    /// go through by position 
    int it_mom = 0 ; 
    int it_dad = 0 ; 

    // if mom == 1, we are currently on mom's chromosome
    float position = 0 ; 
    bool mom = gsl_ran_bernoulli( rng, 0.5 ) ; 

    /// while the positiion is lt both 
    while ( it_mom < parent.getPaternal_trnas().size() && it_dad < parent.getPaternal_trnas().size() ) {

        /// both have it at the same site 
        if ( (*parent.getMaternal_trnas()[it_mom]).getLocus() == (*parent.getPaternal_trnas()[it_dad]).getLocus() ) {

            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.getMaternal_trnas()[it_mom]).getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( mom ) { 
                new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
            }
            else {
                new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
            }
            position = (*parent.getMaternal_trnas()[it_mom]).getLocus() ;
            it_mom ++ ;
            it_dad ++ ; 
        } 
        /// if it's dad that's next
        else if ( (*parent.getMaternal_trnas()[it_mom]).getLocus() > (*parent.getPaternal_trnas()[it_dad]).getLocus() ) {
            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.getPaternal_trnas()[it_dad]).getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( !mom ) { 
                new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
            }
            position = (*parent.getPaternal_trnas()[it_dad]).getLocus() ;
            it_dad ++ ; 
        }
        /// otherwise, mom is next 
        else  { 
            if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.getMaternal_trnas()[it_mom]).getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( mom ) { 
                new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
            }
            position = (*parent.getMaternal_trnas()[it_mom]).getLocus() ;

            it_mom ++ ; 
        }
    }

    /// if it's dad that's last
    while ( it_dad < parent.getPaternal_trnas().size() ) {
        if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.getPaternal_trnas()[it_dad]).getLocus() - position ) ) / 2 ) ) mom = !mom ;
        if ( !mom ) { 
            new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
        }
        position = (*parent.getPaternal_trnas()[it_dad]).getLocus() ;
        it_dad ++ ; 
    }
    /// otherwise, mom is next 
    while ( it_mom < parent.getMaternal_trnas().size() ) {
        if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.getMaternal_trnas()[it_mom]).getLocus() - position ) ) / 2 ) ) mom = !mom ;
        if ( mom ) { 
            new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
        }
        position = (*parent.getMaternal_trnas()[it_mom]).getLocus() ;
        it_mom ++ ; 
    }

}

/// compute fitness
//void reproduce( const double *fitness, vector<individual*> &population, vector<individual*> &new_population, cmd_line &options ) {
void reproduce( const double *fitness, vector<Individual> &population, vector<Individual> &new_population, cmd_line &options ) {

    // populate parent multinomial samplings
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc( population.size(), fitness ) ;

    /// iterate through all individuals and draw parents + recomb
    for ( int i = 0 ; i < population.size() ; i ++ ) {
    
        /// new individual
        Individual new_ind ;
        
        /// grab mom and dad
        int mom = gsl_ran_discrete( rng, g ) ;
        int dad = gsl_ran_discrete( rng, g ) ;

        // NOW sort, and by locus, so this might fix the recombination issue:
        vector<Gene*> momTrnas = population[i].getMaternal_trnas();
        vector<Gene*> dadTrnas = population[i].getPaternal_trnas();
        std::sort(momTrnas.begin() , momTrnas.end(), sortByLocus);
        std::sort(dadTrnas.begin() , dadTrnas.end(), sortByLocus);
               
        // get their chromosomes 
       /* transmit_chromosome( *population[mom], new_ind.getMaternal_trnas() ) ;
        *
        transmit_chromosome( *population[dad], new_ind.getPaternal_trnas() ) ; */
        vector<Gene*> newMomTrnas = population[i].getMaternal_trnas();
        vector<Gene*> newDadTrnas = population[i].getPaternal_trnas();
        transmit_chromosome( population[mom], newMomTrnas ) ;
        transmit_chromosome( population[dad], newDadTrnas ) ;

        // swap individual in
        //swap( *new_population[i], new_ind ) ;
        swap( new_population[i], new_ind ) ;

    }

    // de-allocate space for multinomial samplings
    if (g) {
        gsl_ran_discrete_free( g ) ;
    }
}

#endif
