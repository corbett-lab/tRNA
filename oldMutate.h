#ifndef __MUTATE_H
#define __MUTATE_H


// sort container by locus
bool sortByLocus(gene* a, gene* b) { return (a->locus < b->locus); }

void mutate( vector<individual> &population, cmd_line &options, list<gene*> &trna_bank ) {
    
    // add germline mutations

	/// obviously we will want these mutations to have a range of possible effects on function
	/// bryan expect to modify this. 
	/// the current version simply makes the function 0 if there's a mutation 
    //
    // BRYAN NOTES:
    // - added assign_function.h, which assigns new function to tRNAs with mutations
    // - added neighborhood.h to give new attributes to tRNAs from jumping duplications
    // TODO: add in distribution of mutations
    // - changed locus difference on duplications to +1e-3 for visibility in printing
    // - changed sort function to hopefully resolve bug in recombination
    // - removing tRNA from the genome if its function is < 0.5 -- still there but no longer a functional tRNA
    //
    // TODO: update fitness function such
    //
    //
    // ISSUES:
    // - fitness must have something against total duplicate genes, otherwise will just keep growing
    // - new class of pseudogenes? tRNAs that are there but do not add to fitness. important for comparing to real data!
    // ^ could just set function to zero when it's below a certain point and make fitness function take that into account
    // - *** find consensus fitness function or come up with good one
    // - locus should never go negative or above maximum! must account for this
    // - each duplicate gene should add some fitness, to an extent
    // - how to model without biasing some ideal number of tRNAs to have
    // 

    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).germline ) ) {
            	
            	gene* new_trna = new gene ; 
            	new_trna->locus = (*population[i].maternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].maternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].maternal_trnas[g]).germline ;
            	new_trna->neighborhood = (*population[i].maternal_trnas[g]).neighborhood ;
            	// new_trna->function = 0 ; //(*population[i].maternal_trnas[g]).function ; 
                assign_function( population[i].maternal_trnas[g], new_trna ) ;

                // only add tRNA to bank if it's functional, otherwise not a tRNA
                trna_bank.push_back( new_trna ) ;
                population[i].paternal_trnas[g] = new_trna ; 
            }
        }

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).germline ) ) {
            	gene* new_trna = new gene ; 
            	new_trna->locus = (*population[i].paternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].paternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].paternal_trnas[g]).germline ;
            	new_trna->neighborhood = (*population[i].paternal_trnas[g]).neighborhood ;
            	// new_trna->function = 0 ; //(*population[i].maternal_trnas[g]).function ; 
                assign_function( population[i].paternal_trnas[g], new_trna ) ;

                // only add tRNA to bank if it's functional, otherwise not a tRNA
                trna_bank.push_back( new_trna ) ;
                population[i].paternal_trnas[g] = new_trna ; 
            }
        }
    }

    // delete tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = population[i].maternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + g ) ;
            }
        }
        for ( int g = population[i].paternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + g ) ;
            }
        }
    }

    // duplicate tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                gene* new_trna = new gene ; 
                //if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
                    new_trna->locus = (*population[i].maternal_trnas[g]).locus + gsl_rng_uniform( rng ) ;
                //}
                //else {
                //    new_trna->locus = (*population[i].maternal_trnas[g]).locus - gsl_rng_uniform( rng ) ;
                //}
                // cout << "DUPLICATION EVENT\t" << (population[i].maternal_trnas[g]) << "\t" << new_trna << "\t" ;
                // cout << (population[i].maternal_trnas[g])->locus << "\t" << new_trna->locus << endl ;
            	new_trna->somatic = (*population[i].maternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].maternal_trnas[g]).germline ;
            	new_trna->function = (*population[i].maternal_trnas[g]).function ;
            	new_trna->neighborhood = (*population[i].maternal_trnas[g]).neighborhood ; 

				/// if the trna jumps, need to reassigned neighborhood etc
            	if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    neighborhood( population[i].maternal_trnas[g], new_trna, options ) ;


            		//new_trna->locus = gsl_rng_uniform( rng ) * options.map_length ; 
            		//new_trna->function
            		//new_trna->neighborhood
            		//new_trna->somatic 
            		//new_trna->germline 
            	}

            	/// may want to have this be reset + reset locus with some probability
            	trna_bank.push_back( new_trna ) ; 
                (population[i].maternal_trnas).push_back( new_trna ) ; 
            }
        }

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                gene* new_trna = new gene ; 
            	//if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
                    new_trna->locus = (*population[i].maternal_trnas[g]).locus + gsl_rng_uniform( rng ) ;
                //}
                //else {
                //    new_trna->locus = (*population[i].maternal_trnas[g]).locus - gsl_rng_uniform( rng ) ;
                //}
            	new_trna->somatic = (*population[i].paternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].paternal_trnas[g]).germline ;
            	new_trna->function = (*population[i].paternal_trnas[g]).function ; 
            	new_trna->neighborhood = (*population[i].paternal_trnas[g]).neighborhood ; 

            	/// if the trna jumps, need to reassigned neighborhood etc
            	if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    neighborhood( population[i].paternal_trnas[g], new_trna, options ) ;
            	}

            	/// may want to have this be reset + reset locus with some probability
            	trna_bank.push_back( new_trna ) ; 
                (population[i].paternal_trnas).push_back( new_trna ) ; 
            }
        }
    }



    for ( int i = 0 ; i < population.size() ; i ++ ) {

        // NOW sort, and by locus, so this might fix the recombination issue:
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);

    }
}

#endif
