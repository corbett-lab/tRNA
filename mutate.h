#ifndef __MUTATE_H
#define __MUTATE_H

void update_gene(gene* new_trna,vector<gene*> trnas, int i, int g,int current_gen, cmd_line &options,int &trna_counter,list<gene*> &trna_bank, bool duplicate) {
	if (duplicate) {
        if ( gsl_ran_bernoulli( rng, 0.5 ) or ((trnas[g])->locus + 1 > options.map_length) ) {
              new_trna->locus = (trnas[g])->locus - gsl_rng_uniform( rng ) ;
          }
          else{
              new_trna->locus = (trnas[g])->locus + gsl_rng_uniform( rng ) ;
          }
    	/// if the trna jumps, need to reassigned neighborhood etc
        if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
        	neighborhood(trnas[g], new_trna, options ) ;
        }
	}
	else {
		new_trna->locus = trnas[g]->locus ;
	}
    new_trna->somatic = trnas[g]->somatic ;
    new_trna->germline = trnas[g]->germline ;
    new_trna->neighborhood = trnas[g]->neighborhood ;
    new_trna->birth = current_gen ;
    new_trna->frequency.push_back( 0 ) ;
    new_trna->progenitor = trnas[g]->name ;
    assign_function( trnas[g], new_trna, options ) ;
    trna_counter ++ ;
    new_trna->name = trna_counter ;
    trna_bank.push_back(new_trna ) ;
    trnas[g] = new_trna ;
}

// sort container by locus
bool sortByLocus(gene* a, gene* b) { return (a->locus < b->locus); }

void mutate( vector<individual> &population, cmd_line &options, list<gene*> &trna_bank, int current_gen, int &trna_counter) {
    
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
    // - mutations that make things NOT a tRNA despite bit score change? advantageous mutations?
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

	// TODO merge the duplicated code in 4 for loops!



    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).germline ) ) {
            	gene* new_trna = new gene ;
            	update_gene(new_trna, population[i].maternal_trnas,  i,  g, current_gen, options,trna_counter,trna_bank, false);
            }
        }

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).germline ) ) {
            	gene* new_trna = new gene ; 
            	update_gene(new_trna, population[i].paternal_trnas,  i,  g, current_gen, options,trna_counter,trna_bank, false);

            }
        }
    }

    // duplicate tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                gene* new_trna = new gene ; 
            	update_gene(new_trna, population[i].maternal_trnas,  i,  g, current_gen, options,trna_counter,trna_bank, true);

            }
        }

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                gene* new_trna = new gene ; 
            	update_gene(new_trna, population[i].paternal_trnas,  i,  g, current_gen, options,trna_counter,trna_bank, true);

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



    for ( int i = 0 ; i < population.size() ; i ++ ) {

        // NOW sort, and by locus, so this might fix the recombination issue:
        // NO, next is fitness so we are going to re-sort by function anyway
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);
    }
}

#endif
