#ifndef __MUTATE_H
#define __MUTATE_H

void update_gene(Gene* new_trna,vector<Gene*> trnas, int i, int g,int current_gen, cmd_line &options,int &trna_counter,list<Gene*> &trna_bank, bool duplicate) {
	if (duplicate) {
        if ( gsl_ran_bernoulli( rng, 0.5 ) or ((trnas[g])->getLocus() + 1 > options.map_length) ) {
              new_trna->setLocus((trnas[g])->getLocus() - gsl_rng_uniform( rng )) ;
          }
          else{
              new_trna->setLocus((trnas[g])->getLocus() + gsl_rng_uniform( rng )) ;
          }
    	/// if the trna jumps, need to reassigned neighborhood etc
        if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
        	neighborhood(trnas[g], new_trna, options ) ;
        }
	}
	else {
		new_trna->setLocus(trnas[g]->getLocus()) ;
	}
    new_trna->setSomatic(trnas[g]->getSomatic()) ;
    new_trna->setGermline(trnas[g]->getGermline()) ;
    new_trna->setNeighborhood(trnas[g]->getNeighborhood()) ;
    new_trna->setBirth(current_gen) ;
    /*vector<int> f = new_trna->getFrequency();
    f.push_back(0);
    new_trna->setFrequency(f);*/
    new_trna->getFrequency().push_back(0);
    new_trna->setProgenitor(trnas[g]->getName()) ;
    assign_function( trnas[g], new_trna, options ) ;
    trna_counter ++ ;
    new_trna->setName(trna_counter) ;
    trna_bank.push_back(new_trna ) ;
    trnas[g] = new_trna ;
}

// sort container by locus
bool sortByLocus(Gene* a, Gene* b) { return (a->getLocus() < b->getLocus()); }

//void mutate( vector<individual*> &population, cmd_line &options, list<gene*> &trna_bank, int current_gen, int &trna_counter) {
void mutate( vector<Individual> &population, cmd_line &options, list<Gene*> &trna_bank, int current_gen, int &trna_counter) {

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



    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].getMaternal_trnas().size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].getMaternal_trnas()[g]).getGermline() ) ) {
            	Gene* new_trna = new Gene ;
            	update_gene(new_trna, population[i].getMaternal_trnas(),  i,  g, current_gen, options,trna_counter,trna_bank, false);
            }
        }

        for ( int g = 0 ; g < population[i].getPaternal_trnas().size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].getPaternal_trnas()[g]).getGermline() ) ) {
            	Gene* new_trna = new Gene ; 
            	update_gene(new_trna, population[i].getPaternal_trnas(),  i,  g, current_gen, options,trna_counter,trna_bank, false);

            }
        }
    }

    // duplicate tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = 0 ; g < population[i].getMaternal_trnas().size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                Gene* new_trna = new Gene ; 
            	update_gene(new_trna, population[i].getMaternal_trnas(),  i,  g, current_gen, options,trna_counter,trna_bank, true);

            }
        }

        for ( int g = 0 ; g < population[i].getPaternal_trnas().size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                
                /// copy current gene copy
                Gene* new_trna = new Gene ; 
            	update_gene(new_trna, population[i].getPaternal_trnas(),  i,  g, current_gen, options,trna_counter,trna_bank, true);

            }
        }
    }

    // delete tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        for ( int g = population[i].getMaternal_trnas().size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].getMaternal_trnas().erase( population[i].getMaternal_trnas().begin() + g ) ;
            }
        }
        for ( int g = population[i].getPaternal_trnas().size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].getPaternal_trnas().erase( population[i].getPaternal_trnas().begin() + g ) ;
            }
        }
    }



    for ( int i = 0 ; i < population.size() ; i ++ ) {

        // NOW sort, and by locus, so this might fix the recombination issue:
        // NO, next is fitness so we are going to re-sort by function anyway
        std::sort(population[i].getMaternal_trnas().begin() , population[i].getMaternal_trnas().end(), sortByLocus);
        std::sort(population[i].getPaternal_trnas().begin() , population[i].getPaternal_trnas().end(), sortByLocus);
    }
}

#endif
