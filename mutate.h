#ifndef __MUTATE_H
#define __MUTATE_H


// sort container by locus
bool sortByLocus(gene* a, gene* b) { return (a->locus < b->locus); }

void mutate( vector<individual> &population, cmd_line &options, list<gene*> &trna_bank, int current_gen, int &trna_counter, vector<double> &mutation_penalties ) {
    
    // add germline mutations
    //
    /////////////////////////// BRYAN NOTES:
    // - added assign_sequence.h, which assigns new sequence to tRNAs with mutations
    // - added nonlocal.h and local.h to give new attributes to tRNAs from jumping duplications
    /////////////// TODO: add in distribution of mutations: done
    // - changed sort sequence to hopefully resolve bug in recombination: done
    // - mutations that make things NOT a tRNA despite bit score change? advantageous mutations?
    //
    ///////////// MISC EARLY THOUGHTS I DON'T WANT TO COMPLETELY DELETE:
    // - fitness must have something against total duplicate genes, otherwise will just keep growing
    // - new class of pseudogenes? tRNAs that are there but do not add to fitness. important for comparing to real data!
    // ^ could just set sequence to zero when it's below a certain point and make fitness sequence take that into account
    // - *** find consensus fitness sequence or come up with good one
    // - locus should never go negative or above maximum! must account for this
    // - each duplicate gene should add some fitness, to an extent
    // - how to model without biasing some ideal number of tRNAs to have
    // 

    // germline mutations
    for ( int i = 0 ; i < population.size() ; i ++ ) {

        // maternal germline block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).germline ) ) {
                //// new gene has the same attributes as progenitor except sequence
                //// assign_sequence resolves this with penalties
            	gene* new_trna = new gene ; 
            	new_trna->locus = (*population[i].maternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].maternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].maternal_trnas[g]).germline ;
            	new_trna->expression = (*population[i].maternal_trnas[g]).expression ;
                new_trna->birth = current_gen ;
                new_trna->frequency.push_back( 0 ) ;
                new_trna->progenitor = (*population[i].maternal_trnas[g]).name ;
                assign_sequence( population[i].maternal_trnas[g], new_trna, mutation_penalties, options ) ;
                trna_counter ++ ;
                new_trna->name = trna_counter ;
                trna_bank.push_back( new_trna ) ; 
                population[i].maternal_trnas[g] = new_trna ;
            }
        }

        // paternal germline block
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).germline ) ) {
                //// new gene has the same attributes as progenitor except sequence
                //// assign_sequence resolves this with penalties
            	gene* new_trna = new gene ; 
            	new_trna->locus = (*population[i].paternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].paternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].paternal_trnas[g]).germline ;
            	new_trna->expression = (*population[i].paternal_trnas[g]).expression ;
                new_trna->birth = current_gen ;
                new_trna->frequency.push_back( 0 ) ;
                new_trna->progenitor = (*population[i].paternal_trnas[g]).name ;
                assign_sequence( population[i].paternal_trnas[g], new_trna, mutation_penalties, options ) ;
                trna_counter ++ ;
                new_trna->name = trna_counter ;
                trna_bank.push_back( new_trna ) ;
                population[i].paternal_trnas[g] = new_trna ; 
            }
        }
    }

    // duplications
    for ( int i = 0 ; i < population.size() ; i ++ ) {

        // maternal duplication block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                gene* new_trna = new gene ; 
                // locus, expression, somatic, germline and sequence are all set in local.h and nonlocal.h
                new_trna->birth = current_gen ;
                new_trna->frequency.push_back( 0 ) ;
                new_trna->progenitor = (*population[i].maternal_trnas[g]).name ;

				/// if the trna jumps, need to reassigned neighborhood etc
            	if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    nonlocal( population[i].maternal_trnas[g], new_trna, options ) ;
            	}
                else {
                    local( population[i].maternal_trnas[g], new_trna, options ) ;
                }
                trna_counter ++ ;
                new_trna->name = trna_counter ;
            	trna_bank.push_back( new_trna ) ; 
                (population[i].maternal_trnas).push_back( new_trna ) ; 
            }
        }

        // paternal duplication block
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                gene* new_trna = new gene ; 
                // locus, expression, somatic, germline and sequence are all set in local.h and nonlocal.h
                new_trna->birth = current_gen ;
                new_trna->frequency.push_back( 0 ) ;
                new_trna->progenitor = (*population[i].paternal_trnas[g]).name ;

                /// if the trna jumps, need to reassigned neighborhood etc
                if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    nonlocal( population[i].paternal_trnas[g], new_trna, options ) ;
                }
                else {
                    local( population[i].paternal_trnas[g], new_trna, options ) ;
                }
                trna_counter ++ ;
                new_trna->name = trna_counter ;
                trna_bank.push_back( new_trna ) ; 
                (population[i].paternal_trnas).push_back( new_trna ) ; 
            }
        }

    }

    // delete tRNAs
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        // maternal duplication block
        for ( int g = population[i].maternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + g ) ;
            }
        }
        // paternal duplication block
        for ( int g = population[i].paternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + g ) ;
            }
        }
    }

    // sort tRNAs by locus for ease in recombination
    for ( int i = 0 ; i < population.size() ; i ++ ) {
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);
    }
}

#endif
