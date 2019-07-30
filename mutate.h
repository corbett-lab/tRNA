#ifndef __MUTATE_H
#define __MUTATE_H


// sort container by locus
bool sortByLocus(gene* a, gene* b) { return (a->locus < b->locus); }

void mutate( vector<individual> &population, cmd_line &options, vector<gene*> &trna_bank, int current_gen, int &trna_counter, std::map<int,vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses ) {
    
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
    for ( int i = 0 ; i < population.size() ; ++i ) {
    //tbb::parallel_for ( size_t(0) , size_t(population.size()), [&] (size_t i) ) {

        // maternal germline block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).germline ) ) {
                //// new gene has the same attributes as progenitor except sequence
                //// assign_sequence resolves this with penalties
            	gene* new_trna = ::new gene ; 
            	new_trna->locus = (*population[i].maternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].maternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].maternal_trnas[g]).germline ;
            	new_trna->expression = (*population[i].maternal_trnas[g]).expression ;
                // germline mutations not truly "new" tRNAs -- take birth of progenitor
                new_trna->birth = (*population[i].maternal_trnas[g]).birth ;
                new_trna->progenitor = (*population[i].maternal_trnas[g]).name ;
                new_trna->birth_mode = 'g' ;
                if ( options.mutation_pathways == false ){
                    assign_sequence( population[i].maternal_trnas[g], new_trna, mutations_to_function, options ) ;
                }
                else {
                    assign_genotype( population[i].maternal_trnas[g], new_trna, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
                }
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
            	gene* new_trna = ::new gene ; 
            	new_trna->locus = (*population[i].paternal_trnas[g]).locus ; 
            	new_trna->somatic = (*population[i].paternal_trnas[g]).somatic ; 
            	new_trna->germline = (*population[i].paternal_trnas[g]).germline ;
            	new_trna->expression = (*population[i].paternal_trnas[g]).expression ;
                // germline mutations not truly "new" tRNAs -- take birth of progenitor
                new_trna->birth = (*population[i].paternal_trnas[g]).birth ;
                new_trna->progenitor = (*population[i].paternal_trnas[g]).name ;
                new_trna->birth_mode = 'g' ;
                if ( options.mutation_pathways == false ){
                    assign_sequence( population[i].paternal_trnas[g], new_trna, mutations_to_function, options ) ;
                }
                else {
                    assign_genotype( population[i].paternal_trnas[g], new_trna, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
                }
                trna_counter ++ ;
                new_trna->name = trna_counter ;
                trna_bank.push_back( new_trna ) ;
                population[i].paternal_trnas[g] = new_trna ; 
            }
        }
    }

    // duplications
    for ( int i = 0 ; i < population.size() ; ++i ) {

        // maternal duplication block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                gene* new_trna = ::new gene ; 
                // locus, expression, somatic, germline and sequence are all set in local.h and nonlocal.h
                new_trna->birth = current_gen ;
                new_trna->progenitor = (*population[i].maternal_trnas[g]).name ;

                /// ensure that loci don't double up
                std::map<double, int> temp_loci ;
                for ( int m = 0 ; m < population[i].maternal_trnas.size() ; m ++ ){
                    if ( !temp_loci.count( (population[i].maternal_trnas[m])->locus ) ) {
                        temp_loci.insert( std::pair<float,int> ( (population[i].maternal_trnas[m])->locus, 0 ) ) ;
                    }
                    else {
                        cout << "\t" << "ONE LOCUS HAS MULTIPLE TRNAS IN A SINGLE INDIVIDUAL. EXITING PROGRAM." << endl ;
                        cout << "GENERATION: " << current_gen << endl ;
                        for ( int z = 0 ; z < population[i].maternal_trnas.size() ; z ++ ){
                            cout.precision(10) ;
                            cout << "\t" << (population[i].maternal_trnas[z])->name << "_" << (population[i].maternal_trnas[z])->locus << "_" << (population[i].maternal_trnas[z])->sequence ;
                            cout << "_" << (population[i].maternal_trnas[z])->expression << "_" << (population[i].maternal_trnas[z])->muts << "_" << (population[i].maternal_trnas[z])->birth ;
                            cout << "_" << (population[i].maternal_trnas[z])->progenitor << "_" << (population[i].maternal_trnas[z])->birth_mode ;
                        }
                        exit(0);
                    }
                }

				/// 3 possible dups: nonlocal, local and seg dup -- all 3 maintain sequence and muts of original trna
                // nonlocal - retro-transposition - jump to a random locus in genome, get new expression randomly as function of new locus
                // local - unequal crossing-over - new locus and expression close to progenitor
                // segdup - locus very close to progenitor and expression exact same as progenitor
                if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    new_trna->birth_mode = 'n' ;
                    nonlocal( population[i].maternal_trnas[g], new_trna, temp_loci, options ) ;
                }
                else {
                    if ( gsl_ran_bernoulli( rng, 1 - options.prob_segdup ) ) {
                        new_trna->birth_mode = 'l' ;
                        local( population[i].maternal_trnas[g], new_trna, temp_loci, options ) ;
                    }
                    else {
                        new_trna->birth_mode = 's' ;
                        segdup( population[i].maternal_trnas[g], new_trna, temp_loci, options ) ;
                    }
                }
                temp_loci.clear() ;
                trna_counter ++ ;
                new_trna->name = trna_counter ;
            	trna_bank.push_back( new_trna ) ; 
                (population[i].maternal_trnas).push_back( new_trna ) ; 
            }
        }

        // paternal duplication block
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) {
                gene* new_trna = ::new gene ; 
                // locus, expression, somatic, germline and sequence are all set in local.h and nonlocal.h
                new_trna->birth = current_gen ;
                new_trna->progenitor = (*population[i].paternal_trnas[g]).name ;

                /// ensure that loci don't double up
                std::map<double, int> temp_loci ;
                for ( int p = 0 ; p < population[i].paternal_trnas.size() ; p ++ ){
                    if ( !temp_loci.count( (population[i].paternal_trnas[p])->locus ) ) {
                        temp_loci.insert( std::pair<float,int> ( (population[i].paternal_trnas[p])->locus, 0 ) ) ;
                    }
                    else {
                        cout << "\t" << "ONE LOCUS HAS MULTIPLE TRNAS IN A SINGLE INDIVIDUAL. EXITING PROGRAM." << endl ;
                        cout << "GENERATION: " << current_gen << endl ;
                        for ( int z = 0 ; z < population[i].paternal_trnas.size() ; z ++ ){
                            cout.precision(10) ;
                            cout << "\t" << (population[i].paternal_trnas[z])->name << "_" << (population[i].paternal_trnas[z])->locus << "_" << (population[i].paternal_trnas[z])->sequence ;
                            cout << "_" << (population[i].paternal_trnas[z])->expression << "_" << (population[i].paternal_trnas[z])->muts << "_" << (population[i].paternal_trnas[z])->birth ;
                            cout << "_" << (population[i].paternal_trnas[z])->progenitor << "_" << (population[i].paternal_trnas[z])->birth_mode ;
                        }
                        exit(0);
                    }
                }

                // same as above
                if ( gsl_ran_bernoulli( rng, 1 - options.prob_cluster ) ) {
                    new_trna->birth_mode = 'n' ;
                    nonlocal( population[i].paternal_trnas[g], new_trna, temp_loci, options ) ;
                }
                else {
                    if ( gsl_ran_bernoulli( rng, 1 - options.prob_segdup ) ) {
                        new_trna->birth_mode = 'l' ;
                        local( population[i].paternal_trnas[g], new_trna, temp_loci, options ) ;
                    }
                    else {
                        new_trna->birth_mode = 's' ;
                        segdup( population[i].paternal_trnas[g], new_trna, temp_loci, options ) ;
                    }
                }
                temp_loci.clear() ;
                trna_counter ++ ;
                new_trna->name = trna_counter ;
                trna_bank.push_back( new_trna ) ; 
                (population[i].paternal_trnas).push_back( new_trna ) ; 
            }
        }

    }

    // delete tRNAs
    for ( int i = 0 ; i < population.size() ; ++i ) {
        // maternal deletion block
        for ( int g = population[i].maternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( population[i].maternal_trnas[g]->sequence == 0.0 ) {
                population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + g ) ;
            }
            else if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + g ) ;
            }
        }
        // paternal deletion block
        for ( int g = population[i].paternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( population[i].paternal_trnas[g]->sequence == 0.0 ) {
                population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + g ) ;
            }
            else if ( gsl_ran_bernoulli( rng, options.deletion_rate ) ) {
                population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + g ) ;
            }
        }
    }

    // sort tRNAs by locus for ease in recombination
    for ( int i = 0 ; i < population.size() ; ++i ) {
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);
    }
}

#endif
