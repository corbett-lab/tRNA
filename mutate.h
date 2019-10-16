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

    /// NEW IDEA: only tRNAs that are not orthologs will get a new name

    //////////////////////////////
    ///// GERMLINE MUTATIONS /////
    //////////////////////////////

    for ( int i = 0 ; i < population.size() ; ++i ) {

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
                // germline mutations not truly "new" tRNAs -- same name
                new_trna->birth = (*population[i].maternal_trnas[g]).birth ;
                new_trna->progenitor = (*population[i].maternal_trnas[g]).progenitor ;
                new_trna->name = (*population[i].maternal_trnas[g]).name ;
                new_trna->birth_mode = 'g' ;
                if ( options.dual_rates == true ){
                    assign_genotype_gamma( population[i].maternal_trnas[g], new_trna, options ) ;
                }
                else if ( options.mutation_pathways == false ){
                    assign_genotype_model( population[i].maternal_trnas[g], new_trna, mutations_to_function, options ) ;
                }
                else {
                    assign_genotype_pathways( population[i].maternal_trnas[g], new_trna, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
                }
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
                // germline mutations not truly "new" tRNAs -- same name
                new_trna->birth = (*population[i].paternal_trnas[g]).birth ;
                new_trna->progenitor = (*population[i].paternal_trnas[g]).progenitor ;
                new_trna->name = (*population[i].paternal_trnas[g]).name ;
                new_trna->birth_mode = 'g' ;
                if ( options.dual_rates == true ){
                    assign_genotype_gamma( population[i].paternal_trnas[g], new_trna, options ) ;
                }
                else if ( options.mutation_pathways == false ){
                    assign_genotype_model( population[i].paternal_trnas[g], new_trna, mutations_to_function, options ) ;
                }
                else {
                    assign_genotype_pathways( population[i].paternal_trnas[g], new_trna, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
                }
                trna_bank.push_back( new_trna ) ;
                population[i].paternal_trnas[g] = new_trna ; 
            }
        }
    }

    ////////////////////////
    ///// DUPLICATIONS /////
    ////////////////////////

    for ( int i = 0 ; i < population.size() ; ++i ) {

        // duplicated genes should have a new name

        // segmental duplications
        // 3 thousand centimorgans in map length unscaled
        // 30 morgans map length unscaled
        // geometric will give length in base pairs
        // whatever exponential throws out over 3 billion = what we want over map_length
        // so take exponential, divide by 3 billion and multiply by map_length
        // rate for non-local duplications is per gene, so multiply by number of genes in that genome
        // MATERNAL BLOCK:

        // maternal duplication block
        int temp_size_m = population[i].maternal_trnas.size() ;
        for ( int g = 0 ; g < temp_size_m ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) { 
                if ( gsl_ran_bernoulli( rng, options.prob_cluster ) ) {
                    double my_start = population[i].maternal_trnas[g]->locus - (options.map_length * 0.0000001) ;
                    double my_end = my_start + gsl_ran_exponential(rng, options.segdup_exp) ;
                    for ( int t = 0 ; t < temp_size_m ; t ++ ) {
                        /// SEGMENTAL DUPLICATION:
                        /// if not the tRNA in question, leave it alone
                        /// if tRNA(s) being duplicated, add duplicate copy at new location
                        /// if it is a tRNA from a previous seg dup, move it up
                        if ( population[i].maternal_trnas[t]->locus > my_start ){
                            if ( population[i].maternal_trnas[t]->locus < my_end ){
                                gene* new_trna = ::new gene ; 
                                new_trna->locus = (population[i].maternal_trnas[t])->locus + (my_end-my_start) ; 
                                new_trna->somatic = (population[i].maternal_trnas[t])->somatic ; 
                                new_trna->germline = (population[i].maternal_trnas[t])->germline ;
                                new_trna->sequence = (population[i].maternal_trnas[t])->sequence ;
                                new_trna->genotype = (population[i].maternal_trnas[t])->genotype ;
                                new_trna->muts = (population[i].maternal_trnas[t])->muts ;
                                new_trna->expression = (population[i].maternal_trnas[t])->expression ;
                                new_trna->birth = (population[i].maternal_trnas[t])->birth ;
                                new_trna->birth_mode = 'S' ;
                                trna_bank.push_back( new_trna ) ;
                                new_trna->progenitor = (population[i].maternal_trnas[t])->name ;
                                trna_counter ++ ;
                                new_trna->name = trna_counter ;
                                population[i].maternal_trnas.push_back( new_trna ) ; 
                            }
                            // don't bother moving up unless the tRNA in question was itself produced by a seg dup
                            else if ( population[i].maternal_trnas[t]->birth_mode == 'S' ){
                                gene* new_trna = ::new gene ; 
                                new_trna->locus = (population[i].maternal_trnas[t])->locus + (my_end-my_start) ; 
                                new_trna->somatic = (population[i].maternal_trnas[t])->somatic ; 
                                new_trna->germline = (population[i].maternal_trnas[t])->germline ;
                                new_trna->sequence = (population[i].maternal_trnas[t])->sequence ;
                                new_trna->genotype = (population[i].maternal_trnas[t])->genotype ;
                                new_trna->muts = (population[i].maternal_trnas[t])->muts ;
                                new_trna->expression = (population[i].maternal_trnas[t])->expression ;
                                new_trna->birth = (population[i].maternal_trnas[t])->birth ;
                                new_trna->birth_mode = 'S' ;
                                trna_bank.push_back( new_trna ) ;
                                new_trna->progenitor = (population[i].maternal_trnas[t])->name ;
                                trna_counter ++ ;
                                new_trna->name = trna_counter ;
                                population[i].maternal_trnas[t] = new_trna ; 
                            }
                        }
                    }
                }

                else {
                    /// otherwise, a nonlocal duplication
                    gene* new_trna = ::new gene ; 
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
                                cout.precision(15) ;
                                cout << "\t" << (population[i].maternal_trnas[z])->name << "_" << (population[i].maternal_trnas[z])->locus << "_" << (population[i].maternal_trnas[z])->sequence ;
                                cout << "_" << (population[i].maternal_trnas[z])->expression << "_" << (population[i].maternal_trnas[z])->muts << "_" << (population[i].maternal_trnas[z])->birth ;
                                cout << "_" << (population[i].maternal_trnas[z])->progenitor << "_" << (population[i].maternal_trnas[z])->birth_mode ;
                            }
                            exit(0);
                        }
                    }
                    new_trna->birth_mode = 'n' ;
                    nonlocal( population[i].maternal_trnas[g], new_trna, temp_loci, options ) ; 
                    temp_loci.clear() ;
                    trna_counter ++ ;
                    new_trna->name = trna_counter ;
                	trna_bank.push_back( new_trna ) ; 
                    (population[i].maternal_trnas).push_back( new_trna ) ; 
                }
            }
        }

        // paternal duplication block
        int temp_size_p = population[i].paternal_trnas.size() ;
        for ( int g = 0 ; g < temp_size_p ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, options.duplication_rate ) ) { 
                if ( gsl_ran_bernoulli( rng, options.prob_cluster ) ) {
                    double my_start = population[i].paternal_trnas[g]->locus - (options.map_length * 0.0000001) ;
                    double my_end = my_start + gsl_ran_exponential(rng, options.segdup_exp) ;
                    for ( int t = 0 ; t < temp_size_p ; t ++ ) {
                        if ( population[i].paternal_trnas[t]->locus > my_start ){
                            gene* new_trna = ::new gene ; 
                            new_trna->locus = (population[i].paternal_trnas[t])->locus + (my_end-my_start) ; 
                            new_trna->somatic = (population[i].paternal_trnas[t])->somatic ; 
                            new_trna->germline = (population[i].paternal_trnas[t])->germline ;
                            new_trna->sequence = (population[i].paternal_trnas[t])->sequence ;
                            new_trna->genotype = (population[i].paternal_trnas[t])->genotype ;
                            new_trna->muts = (population[i].paternal_trnas[t])->muts ;
                            new_trna->expression = (population[i].paternal_trnas[t])->expression ;
                            new_trna->birth = (population[i].paternal_trnas[t])->birth ;
                            new_trna->birth_mode = 'S' ;
                            trna_bank.push_back( new_trna ) ;
                            if ( population[i].paternal_trnas[t]->locus < my_end ){
                                new_trna->progenitor = (population[i].paternal_trnas[t])->name ;
                                trna_counter ++ ;
                                new_trna->name = trna_counter ;
                                population[i].paternal_trnas.push_back( new_trna ) ; 
                            }
                            else {
                                new_trna->progenitor = (population[i].paternal_trnas[t])->progenitor ;
                                new_trna->name = (population[i].paternal_trnas[t])->name ;
                                population[i].paternal_trnas[t] = new_trna ;
                            }
                        }
                    }
                }

                else {
                    gene* new_trna = ::new gene ; 
                    new_trna->birth = current_gen ;
                    new_trna->progenitor = (*population[i].paternal_trnas[g]).name ;

                    /// ensure that loci don't double up
                    std::map<double, int> temp_loci ;
                    for ( int m = 0 ; m < population[i].paternal_trnas.size() ; m ++ ){
                        if ( !temp_loci.count( (population[i].paternal_trnas[m])->locus ) ) {
                            temp_loci.insert( std::pair<float,int> ( (population[i].paternal_trnas[m])->locus, 0 ) ) ;
                        }
                        else {
                            cout << "\t" << "ONE LOCUS HAS MULTIPLE TRNAS IN A SINGLE INDIVIDUAL. EXITING PROGRAM." << endl ;
                            cout << "GENERATION: " << current_gen << endl ;
                            for ( int z = 0 ; z < population[i].paternal_trnas.size() ; z ++ ){
                                cout.precision(15) ;
                                cout << "\t" << (population[i].paternal_trnas[z])->name << "_" << (population[i].paternal_trnas[z])->locus << "_" << (population[i].paternal_trnas[z])->sequence ;
                                cout << "_" << (population[i].paternal_trnas[z])->expression << "_" << (population[i].paternal_trnas[z])->muts << "_" << (population[i].paternal_trnas[z])->birth ;
                                cout << "_" << (population[i].paternal_trnas[z])->progenitor << "_" << (population[i].paternal_trnas[z])->birth_mode ;
                            }
                            exit(0);
                        }
                    }
                    new_trna->birth_mode = 'n' ;
                    nonlocal( population[i].paternal_trnas[g], new_trna, temp_loci, options ) ; 
                    temp_loci.clear() ;
                    trna_counter ++ ;
                    new_trna->name = trna_counter ;
                    trna_bank.push_back( new_trna ) ; 
                    (population[i].paternal_trnas).push_back( new_trna ) ; 
                }
            }
        }
    }


    /////////////////////
    ///// DELETIONS /////
    /////////////////////

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
