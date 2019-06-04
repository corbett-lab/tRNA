#ifndef __REDUND_FITNESS_H
#define __REDUND_FITNESS_H

/// compute fitness
void compute_redund_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

    //////////////////////////////////
    //////////////////////////////////
    ///////// REGULAR METHOD ///////// 
    //////////////////////////////////
    //////////////////////////////////
    //
    // tRNA MODEL BASED ON OBSERVED DATA:
    // fitness = fitness_seq * fitness_exp
    // fitness_seq = sequence (already based off of fitness values!)
    // fitness_exp = 1 - e^(-10.72 * x) // this is from the SPT15 data from yeast mapping expression to fitness
    // ^^ note: this assumes that one tRNA is enough. this is worth testing as an alternative hypothesis, where the gaussian is our null.

    /// needs to be redone:
    // calculate SUM of expression across entire tRNA gene set
    //
    // start with fitness 0.0
    // go across all tRNAs
    // find relationship between breadth of expression to depth
    // -- both are important, should not only use depth
    // -- also do not have that many tissues with definitive depth
    // multiply depth of expression by sequence and add that to your fitness?

    for ( int i = 0 ; i < population.size() ; ++i ) {
        double total_function = 0.0 ;

        // maternal fitness block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                    if ( options.mutation_pathways == false ){
                        int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
                        total_function += (((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * (*population[i].maternal_trnas[g]).expression) ;
                    }
                    else{
                        int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                        total_function += ( ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
                    }
                }
            }
            else {
                total_function += (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ;
            }
        }

        // paternal fitness block
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                    if ( options.mutation_pathways == false ){
                        int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
                        total_function += (((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * (*population[i].paternal_trnas[g]).expression) ;
                    }
                    else{
                        int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                        total_function += ( ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
                    }
                }
            }
            else {
                total_function += (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ;
            }
        }

        if ( total_function >= options.min_fitness ){
            fitness[i] = 1.0 ;
        }
        else {
            fitness[i] = total_function / options.min_fitness ;
        }
        //cout << total_function << endl ;
        if ( fitness[i] > 1.0 ) {
            cout << "FITNESS > 1.0:\t" << fitness[i] << "\tFITNESS NOT ALLOWED TO BE > 1.0. EXITING PROGRAM." << endl ;
            exit(0);
        } 
    }
    
}

#endif