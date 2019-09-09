#ifndef __FITNESS_H
#define __FITNESS_H

/// compute fitness

// actual functions
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

    for ( int i = 0 ; i < population.size() ; i ++ ) {
        double total_function = 0.0 ;

        // maternal fitness block
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                // somatic SNP
                if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                // if mutation destroys tRNA, draw number of tissues NOT affected (1 - uniform = uniform) times normal fitness effect of tRNA
                                total_function += (gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ;
                            }
                            else {
                                // if mutation doesn't destroy the tRNA, draw a gamma effect
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                // fitness effect += proportion of cells affected times tRNA with gamma effect
                                // also add in the rest of the tissues, which have a genotype without the somatic SNP
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].maternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].maternal_trnas[g]).expression) ;
                                total_function += temp_function + ((1.0 - temp_prop) * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            // draw mutation effect and proportion of tissues affected
                            int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            // add in tissues affected * mutation effect
                            // also add in rest of tissues times normal genotype
                            total_function += (temp_prop * ((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * (*population[i].maternal_trnas[g]).expression) ;
                            total_function += ((1.0 - temp_prop) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ;
                        }
                        else{
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            // add in tissues affected * mutation effect
                            // also add in rest of tissues times normal phenotype
                            total_function += (temp_prop * ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * ((*population[i].maternal_trnas[g]).expression)) ;
                            total_function += ((1.0 - temp_prop) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ;
                        }
                    }
                }
                // no somatic SNP
                else {
                    // just add normal genotype
                    total_function += (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ;
                }

                // somatic dup
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    // draw tissues affected
                    double temp_prop = gsl_rng_uniform(rng) ;
                    // assume a local duplication for simplicity, we have no idea how somatic dups are expressed
                    // add effect of dup times tissues affected
                    // don't need to add the rest because you've already accounted for normal genotype (under no somatic SNP)
                    total_function += (temp_prop * ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856))) ;
                }
            }
            // if deletion, same thing as a SNP that destroys all tRNA function
            else {
                total_function += ( gsl_rng_uniform(rng) * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) ) ;
            }
        }

        // paternal fitness block
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                // somatic SNP
                if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                // if mutation destroys tRNA, draw number of tissues NOT affected (1 - uniform = uniform) times normal fitness effect of tRNA
                                total_function += (gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ;
                            }
                            else {
                                // if mutation doesn't destroy the tRNA, draw a gamma effect
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                // fitness effect += proportion of cells affected times tRNA with gamma effect
                                // also add in the rest of the tissues, which have a genotype without the somatic SNP
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].paternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].paternal_trnas[g]).expression) ;
                                total_function += temp_function + ((1.0 - temp_prop) * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            // draw mutation effect and proportion of tissues affected
                            int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            // add in tissues affected * mutation effect
                            // also add in rest of tissues times normal genotype
                            total_function += (temp_prop * ((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * (*population[i].paternal_trnas[g]).expression) ;
                            total_function += ((1.0 - temp_prop) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ;
                        }
                        else{
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            // add in tissues affected * mutation effect
                            // also add in rest of tissues times normal phenotype
                            total_function += (temp_prop * ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * ((*population[i].paternal_trnas[g]).expression)) ;
                            total_function += ((1.0 - temp_prop) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ;
                        }
                    }
                }
                // no somatic SNP
                else {
                    // just add normal genotype
                    total_function += (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ;
                }

                // somatic dup
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    // draw tissues affected
                    double temp_prop = gsl_rng_uniform(rng) ;
                    // assume a local duplication for simplicity, we have no idea how somatic dups are expressed
                    // add effect of dup times tissues affected
                    // don't need to add the rest because you've already accounted for normal genotype (under no somatic SNP)
                    total_function += (temp_prop * ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856))) ;
                }
            }
            // if deletion, same thing as a SNP that destroys all tRNA function
            else {
                total_function += ( gsl_rng_uniform(rng) * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) ) ;
            }
        }

        fitness[i] = total_function / options.min_fitness ;

        if (fitness[i] > 1.0) {
            fitness[i] = 1.0 ;
        }
        if (fitness[i] < 0.0) {
            fitness[i] = 0.0 ;
        }
    }
}

/// compute fitness
void compute_model_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

    /////////////////////
    ////// MODEL 4 //////
    /////////////////////
    // alternative fitness function described in model4 paper:
    // if at least 1 functional tRNA, fitness = 1 - somatic rate ^ number of tRNAs
    // entirely dependent on sequence -- expression assumed the same
    // linear mapping of function to fitness
    //
    // in the models recall that all mutations set sequence to zero

    if ( options.model4 == true ) {
        for ( int i = 0 ; i < population.size() ; ++i ) {

            // get total function of all of an individual's tRNAs
            double mom_function = 0.0 ;
            double dad_function = 0.0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                mom_function += (*population[i].maternal_trnas[g]).sequence ;
                }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                dad_function += (*population[i].paternal_trnas[g]).sequence ;
                }

            // if they have no tRNAs, set their fitness to 0
            if ( (mom_function == 0) and (dad_function == 0) ){
                fitness[i] = 0.0 ;
            }

            // take whichever chromosome has more total function and take error rate raised to that power
            // (total function == number of tRNAs here because they're all 0.0 and 1.0 in this model)
            // (also originally a haploid model so this our adaptation to diploid)
            else if ( mom_function >= dad_function ){
                fitness[i] = 1.0 - ( pow(options.model4_deverr, mom_function) ) ;
                //cout << options.model4_deverr << "\t" << mom_function << "\t" << pow(options.model4_deverr, mom_function) << endl ;
            }
            else {
                fitness[i] = 1.0 - ( pow(options.model4_deverr, dad_function) ) ;
                //cout << options.model4_deverr << "\t" << dad_function << "\t" << pow(options.model4_deverr, dad_function) << endl ;
            }
        }
    }

    //////////////////////////
    ////// MODELS 1 & 2 //////
    //////////////////////////
    // in models 1 and 2, just take the best functioning tRNA and that's your fitness
    // in the models recall that all mutations set sequence to zero

    else if (( options.model1 == true ) or ( options.model2 == true )) {

        for ( int i = 0 ; i < population.size() ; ++i ) {
            double max_function = 0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                if ((*population[i].maternal_trnas[g]).sequence > max_function){
                    max_function = (*population[i].maternal_trnas[g]).sequence ;
                }
            }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                if ((*population[i].paternal_trnas[g]).sequence > max_function){
                    max_function = (*population[i].paternal_trnas[g]).sequence ;
                }
            }
            fitness[i] = max_function ;
        }
    }


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

    else {
        for ( int i = 0 ; i < population.size() ; ++i ) {
            double max_function = 0.0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                    if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                        if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                            double temp = (gsl_rng_uniform(rng)) ;
                            if ( options.dual_rates == true ){
                                // if mutation destroys tRNA, need proportion of tissues NOT affected times normal genotype
                                if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                    if (temp * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression > max_function){
                                        max_function = temp * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression ;
                                    }
                                }
                                else {
                                    double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                    if ( temp_gamma > 1.0 ){
                                        temp_gamma = 1.0 ;
                                    }
                                    // store function of mutation effect times proportion of tissues affected
                                    double temp_function = (temp * ((*population[i].maternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].maternal_trnas[g]).expression) ;
                                    // also get unaffected tissues * normal genotype
                                    if ( temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) > max_function ){
                                        max_function = (temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                                    }
                                }
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
                                // store function of mutation effect times proportion of tissues affected
                                double temp_function = (temp * ((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * (*population[i].maternal_trnas[g]).expression) ;
                                // also get unaffected tissues * normal genotype
                                if ( temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) > max_function ){
                                    max_function = (temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                                }
                            }
                            else{
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                                // store function of mutation effect times proportion of tissues affected
                                double temp_function = (temp * ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * (*population[i].maternal_trnas[g]).expression) ;
                                // also get unaffected tissues * normal genotype
                                if ( temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression) > max_function ){
                                    max_function = (temp_function + ((1.0 - temp) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                                }
                            }
                        } 
                    }
                    // no somatic SNP, just use normal genotype
                    else {
                        if (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression > max_function){
                            max_function = ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression ;
                        }
                    }
                    // somatic dup doesn't make much sense in this fitness function
                    // if we're only looking at your best tRNA, a copy of that tRNA won't be better than the one it is a copy of
                }
                // if deletion, same thing as a SNP that destroys all tRNA function
                else {
                    double temp = gsl_rng_uniform(rng) ;
                    if (temp * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression > max_function ) {
                        max_function = (temp * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression ) ;
                    }
                }
            }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                    if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                        if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                            double temp = (gsl_rng_uniform(rng)) ;
                            if ( options.dual_rates == true ){
                                // if mutation destroys tRNA, need proportion of tissues NOT affected times normal genotype
                                if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                    if (temp * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression > max_function){
                                        max_function = temp * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression ;
                                    }
                                }
                                else {
                                    double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                    if ( temp_gamma > 1.0 ){
                                        temp_gamma = 1.0 ;
                                    }
                                    // store function of mutation effect times proportion of tissues affected
                                    double temp_function = (temp * ((*population[i].paternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].paternal_trnas[g]).expression) ;
                                    // also get unaffected tissues * normal genotype
                                    if ( temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) > max_function ){
                                        max_function = (temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                                    }
                                }
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
                                // store function of mutation effect times proportion of tissues affected
                                double temp_function = (temp * ((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * (*population[i].paternal_trnas[g]).expression) ;
                                // also get unaffected tissues * normal genotype
                                if ( temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) > max_function ){
                                    max_function = (temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                                }
                            }
                            else{
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                                // store function of mutation effect times proportion of tissues affected
                                double temp_function = (temp * ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * (*population[i].paternal_trnas[g]).expression) ;
                                // also get unaffected tissues * normal genotype
                                if ( temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression) > max_function ){
                                    max_function = (temp_function + ((1.0 - temp) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                                }
                            }
                        } 
                    }
                    // no somatic SNP, just use normal genotype
                    else {
                        if (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression > max_function){
                            max_function = ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression ;
                        }
                    }
                    // somatic dup doesn't make much sense in this fitness function
                    // if we're only looking at your best tRNA, a copy of that tRNA won't be better than the one it is a copy of
                }
                // if deletion, same thing as a SNP that destroys all tRNA function
                else {
                    double temp = gsl_rng_uniform(rng) ;
                    if (temp * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression > max_function ) {
                        max_function = (temp * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression ) ;
                    }
                }
            }

            fitness[i] = max_function ;
            if (fitness[i] > 1.0) {
                fitness[i] = 1.0 ;
            }
            if (fitness[i] < 0.0) {
                fitness[i] = 0.0 ;
            }
        }
    }
}

/// compute fitness
void compute_gaussian_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, double opt_fit, cmd_line &options ) {
    // gaussian with set mean and sd
    // float x = sum_all_tRNAs(tRNA.sequence * tRNA.expression)
    // fitness = ( 1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) ) ) * exp(-1 * ( (pow((x - fitness.mean), 2)) / (2 * pow(options.fitness_sd, 2)) )
    // ^ normalize this so that highest is 1.0 and all others are a fraction of that

    double temp_con = sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) ) ;
    for ( int i = 0 ; i < population.size() ; ++i ) {

        double x_ind = 0.0 ;
        // maternal x calculation block:
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            // if mutation destroys tRNA, get prop tissues NOT affected times normal genotype
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                x_ind += gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression ;
                            }
                            else {
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                // get prop tissues affected times effect
                                // plus also 1-tissues * normal genotype
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].maternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].maternal_trnas[g]).expression) ;
                                x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = ( temp_prop * ((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                        }
                        else {
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = ( temp_prop * ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                        }
                    }
                    // else would mean 10 muts so sequence value is 0, so add 0, so do nothing
                }
                else {
                    x_ind += ( ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression) ) ;
                }
                // for dups just add prop tissues affected times normal genotype
                // already added first gene above
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    x_ind += ( gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856))) ;
                }
            }
            // somatic deletion same as mutation destroying tRNA
            else {
                x_ind += ( gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression) ) ;
            }
        }

        // paternal x calculation block:
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            // if mutation destroys tRNA, get prop tissues NOT affected times normal genotype
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                x_ind += gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression ;
                            }
                            else {
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                // get prop tissues affected times effect
                                // plus also 1-tissues * normal genotype
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].paternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].paternal_trnas[g]).expression) ;
                                x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = ( temp_prop * ((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                        }
                        else {
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = ( temp_prop * ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + (1.0 - temp_prop * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                        }
                    }
                    // else would mean 10 muts so sequence value is 0, so add 0, so do nothing
                }
                else {
                    x_ind += ( ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
                }
                // for dups just add prop tissues affected times normal genotype
                // already added first gene above
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    x_ind += ( gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856))) ;
                }
            }
            // somatic deletion same as mutation destroying tRNA
            else {
                x_ind += ( gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
            }
        }

        fitness[i] = (1 / temp_con ) * ( exp(-1 * ( (pow((x_ind - options.fitness_mean), 2)) / (2 * (pow(options.fitness_sd, 2))) ) ) ) / opt_fit ;
        if (fitness[i] > 1.0) {
            fitness[i] = 1.0 ;
        }
        if (fitness[i] < 0.0) {
            fitness[i] = 0.0 ;
        }
    }
}

/// compute fitness
void compute_exp_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

    // gaussian with set mean and sd
    // float x = sum_all_tRNAs(tRNA.sequence * tRNA.expression)
    // fitness = ( 1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) ) ) * exp(-1 * ( (pow((x - fitness.mean), 2)) / (2 * pow(options.fitness_sd, 2)) )
    // ^ normalize this so that highest is 1.0 and all others are a fraction of that

    double temp_fitness [options.n] ;
    double max_fitness = 0.0 ;

    for ( int i = 0 ; i < population.size() ; ++i ) {

        double x_ind = 0.0 ;

        // maternal x calculation block:
        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].maternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            // mutation destroying tRNA = prop tissues NOT affected * normal genotype
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                x_ind += gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression ;
                            }
                            else {
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].maternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].maternal_trnas[g]).expression) ;
                                x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = (((mutations_to_function[((*population[i].maternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].maternal_trnas[g]).expression)) ;
                            x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                        }
                        else{
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = (temp_prop * ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) * ((*population[i].maternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression))) ;
                        }
                    }
                    // else would mean 10 muts so sequence value is 0, so add 0, so do nothing
                }
                else {
                    x_ind += ( ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression) ) ;
                }
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    x_ind += ( gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856)) ) ;
                }
            }
            else {
                x_ind += ( gsl_rng_uniform(rng) * ((*population[i].maternal_trnas[g]).sequence) * ((*population[i].maternal_trnas[g]).expression) ) ;
            }
        }

        // paternal x calculation block:
        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            if ( !gsl_ran_bernoulli( rng, options.somatic_del ) ){
                if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                    if ( (*population[i].paternal_trnas[g]).muts < options.max_mutations ) {
                        if ( options.dual_rates == true ){
                            // mutation destroying tRNA = prop tissues NOT affected * normal genotype
                            if ( gsl_ran_bernoulli( rng, options.prop_destroy )){
                                x_ind += gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression ;
                            }
                            else {
                                double temp_gamma = gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ;
                                if ( temp_gamma > 1.0 ){
                                    temp_gamma = 1.0 ;
                                }
                                double temp_prop = gsl_rng_uniform(rng) ;
                                double temp_function = (temp_prop * ((*population[i].paternal_trnas[g]).sequence - (1.0 - temp_gamma)) * (*population[i].paternal_trnas[g]).expression) ;
                                x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                            }
                        }
                        else if ( options.mutation_pathways == false ){
                            int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = (((mutations_to_function[((*population[i].paternal_trnas[g]).muts+1)])[random_index]) * ((*population[i].paternal_trnas[g]).expression)) ;
                            x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                        }
                        else{
                            int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                            double temp_prop = gsl_rng_uniform(rng) ;
                            double temp_function = (temp_prop * ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) * ((*population[i].paternal_trnas[g]).expression) ) ;
                            x_ind += (temp_function + ((1.0-temp_prop)*(((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression))) ;
                        }
                    }
                    // else would mean 10 muts so sequence value is 0, so add 0, so do nothing
                }
                else {
                    x_ind += ( ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
                }
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    x_ind += ( gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression + gsl_ran_gaussian(rng, 0.15856)) ) ;
                }
            }
            else {
                x_ind += ( gsl_rng_uniform(rng) * ((*population[i].paternal_trnas[g]).sequence) * ((*population[i].paternal_trnas[g]).expression) ) ;
            }
        }

        temp_fitness[i] = x_ind ;
        if ( x_ind > max_fitness ){
            max_fitness = x_ind ;
            //cout << max_fitness << endl ;
        }
    }

    for ( int i = 0 ; i < population.size() ; ++i ) {
        fitness[i] = 1.0 - exp(options.fitness_lambda * (temp_fitness[i] / max_fitness)) ;
        if (fitness[i] > 1.0) {
            fitness[i] = 1.0 ;
        }
        if (fitness[i] < 0.0) {
            fitness[i] = 0.0 ;
        }
    }
}


// function to call other models:
void compute_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, double opt_fit, cmd_line &options ) {
    if ( options.fitness_func == "redundant" ) {
        compute_redund_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
    }
    else if ( options.fitness_func == "model" ) {
        compute_model_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
    }
    else if ( ( options.fitness_func == "exp" ) or ( options.fitness_func == "exponential" ) ){
        compute_exp_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
    }
    else {
        compute_gaussian_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, opt_fit, options ) ;
    }
}



#endif