#ifndef __FITNESS_H
#define __FITNESS_H

/// compute fitness

// actual functions
void compute_redund_fitness( double fitness[], vector<individual> &population, std::map<int, vector<double>> &mutations_to_function, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

    // fitness = cumulative ( trna_seq * trna_exp ) / min_fitness

    for ( int i = 0 ; i < population.size() ; i ++ ) {
        double total_function = 0.0 ;

        //////////////////////////////
        /// MATERNAL FITNESS BLOCK ///
        //////////////////////////////

        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            double trna_function = 0.0 ;
            
            double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
            }

            double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
            }
            
            int num_SNPs = gsl_ran_poisson( rng, ((*population[i].maternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
            vector<double> props ; 
            if ( num_SNPs > 0 ){
                for ( int m = 0 ; m < num_SNPs ; m ++ ){
                    props.push_back( gsl_rng_uniform(rng) ) ;
                }
                std::sort( props.begin(), props.end() ) ;
                double prior_prop = 0.0 ;
                double current_seq = (*population[i].maternal_trnas[g]).sequence ;
                int counter = 0 ;
                for ( auto p : props ){
                    if ( counter == 0 ){
                        trna_function += ((p - prior_prop) * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                    }
                    else {
                        if ( (*population[i].maternal_trnas[g]).muts+counter < options.max_mutations ){
                            if ( options.dual_rates == true ){
                                current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)]).size() ;
                                current_seq -= (((mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)])[random_index])) ;
                            }
                            else {
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                                current_seq -= ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) ;
                            }
                            if ( current_seq < 0.0 ){
                                current_seq = 0.0 ;
                            }
                            trna_function += ((p - prior_prop) * current_seq * (*population[i].maternal_trnas[g]).expression) ;
                        }
                    }
                    prior_prop = p ;
                    counter ++ ;
                } 

            }

            else {
                trna_function = (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ;
            }
            trna_function += ( prop_dup * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
            trna_function -= ( prop_del * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
            total_function += trna_function ;
        }

        //////////////////////////////
        /// PATERNAL FITNESS BLOCK ///
        //////////////////////////////

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            double trna_function = 0.0 ;
            
            double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
            }

            double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
            }
            
            int num_SNPs = gsl_ran_poisson( rng, ((*population[i].paternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
            vector<double> props ; 
            if ( num_SNPs > 0 ){
                for ( int m = 0 ; m < num_SNPs ; m ++ ){
                    props.push_back( gsl_rng_uniform(rng) ) ;
                }
                std::sort( props.begin(), props.end() ) ;
                double prior_prop = 0.0 ;
                double current_seq = (*population[i].paternal_trnas[g]).sequence ;
                int counter = 0 ;
                for ( auto p : props ){
                    if ( counter == 0 ){
                        trna_function += ((p - prior_prop) * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                    }
                    else {
                        if ( (*population[i].paternal_trnas[g]).muts+counter < options.max_mutations ){
                            if ( options.dual_rates == true ){
                                current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)]).size() ;
                                current_seq -= (((mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)])[random_index])) ;
                            }
                            else {
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                                current_seq -= ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) ;
                            }
                            if ( current_seq < 0.0 ){
                                current_seq = 0.0 ;
                            }
                            trna_function += ((p - prior_prop) * current_seq * (*population[i].paternal_trnas[g]).expression) ;
                        }
                    }
                    prior_prop = p ;
                    counter ++ ;
                } 

            }
            else {
                trna_function = (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ;
            }

            trna_function += ( prop_dup * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
            trna_function -= ( prop_del * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
            total_function += trna_function ;
            
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
    //
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
            }
            else {
                fitness[i] = 1.0 - ( pow(options.model4_deverr, dad_function) ) ;
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

    else {

        // fitness = max ( trna_seq * trna_exp )

        for ( int i = 0 ; i < population.size() ; i ++ ) {
            double max_function = 0.0 ;

            //////////////////////////////
            /// MATERNAL FITNESS BLOCK ///
            //////////////////////////////

            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                double trna_function = 0.0 ;
                
                double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
                }

                double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
                if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                    double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
                }
                
                int num_SNPs = gsl_ran_poisson( rng, ((*population[i].maternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
                vector<double> props ; 
                if ( num_SNPs > 0 ){
                    for ( int m = 0 ; m < num_SNPs ; m ++ ){
                        props.push_back( gsl_rng_uniform(rng) ) ;
                    }
                    std::sort( props.begin(), props.end() ) ;
                    double prior_prop = 0.0 ;
                    double current_seq = (*population[i].maternal_trnas[g]).sequence ;
                    int counter = 0 ;
                    for ( auto p : props ){
                        if ( counter == 0 ){
                            trna_function += ((p - prior_prop) * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                        }
                        else {
                            if ( (*population[i].maternal_trnas[g]).muts+counter < options.max_mutations ){
                                if ( options.dual_rates == true ){
                                    current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                                }
                                else if ( options.mutation_pathways == false ){
                                    int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)]).size() ;
                                    current_seq -= (((mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)])[random_index])) ;
                                }
                                else {
                                    int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                                    current_seq -= ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) ;
                                }
                                if ( current_seq < 0.0 ){
                                    current_seq = 0.0 ;
                                }
                                trna_function += ((p - prior_prop) * current_seq * (*population[i].maternal_trnas[g]).expression) ;
                            }
                        }
                        prior_prop = p ;
                        counter ++ ;
                    } 
                }

                else {
                    trna_function = (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ;
                }

                trna_function += ( prop_dup * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
                trna_function -= ( prop_del * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
                if ( trna_function > max_function ){
                    max_function = trna_function ;
                }
            }

            //////////////////////////////
            /// PATERNAL FITNESS BLOCK ///
            //////////////////////////////

            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                double trna_function = 0.0 ;
                
                double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
                if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                    prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
                }

                double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
                if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                    double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
                }
                
                int num_SNPs = gsl_ran_poisson( rng, ((*population[i].paternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
                vector<double> props ; 
                if ( num_SNPs > 0 ){
                    for ( int m = 0 ; m < num_SNPs ; m ++ ){
                        props.push_back( gsl_rng_uniform(rng) ) ;
                    }
                    std::sort( props.begin(), props.end() ) ;
                    double prior_prop = 0.0 ;
                    double current_seq = (*population[i].paternal_trnas[g]).sequence ;
                    int counter = 0 ;
                    for ( auto p : props ){
                        if ( counter == 0 ){
                            trna_function += ((p - prior_prop) * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                        }
                        else {
                            if ( (*population[i].paternal_trnas[g]).muts+counter < options.max_mutations ){
                                if ( options.dual_rates == true ){
                                    current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                                }
                                else if ( options.mutation_pathways == false ){
                                    int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)]).size() ;
                                    current_seq -= (((mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)])[random_index])) ;
                                }
                                else {
                                    int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                                    current_seq -= ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) ;
                                }
                                if ( current_seq < 0.0 ){
                                    current_seq = 0.0 ;
                                }
                                trna_function += ((p - prior_prop) * current_seq * (*population[i].paternal_trnas[g]).expression) ;
                            }
                        }
                        prior_prop = p ;
                        counter ++ ;
                    } 
                }
                else {
                    trna_function = (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ;
                }

                trna_function += ( prop_dup * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
                trna_function -= ( prop_del * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
                if ( trna_function > max_function ){
                    max_function = trna_function ;
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

    // fitness = cumulative ( trna_seq * trna_exp ) / min_fitness

    for ( int i = 0 ; i < population.size() ; i ++ ) {
        double total_function = 0.0 ;

        //////////////////////////////
        /// MATERNAL FITNESS BLOCK ///
        //////////////////////////////

        for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
            double trna_function = 0.0 ;
            
            double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
            }

            double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
            }
            
            int num_SNPs = gsl_ran_poisson( rng, ((*population[i].maternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
            vector<double> props ; 
            if ( num_SNPs > 0 ){
                for ( int m = 0 ; m < num_SNPs ; m ++ ){
                    props.push_back( gsl_rng_uniform(rng) ) ;
                }
                std::sort( props.begin(), props.end() ) ;
                double prior_prop = 0.0 ;
                double current_seq = (*population[i].maternal_trnas[g]).sequence ;
                int counter = 0 ;
                for ( auto p : props ){
                    if ( counter == 0 ){
                        trna_function += ((p - prior_prop) * (((*population[i].maternal_trnas[g]).sequence) * (*population[i].maternal_trnas[g]).expression)) ;
                    }
                    else {
                        if ( (*population[i].maternal_trnas[g]).muts+counter < options.max_mutations ){
                            if ( options.dual_rates == true ){
                                current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)]).size() ;
                                current_seq -= (((mutations_to_function[((*population[i].maternal_trnas[g]).muts+counter)])[random_index])) ;
                            }
                            else {
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)]).size() ;
                                current_seq -= ((genotype_to_fitnesses[((*population[i].maternal_trnas[g]).genotype)])[random_index]) ;
                            }
                            if ( current_seq < 0.0 ){
                                current_seq = 0.0 ;
                            }
                            trna_function += ((p - prior_prop) * current_seq * (*population[i].maternal_trnas[g]).expression) ;
                        }
                    }
                    prior_prop = p ;
                    counter ++ ;
                } 

            }

            else {
                trna_function = (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ;
            }

            trna_function += ( prop_dup * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
            trna_function -= ( prop_del * (*population[i].maternal_trnas[g]).sequence * (*population[i].maternal_trnas[g]).expression ) ;
            total_function += trna_function ;
            
        }

        //////////////////////////////
        /// PATERNAL FITNESS BLOCK ///
        //////////////////////////////

        for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
            double trna_function = 0.0 ;
            
            double prop_dup = 0.0 ; // somatic dup as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_dup ) ){
                prop_dup = gsl_rng_uniform(rng) ; // will add that proportion times the sequence and expression of the given tRNA
            }

            double prop_del = 0.0 ; // somatic del as bernoulli draw then uniform for number of cells affected
            if ( gsl_ran_bernoulli( rng, options.somatic_del ) ){
                double prop_del = gsl_rng_uniform(rng) ; // will then subtract that proportion times the sequence and expression of the given tRNA
            }
            
            int num_SNPs = gsl_ran_poisson( rng, ((*population[i].paternal_trnas[g]).somatic) ) ; // get somatic SNPs as poisson process:
            vector<double> props ; 
            if ( num_SNPs > 0 ){
                for ( int m = 0 ; m < num_SNPs ; m ++ ){
                    props.push_back( gsl_rng_uniform(rng) ) ;
                }
                std::sort( props.begin(), props.end() ) ;
                double prior_prop = 0.0 ;
                double current_seq = (*population[i].paternal_trnas[g]).sequence ;
                int counter = 0 ;
                for ( auto p : props ){
                    if ( counter == 0 ){
                        trna_function += ((p - prior_prop) * (((*population[i].paternal_trnas[g]).sequence) * (*population[i].paternal_trnas[g]).expression)) ;
                    }
                    else {
                        if ( (*population[i].paternal_trnas[g]).muts+counter < options.max_mutations ){
                            if ( options.dual_rates == true ){
                                current_seq -= (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale ) ) ;
                            }
                            else if ( options.mutation_pathways == false ){
                                int random_index = rand() % (mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)]).size() ;
                                current_seq -= (((mutations_to_function[((*population[i].paternal_trnas[g]).muts+counter)])[random_index])) ;
                            }
                            else {
                                int random_index = rand() % (genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)]).size() ;
                                current_seq -= ((genotype_to_fitnesses[((*population[i].paternal_trnas[g]).genotype)])[random_index]) ;
                            }
                            if ( current_seq < 0.0 ){
                                current_seq = 0.0 ;
                            }
                            trna_function += ((p - prior_prop) * current_seq * (*population[i].paternal_trnas[g]).expression) ;
                        }
                    }
                    prior_prop = p ;
                    counter ++ ;
                } 

            }
            else {
                trna_function = (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ;
            }

            total_function += trna_function ;
            total_function += ( prop_dup * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
            total_function -= ( prop_del * (*population[i].paternal_trnas[g]).sequence * (*population[i].paternal_trnas[g]).expression ) ;
        }

        fitness[i] = (1 / temp_con ) * ( exp(-1 * ( (pow((total_function - options.fitness_mean), 2)) / (2 * (pow(options.fitness_sd, 2))) ) ) ) / opt_fit ;
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
    else {
        compute_gaussian_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, opt_fit, options ) ;
    }
}



#endif